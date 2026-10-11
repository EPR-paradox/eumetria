"""Q-10 系列实验共用的估计器定义（避免每个脚本各写一份、彼此漂移）。

模型（§3.1）：

    dx(x,y) = Σ_j a_j·Z_j(ρ,θ)  +  Σ_{p,q} b_pq·ξ^p η^q + ε
              └── wafer-level ──┘   └── field-level ──┘

场内基用 order=2 且**去掉常数项**（D-09）：terms = (1,0)(0,1)(2,0)(1,1)(0,2)。
常数项与 wafer 的 Z1 是同一列，归 wafer-level。

七个估计器覆盖"所有会把系数算出来的合理做法"，分四族：

  联合族      A  一次 lstsq 解 [Z|C]                          ← 项目采用（D-09）
              A2 backfitting 交替迭代到收敛
  聚合族      B  完整两级：逐场拟合 → 逐场量做均值/变差分解 → 全片拟合 → 场内高阶
              B' 退化的两级：逐场拟合后**只把平移**送去场间回归（丢掉线性项）
  池化顺序族  C  池化点 + 先场内后全片（一次扫描）
              C2 池化点 + 先全片后场内（一次扫描）
  逐场族      D  逐场各自拟合(二阶) → 对残差拟合全片
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

FieldTerms = list[str]


@dataclass(frozen=True)
class Setup:
    """一次实验的全部输入。"""

    y: np.ndarray               # (n,) 单通道测量（无噪声时为真值）
    zw: np.ndarray              # (n, n_wafer) wafer 设计矩阵
    cf: np.ndarray              # (n, n_field) 场内设计矩阵（不含常数项）
    cf1: np.ndarray             # (n, 3) 场内一阶 [1, ξ, η]（逐场拟合用，含常数项）
    cf_all: np.ndarray          # (n, 6) 场内全二阶（含常数项，逐场拟合用）
    zc: np.ndarray              # (n_fields, n_wafer) 场中心坐标上的 wafer 设计矩阵
    zdx: np.ndarray             # (n_fields, n_wafer) ∂Z_j/∂x(场中心)·(W/2)：全片梯度在单场内的一阶表现
    zdy: np.ndarray             # (n_fields, n_wafer) ∂Z_j/∂y(场中心)·(H/2)
    fid: np.ndarray             # (n,) field_id
    fids: list[int]             # 有观测的 field_id（升序）
    field_terms: FieldTerms     # 与 cf 的列一一对应


@dataclass(frozen=True)
class Result:
    name: str
    a: np.ndarray               # wafer 系数（按 zw 的列）
    b: np.ndarray               # 场内系数（按 field_terms）


def _lsq(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef


def est_joint(s: Setup) -> Result:
    """A：一次联合求解。"""
    n_w = s.zw.shape[1]
    sol = _lsq(np.hstack([s.zw, s.cf]), s.y)
    return Result("A 联合求解", sol[:n_w], sol[n_w:])


def est_backfitting(s: Setup, sweeps: int = 200) -> Result:
    """A2：交替迭代到收敛（线性模型下与 A 等价）。"""
    a = np.zeros(s.zw.shape[1])
    b = np.zeros(s.cf.shape[1])
    for _ in range(sweeps):
        a = _lsq(s.zw, s.y - s.cf @ b)
        b = _lsq(s.cf, s.y - s.zw @ a)
    return Result("A2 backfitting 迭代到收敛", a, b)


def est_lumped_full(s: Setup) -> Result:
    """B：完整聚合两级（业界叙述的完整版）。

    阶段 1  **逐场**拟合 [1, ξ, η] → 每场的平移 a0_k 与线性项 a1_k, a2_k
    阶段 2  把逐场量对**场中心坐标**回归，但关键是要认识到：

               a0_k ≈ W(场中心)                                  全片指纹在该场中心的值
               a1_k ≈ b_shared(1,0) + (∂W/∂x)(场中心)·(W_field/2)  场内线性指纹  +  全片梯度
               a2_k ≈ b_shared(0,1) + (∂W/∂y)(场中心)·(H_field/2)

            所以把三个量堆起来对 [Zc ; dZdx ; dZdy] 做一次 LS，就同时解出
            **全片低阶系数 a** 和**共享的场内线性指纹 b_lin** —— 这正是"均值/变差分解"：
            变差（随场中心变化的那部分）归全片，均值归场内。
    阶段 3  扣掉场内线性模型与全片低阶，余项里拟合共享的**场内高阶**项。

    这是"先场间（场中心坐标）再场内"的忠实实现 —— 注意它用的是**场中心坐标及其局部梯度**，
    不是把同一批点池化后按顺序扫描（那是 C，性质完全不同）。
    """
    per = {k: _lsq(s.cf1[s.fid == k], s.y[s.fid == k]) for k in s.fids}
    obs = np.array([per[k] for k in s.fids])          # (n_fields, 3)
    n_f, n_w = obs.shape[0], s.zc.shape[1]

    # 阶段 2 的块设计矩阵：(3·n_fields) × (n_w + 2)
    design = np.zeros((3 * n_f, n_w + 2))
    design[0:n_f, 0:n_w] = s.zc
    design[n_f:2 * n_f, 0:n_w] = s.zdx           # ξ 的系数来自全片梯度（内部已乘 W/2）
    design[2 * n_f:3 * n_f, 0:n_w] = s.zdy
    design[n_f:2 * n_f, n_w] = 1.0               # 共享场内线性指纹 (1,0)
    design[2 * n_f:3 * n_f, n_w + 1] = 1.0       # 共享场内线性指纹 (0,1)

    # obs 是 (n_fields, 3)，要按**列**展平才能对齐「先所有 a0、再所有 a1、再所有 a2」的块结构；
    # 默认的 C 序展平会把它交错成 (a0_0, a1_0, a2_0, a0_1, ...)，块就全错位了。
    sol = _lsq(design, obs.reshape(-1, order="F"))
    a = sol[:n_w]
    b_lin = {"1,0": sol[n_w], "0,1": sol[n_w + 1]}

    lin_cols = [i for i, t in enumerate(s.field_terms) if _deg(t) == 1]
    high_cols = [i for i, t in enumerate(s.field_terms) if _deg(t) >= 2]
    b = np.zeros(len(s.field_terms))
    for i in lin_cols:
        b[i] = b_lin[s.field_terms[i]]

    resid = s.y - s.zw @ a - s.cf[:, lin_cols] @ b[lin_cols]
    b[high_cols] = _lsq(s.cf[:, high_cols], resid)
    return Result("B 聚合两级（均值/变差分解）", a, b)


def est_lumped_translation_only(s: Setup) -> Result:
    """B'：退化的聚合两级 —— 阶段 1 的线性项丢掉，只用平移做场间回归。"""
    per = {k: _lsq(s.cf1[s.fid == k], s.y[s.fid == k]) for k in s.fids}
    trans = np.array([per[k][0] for k in s.fids])
    a = _lsq(s.zc, trans)
    resid = np.empty_like(s.y)
    for k in s.fids:
        m = s.fid == k
        resid[m] = s.y[m] - s.cf1[m] @ per[k]
    b = _lsq(s.cf, resid)
    return Result("B' 聚合两级（丢线性项）", a, b)


def est_pooled_field_first(s: Setup) -> Result:
    """C：池化点 + 先场内后全片（一次扫描）。"""
    b = _lsq(s.cf, s.y)
    a = _lsq(s.zw, s.y - s.cf @ b)
    return Result("C 池化顺序（先场内后全片）", a, b)


def est_pooled_wafer_first(s: Setup) -> Result:
    """C2：池化点 + 先全片后场内（一次扫描）。"""
    a = _lsq(s.zw, s.y)
    b = _lsq(s.cf, s.y - s.zw @ a)
    return Result("C2 池化顺序（先全片后场内）", a, b)


def est_per_field(s: Setup) -> Result:
    """D：逐场各自拟合二阶，再对残差拟合全片。"""
    resid = s.y.copy()
    for k in s.fids:
        m = s.fid == k
        resid[m] = s.y[m] - s.cf_all[m] @ _lsq(s.cf_all[m], s.y[m])
    a = _lsq(s.zw, resid)
    b = _lsq(s.cf, resid)
    return Result("D 逐场拟合 + 残差全片", a, b)


def _deg(term: str) -> int:
    p, q = term.split(",")
    return int(p) + int(q)


ALL_ESTIMATORS = (
    est_joint,
    est_backfitting,
    est_lumped_full,
    est_lumped_translation_only,
    est_pooled_field_first,
    est_pooled_wafer_first,
    est_per_field,
)
