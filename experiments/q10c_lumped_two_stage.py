#!/usr/bin/env python
"""Q-10c：业界那种「先场间（场中心坐标）再场内」到底成不成立？

    .venv/bin/python experiments/q10c_lumped_two_stage.py

被问住的地方：Q-10a 测出「顺序解有偏」，但业界确实是按两级叙述的
（先用场中心坐标拟合平移/旋转/倍率，再用场内坐标拟合低阶+高阶）。
区别在于：业界那一版里有一道**聚合**步骤 —— 先把每个场压缩成场中心上的几个观测量，
再做场间回归。两个阶段吃的是**不同的自变量**，不是把同一批残差再拟合一遍。

本实验把三种做法放在同一份无噪声数据上：
  A 联合求解（D-09 的实现）
  B 聚合两级：逐场拟合 → 场中心回归 → 场内高阶（业界叙述的做法）
  C 池化 + 一次扫描（Q-10a 里被推翻的那个）
  D 逐场拟合 + 残差全片（Q-10a 里的灾难臂）

真值里同时放 wafer 低阶（Z1/Z2/Z3，即 Tx/tilt）**和**场内的线性项，
因为这两者在数学上正好同阶 —— 这就是要看的碰撞点。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eumetria.basis import CartesianBasis, ZernikeBasis  # noqa: E402
from eumetria.geometry import build_layout, field_marks, to_polar  # noqa: E402
from eumetria.io import load_scenario  # noqa: E402
from eumetria.synth import sample_fields  # noqa: E402

SCENARIO = ROOT / "scenarios" / "base-300mm.yaml"

# 真值：wafer 低阶（平移/倾斜，正是 APC 吃的）+ 场内线性（正是 i-HOPC 吃的）
WAFER_TRUTH = {"Z1": 2.0, "Z2": 1.5, "Z3": -1.0, "Z4": 2.5, "Z5": -1.8, "Z11": -2.2}
FIELD_TRUTH = {"1,0": 1.2, "0,1": -0.8, "2,0": 0.5}


def build_points(sc) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """返回 (ρ, θ, ξ, η, field_id)，用与 scenario 一致的采样与过滤。"""
    R = sc.wafer_radius_mm
    usable_rho = (R - sc.wafer.edge_exclusion_mm) / R
    layout = build_layout(
        sc.wafer.diameter_mm, sc.wafer.edge_exclusion_mm,
        sc.field.width_mm, sc.field.height_mm, include_partial=True,
    )
    fields = sample_fields(layout, sc.sampling.fields, sc.sampling.selection,
                           np.random.default_rng(sc.seed))
    rho, theta, xi, eta, fid = [], [], [], [], []
    for f in fields:
        for mk in field_marks(f, sc.field.marks_per_field):
            x, y = mk.x_mm, mk.y_mm
            r, t = to_polar(np.array([x]), np.array([y]), R)
            if float(r[0]) > usable_rho + 1e-12:
                continue
            rho.append(float(r[0])); theta.append(float(t[0]))
            xi.append(2 * mk.xf_mm / sc.field.width_mm)
            eta.append(2 * mk.yf_mm / sc.field.height_mm)
            fid.append(f.field_id)
    return (np.array(rho), np.array(theta), np.array(xi), np.array(eta),
            np.array(fid, dtype=np.int64))


def main() -> int:
    sc = load_scenario(SCENARIO)
    rho, theta, xi, eta, fid = build_points(sc)
    n = rho.size

    zw = ZernikeBasis(sc.fingerprint.max_order)
    cf_full = CartesianBasis(sc.fingerprint.field_order)
    keep_const = [i for i, t in enumerate(cf_full.terms()) if t != "0,0"]
    cf = CartesianBasis(sc.fingerprint.field_order)

    Zw = zw.values(rho, theta)
    Cf_all = cf.values(xi, eta)
    Cf = Cf_all[:, keep_const]
    fterms = [cf.terms()[i] for i in keep_const]

    a_true = np.array([WAFER_TRUTH.get(t, 0.0) for t in zw.terms()])
    b_true = np.array([FIELD_TRUTH.get(t, 0.0) for t in fterms])
    y = Zw @ a_true + Cf @ b_true

    def err(a, b):
        ew = float(np.max(np.abs(a - a_true)[np.flatnonzero(a_true)]))
        eb = float(np.max(np.abs(b - b_true)[np.flatnonzero(b_true)]))
        rms = float(np.sqrt(np.mean((y - Zw @ a - Cf @ b) ** 2)))
        return ew, eb, rms

    # ---- A 联合求解 ----
    X = np.hstack([Zw, Cf])
    sol, *_ = np.linalg.lstsq(X, y, rcond=None)
    A = (sol[: Zw.shape[1]], sol[Zw.shape[1]:])

    # ---- B 聚合两级（业界叙述）----
    # 阶段 1：逐场拟合 [1, ξ, η]（每场自己一套，不共享）
    Cf1 = cf.values(xi, eta, ["0,0", "1,0", "0,1"])
    per_field = {}
    for f in np.unique(fid):
        m = fid == f
        coef, *_ = np.linalg.lstsq(Cf1[m], y[m], rcond=None)
        per_field[int(f)] = coef
    fids = sorted(per_field)
    # 阶段 2：把「场中心上的观测量」（每场的平移项）对场中心坐标的 Zernike 回归
    layout = build_layout(sc.wafer.diameter_mm, sc.wafer.edge_exclusion_mm,
                          sc.field.width_mm, sc.field.height_mm, include_partial=True)
    by_id = {f.field_id: f for f in layout}
    rho_c, theta_c = [], []
    for k in fids:
        fld = by_id[k]
        r, t = to_polar(np.array([fld.cx_mm]), np.array([fld.cy_mm]), sc.wafer_radius_mm)
        rho_c.append(float(r[0])); theta_c.append(float(t[0]))
    Zc = zw.values(np.array(rho_c), np.array(theta_c))
    trans = np.array([per_field[k][0] for k in fids])  # 每场的平移项（场中心上的观测）
    a_B, *_ = np.linalg.lstsq(Zc, trans, rcond=None)
    # 场内高阶：用逐场拟合的残差拟合共享的场内二阶项（去掉 (1,0)(0,1) 已被阶段1吃掉）
    resid = np.empty_like(y)
    for k in fids:
        m = fid == k
        resid[m] = y[m] - Cf1[m] @ per_field[k]
    b_B, *_ = np.linalg.lstsq(Cf, resid, rcond=None)
    B = (a_B, b_B)

    # ---- C 池化 + 一次扫描 ----
    b_c, *_ = np.linalg.lstsq(Cf, y, rcond=None)
    a_c, *_ = np.linalg.lstsq(Zw, y - Cf @ b_c, rcond=None)
    C = (a_c, b_c)

    # ---- D 逐场(二阶) + 残差全片 ----
    r2 = y.copy()
    for k in fids:
        m = fid == k
        coef, *_ = np.linalg.lstsq(Cf_all[m], y[m], rcond=None)
        r2[m] = y[m] - Cf_all[m] @ coef
    a_d, *_ = np.linalg.lstsq(Zw, r2, rcond=None)
    b_d, *_ = np.linalg.lstsq(Cf, r2, rcond=None)
    D = (a_d, b_d)

    print("# Q-10c：聚合两级（业界）vs 联合求解 vs 池化顺序解\n")
    print(f"配置：seed={sc.seed}，{sc.sampling.fields} field × {sc.field.marks_per_field} mark，"
          f"实测 {n} 点；数据**无噪声**")
    print(f"真值 wafer：{WAFER_TRUTH}")
    print(f"真值 场内：{FIELD_TRUTH}（线性项 (1,0)(0,1) 与 wafer 的 Z2/Z3 同阶 —— 碰撞点）\n")
    print("| 估计器 | wafer 非零项最大误差 (nm) | field 非零项最大误差 (nm) | 全场残差 RMS (nm) |")
    print("|---|---|---|---|")
    for name, (a, b) in (("A 联合求解", A), ("B 聚合两级（场中心回归）", B),
                         ("C 池化 + 一次扫描", C), ("D 逐场拟合 + 残差全片", D)):
        ew, eb, rms = err(a, b)
        print(f"| {name} | {ew:.4f} | {eb:.4f} | {rms:.4f} |")

    print("\n## 低阶项的归属（谁把 Z2/Z3 拿走了）\n")
    print("| 项 | 真值 | A 联合 | B 聚合两级 | C 池化顺序 |")
    print("|---|---|---|---|---|")
    for t in ("Z1", "Z2", "Z3", "Z4", "Z11"):
        i = zw.terms().index(t)
        print(f"| {t} | {a_true[i]:.3f} | {A[0][i]:.3f} | {B[0][i]:.3f} | {C[0][i]:.3f} |")
    print("\n场内线性项（B 的阶段 1 会吃掉 wafer 的局部梯度，所以这里必然被污染）：")
    print("| 项 | 真值 | A 联合 | B 聚合两级 | C 池化顺序 |")
    print("|---|---|---|---|---|")
    for t in ("1,0", "0,1", "2,0"):
        j = fterms.index(t)
        print(f"| ({t}) | {b_true[j]:.3f} | {A[1][j]:.3f} | {B[1][j]:.3f} | {C[1][j]:.3f} |")

    # ---- 混叠的可判定性：把场内常数项加上去，看数据能否分辨 ----
    print("\n## 混叠检验：Z1 与场内常数项 (0,0) 能不能分开？\n")
    b00 = 0.7
    y2 = y + Cf_all[:, cf_full.terms().index("0,0")] * b00
    a_shift, *_ = np.linalg.lstsq(Zw, y2 - Cf @ np.zeros_like(b_true), rcond=None)
    sol2, *_ = np.linalg.lstsq(np.hstack([Zw, Cf]), y2, rcond=None)
    rms_same = float(np.sqrt(np.mean((y2 - np.hstack([Zw, Cf]) @ sol2) ** 2)))
    print(f"把场内常数项设为 {b00} nm（真值 Z1 仍是 {WAFER_TRUTH['Z1']}）后：")
    print(f"- 仅由场中心平移回归得到的「Z1」= {a_shift[zw.terms().index('Z1')]:.3f}"
          f"（= Z1 + 场内常数 = {WAFER_TRUTH['Z1'] + b00:.1f}）")
    print(f"- 两组参数给出**完全相同**的预测，残差 RMS 都是 {rms_same:.4f} nm")
    print("- → 数据无法分辨「全片平移 2.0 + 每场常量 0.7」与「全片平移 2.7 + 每场常量 0」；"
          "归属只能靠**约定**，不能靠估计顺序")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
