#!/usr/bin/env python
"""Q-10 实验一：两级拟合的「顺序」有影响吗？应该用哪种估计器？

    .venv/bin/python experiments/q10_fit_order.py

背景：详设 §3.1 把模型写成

    dx(x,y) = Σ_j a_j·Z_j(ρ,θ)  +  Σ_{p,q} b_pq·ξ^p η^q  +  ε
              └── wafer-level ──┘   └── field-level ──┘

初版 §7.4 规定「先拟合 field-level，再用残差拟合 wafer-level，顺序不可颠倒」。
这个实验用已知真值检验这句话：把五种估计器放在同一份**无噪声**数据上，
直接量出各自的系数偏差。无噪声是关键 —— 有噪声时偏差会被噪声掩盖。

结论（2026-10，默认场景，seed 20261003）：
  1. 顺序解有偏，且偏差方向取决于先后；「顺序不可颠倒」不能成立
  2. 联合求解与 backfitting 迭代到收敛都精确到 0.0000 nm —— 用哪个都行，别停在一次扫描
  3. 逐场各自拟合再对残差拟合是全片信号被按场吸收，误差 2.7 nm，不可用
  4. 顺带查出模型缺陷：Zernike 的 Z1 与 Cartesian 的 (0,0) 是同一列（常量重复）
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eumetria.basis import CartesianBasis, ZernikeBasis
from eumetria.io import load_scenario
from eumetria.synth import generate

SCENARIO = ROOT / "scenarios" / "base-300mm.yaml"


def coupling(zw: np.ndarray, cf: np.ndarray, zterms, fterms, top: int = 3):
    """两个基块的列间最大 |cos|（列先中心化）—— 衡量「不正交到什么程度」。"""
    pairs = []
    for i, zt in enumerate(zterms):
        zi = zw[:, i] - zw[:, i].mean()
        nz = np.linalg.norm(zi)
        for j, ft in enumerate(fterms):
            cj = cf[:, j] - cf[:, j].mean()
            nc = np.linalg.norm(cj)
            if nz > 0 and nc > 0:
                pairs.append((abs(float(zi @ cj) / (nz * nc)), zt, ft))
    pairs.sort(reverse=True)
    return pairs[:top]


def identical_columns(zw: np.ndarray, cf: np.ndarray, zterms, fterms):
    """找出两个基块里完全相同的列 —— 常量项重复会让设计矩阵降秩且条件数爆炸。"""
    hits = []
    for i, zt in enumerate(zterms):
        for j, ft in enumerate(fterms):
            if np.allclose(zw[:, i], cf[:, j], rtol=0, atol=1e-12):
                hits.append((zt, ft))
    return hits


def main() -> int:
    sc = load_scenario(SCENARIO)
    run = generate(sc)
    pts = run.points
    y = run.truth_values[:, 0]  # 用 x 通道、无噪声

    zw_basis = ZernikeBasis(sc.fingerprint.max_order)
    cf_basis = CartesianBasis(sc.fingerprint.field_order)
    xi, eta = pts.xi_eta()
    Zw = zw_basis.values(pts.rho, pts.theta_rad)
    Cf = cf_basis.values(xi, eta)

    a_true = np.array([sc.fingerprint.terms.get(t, 0.0) for t in zw_basis.terms()])
    b_true = np.array([sc.fingerprint.field_terms.get(t, 0.0) for t in cf_basis.terms()])
    n_w, n_f = Zw.shape[1], Cf.shape[1]
    fid = pts.field_id

    def solve_joint(zw, cf):
        X = np.hstack([zw, cf])
        sol, *_ = np.linalg.lstsq(X, y, rcond=None)
        return sol[: zw.shape[1]], sol[zw.shape[1]:]

    def solve_field_first(zw, cf):
        b, *_ = np.linalg.lstsq(cf, y, rcond=None)
        a, *_ = np.linalg.lstsq(zw, y - cf @ b, rcond=None)
        return a, b

    def solve_wafer_first(zw, cf):
        a, *_ = np.linalg.lstsq(zw, y, rcond=None)
        b, *_ = np.linalg.lstsq(cf, y - zw @ a, rcond=None)
        return a, b

    def solve_backfitting(zw, cf, sweeps: int = 200):
        a = np.zeros(zw.shape[1])
        b = np.zeros(cf.shape[1])
        for _ in range(sweeps):
            a, *_ = np.linalg.lstsq(zw, y - cf @ b, rcond=None)
            b, *_ = np.linalg.lstsq(cf, y - zw @ a, rcond=None)
        return a, b

    def solve_per_field(zw, cf):
        resid = y.copy()
        for f in np.unique(fid):
            m = fid == f
            if int(m.sum()) >= cf.shape[1]:
                coef, *_ = np.linalg.lstsq(cf[m], y[m], rcond=None)
                resid[m] = y[m] - cf[m] @ coef
        a, *_ = np.linalg.lstsq(zw, resid, rcond=None)
        b, *_ = np.linalg.lstsq(cf, resid, rcond=None)
        return a, b

    def score(a, b, zw, cf):
        nz_w, nz_b = np.flatnonzero(a_true), np.flatnonzero(b_true)
        err_w = float(np.max(np.abs(a - a_true)[nz_w]))
        err_b = float(np.max(np.abs(b - b_true)[nz_b]))
        rms = float(np.sqrt(np.mean((y - zw @ a - cf @ b) ** 2)))
        return err_w, err_b, rms

    print("# Q-10 实验一：两级拟合的顺序 / 估计器选择\n")
    print(f"配置：seed={sc.seed}，{sc.wafer.diameter_mm:.0f}mm，"
          f"field {sc.field.width_mm}×{sc.field.height_mm}mm，"
          f"{sc.sampling.fields} field × {sc.field.marks_per_field} mark，"
          f"实测 {len(pts)} 点；数据**无噪声**（真值本身）")
    print(f"wafer 基：Zernike Noll Z1..Z{sc.fingerprint.max_order}（{n_w} 项）；"
          f"field 基：Cartesian order={sc.fingerprint.field_order}（{n_f} 项）\n")

    rows = [
        ("联合求解（一次 lstsq）", solve_joint),
        ("backfitting 迭代到收敛", solve_backfitting),
        ("先场内后全片（一次扫描）", solve_field_first),
        ("先全片后场内（一次扫描）", solve_wafer_first),
        ("逐场各自拟合(二阶)后再拟合全片", solve_per_field),
    ]
    print("| 估计器 | wafer 非零项最大误差 (nm) | field 非零项最大误差 (nm) | 全场残差 RMS (nm) |")
    print("|---|---|---|---|")
    results = {}
    for name, fn in rows:
        a, b = fn(Zw, Cf)
        ew, eb, rms = score(a, b, Zw, Cf)
        results[name] = (a, b, ew, eb, rms)
        print(f"| {name} | {ew:.4f} | {eb:.4f} | {rms:.4f} |")

    print("\n## 基块耦合强度（列中心化后的最大 |cos|）\n")
    print("| 列对 | \\|cos\\| |")
    print("|---|---|")
    for c, zt, ft in coupling(Zw, Cf, zw_basis.terms(), cf_basis.terms()):
        print(f"| Zernike {zt} vs Cartesian ({ft}) | {c:.3f} |")

    dup = identical_columns(Zw, Cf, zw_basis.terms(), cf_basis.terms())
    print("\n## 模型缺陷：完全相同的列（常量项重复）\n")
    if dup:
        for zt, ft in dup:
            print(f"- Zernike **{zt}** 与 Cartesian **({ft})** 是同一列（都恒等于 1）")
    print(f"\n联合设计矩阵条件数：**{np.linalg.cond(np.hstack([Zw, Cf])):.3e}**"
          f"（Zw 单独 {np.linalg.cond(Zw):.1f}，Cf 单独 {np.linalg.cond(Cf):.1f}）")
    print("→ 常量项归属不清，联合模型降秩。取「场内基不含常数项」后：")

    cf_no_const = CartesianBasis(sc.fingerprint.field_order)
    keep = [i for i, t in enumerate(cf_no_const.terms()) if t != "0,0"]
    Cf2 = Cf[:, keep]
    a2, _b2 = solve_joint(Zw, Cf2)
    _ew2, _eb2, _rms2 = score(a2, np.zeros_like(b_true), Zw, Cf)
    print(f"   cond = {np.linalg.cond(np.hstack([Zw, Cf2])):.1f}；"
          f"联合解 wafer 非零项误差 {float(np.max(np.abs(a2 - a_true)[np.flatnonzero(a_true)])):.4f} nm")

    print("\n## 联合解的 wafer 系数（与真值对照）\n")
    print("| 项 | 真值 | 联合解 | 先场内后全片 | 逐场后全片 |")
    print("|---|---|---|---|---|")
    for i, t in enumerate(zw_basis.terms()):
        print(f"| {t} | {a_true[i]:.3f} | {results['联合求解（一次 lstsq）'][0][i]:.3f} | "
              f"{results['先场内后全片（一次扫描）'][0][i]:.3f} | "
              f"{results['逐场各自拟合(二阶)后再拟合全片'][0][i]:.3f} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
