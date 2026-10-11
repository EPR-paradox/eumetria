#!/usr/bin/env python
"""Q-10 汇总：七种估计器的效果一览。

    .venv/bin/python experiments/q10_estimator_comparison.py

这是 Q-10 系列的**汇总表**：估计器定义集中在 `_estimators.py`，本脚本把它们跑在
**两个不同真值**上并排给出 —— 因为同一种做法的表现会随真值组成而变（Q-10c 的教训：
真值里不含低阶项时，"池化顺序解"看着只有 0.03nm 偏差；含低阶项时它直接崩）。

  T1 无低阶：wafer 只有 Z4/Z5/Z6/Z7/Z11（Z1=Z2=Z3=0）
  T2 含低阶：再加 Z1=2.0、Z2=1.5、Z3=-1.0（全片平移与倾斜，正是 APC 吃的）

误差按组拆开报，而不是只报一个"最大误差" —— 低阶和高阶的失败模式完全不同：
低阶项争的是"归属"（混叠），高阶项争的是"统计效率"（观测数）。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments"))

from _estimators import ALL_ESTIMATORS, Result, Setup  # noqa: E402
from eumetria.basis import CartesianBasis, ZernikeBasis  # noqa: E402
from eumetria.geometry import build_layout, field_marks, to_polar  # noqa: E402
from eumetria.io import load_scenario  # noqa: E402
from eumetria.synth import sample_fields  # noqa: E402

SCENARIO = ROOT / "scenarios" / "base-300mm.yaml"

FIELD_TRUTH = {"1,0": 1.2, "0,1": -0.8, "2,0": 0.5}
HIGH_WAFER = {"Z4": 2.5, "Z5": -1.8, "Z6": 1.1, "Z7": 0.9, "Z11": -2.2}
LOW_WAFER = {"Z1": 2.0, "Z2": 1.5, "Z3": -1.0}

LOW_TERMS = ("Z1", "Z2", "Z3")
HIGH_TERMS = ("Z4", "Z5", "Z6", "Z7", "Z11")


def build_setup(sc, wafer_truth: dict[str, float]) -> tuple[Setup, np.ndarray, np.ndarray]:
    R = sc.wafer_radius_mm
    usable_rho = (R - sc.wafer.edge_exclusion_mm) / R
    layout = build_layout(sc.wafer.diameter_mm, sc.wafer.edge_exclusion_mm,
                          sc.field.width_mm, sc.field.height_mm, include_partial=True)
    fields = sample_fields(layout, sc.sampling.fields, sc.sampling.selection,
                           np.random.default_rng(sc.seed))
    rho, theta, xi, eta, fid = [], [], [], [], []
    for f in fields:
        for mk in field_marks(f, sc.field.marks_per_field):
            r, t = to_polar(np.array([mk.x_mm]), np.array([mk.y_mm]), R)
            if float(r[0]) > usable_rho + 1e-12:
                continue
            rho.append(float(r[0])); theta.append(float(t[0]))
            xi.append(2 * mk.xf_mm / sc.field.width_mm)
            eta.append(2 * mk.yf_mm / sc.field.height_mm)
            fid.append(f.field_id)
    rho = np.array(rho); theta = np.array(theta)
    xi = np.array(xi); eta = np.array(eta); fid = np.array(fid, dtype=np.int64)

    zw_b = ZernikeBasis(sc.fingerprint.max_order)
    cf_b = CartesianBasis(sc.fingerprint.field_order)
    keep = [i for i, t in enumerate(cf_b.terms()) if t != "0,0"]
    fterms = [cf_b.terms()[i] for i in keep]

    Zw = zw_b.values(rho, theta)
    Cf = cf_b.values(xi, eta)[:, keep]
    Cf1 = cf_b.values(xi, eta, ["0,0", "1,0", "0,1"])
    Cf_all = cf_b.values(xi, eta)

    fids = sorted({int(k) for k in fid})
    by_id = {f.field_id: f for f in layout}
    cx = np.array([by_id[k].cx_mm for k in fids])
    cy = np.array([by_id[k].cy_mm for k in fids])
    rho_c, theta_c = [], []
    for k in fids:
        r, t = to_polar(np.array([by_id[k].cx_mm]), np.array([by_id[k].cy_mm]), R)
        rho_c.append(float(r[0])); theta_c.append(float(t[0]))
    Zc = zw_b.values(np.array(rho_c), np.array(theta_c))
    # 全片指纹在单场内的局部梯度：∂Z_j/∂x(场中心)·(W/2)（ξ 是半宽归一化坐标）
    zdx = _grad(zw_b, cx, cy, R, sc.field.width_mm / 2, axis=0)
    zdy = _grad(zw_b, cx, cy, R, sc.field.height_mm / 2, axis=1)

    a_true = np.array([wafer_truth.get(t, 0.0) for t in zw_b.terms()])
    b_true = np.array([FIELD_TRUTH.get(t, 0.0) for t in fterms])
    y = Zw @ a_true + Cf @ b_true

    setup = Setup(y=y, zw=Zw, cf=Cf, cf1=Cf1, cf_all=Cf_all, zc=Zc,
                  zdx=zdx, zdy=zdy, fid=fid, fids=fids, field_terms=fterms)
    return setup, a_true, b_true


def _grad(zw_b: ZernikeBasis, cx, cy, R: float, scale: float, axis: int) -> np.ndarray:
    """∂Z_j/∂x（或 ∂y）在场中心的数值中心差分，再乘该方向的半宽。"""
    step = 1.0  # mm
    out = np.zeros((len(cx), len(zw_b.terms())))
    for i, (x, y) in enumerate(zip(cx, cy)):
        opts = (np.array([x + step]), np.array([y])) if axis == 0 else (
            np.array([x]), np.array([y + step]))
        omts = (np.array([x - step]), np.array([y])) if axis == 0 else (
            np.array([x]), np.array([y - step]))
        rp, tp = to_polar(*opts, R)
        rm, tm = to_polar(*omts, R)
        out[i] = (zw_b.values(rp, tp)[0] - zw_b.values(rm, tm)[0]) / (2 * step) * scale
    return out


def errors(res: Result, setup: Setup, a_true, b_true) -> dict[str, float]:
    a_err = np.abs(res.a - a_true)
    b_err = np.abs(res.b - b_true)
    zt = [f"Z{i + 1}" for i in range(len(a_true))]
    low_w = max(a_err[zt.index(t)] for t in LOW_TERMS)
    high_w = max(a_err[zt.index(t)] for t in HIGH_TERMS)
    low_f = max(b_err[setup.field_terms.index(t)] for t in ("1,0", "0,1"))
    high_f = b_err[setup.field_terms.index("2,0")]
    rms = float(np.sqrt(np.mean((setup.y - setup.zw @ res.a - setup.cf @ res.b) ** 2)))
    return {"低阶 wafer": low_w, "高阶 wafer": high_w,
            "低阶场内": low_f, "高阶场内": high_f, "残差 RMS": rms}


def stage2_diagnostics(setup: Setup) -> None:
    """诊断聚合两级的阶段 2 为什么不是零误差：块设计矩阵的病态来源。

    阶段 2 把「每场平移 a0_k、场内线性 a1_k/a2_k」对 [Zc ; dZdx ; dZdy] 回归。
    对 Z2/Z3 这类**梯度在全片上恒定**的项，dZdx/dZdy 的列是常数列 ——
    与代表「共享场内线性指纹」的全 1 列**完全共线**，只能靠平移块上的空间变化区分，
    所以低阶项是弱确定的。
    """
    n_f, n_w = setup.zc.shape[0], setup.zc.shape[1]
    design = np.zeros((3 * n_f, n_w + 2))
    design[0:n_f, 0:n_w] = setup.zc
    design[n_f:2 * n_f, 0:n_w] = setup.zdx
    design[2 * n_f:3 * n_f, 0:n_w] = setup.zdy
    design[n_f:2 * n_f, n_w] = 1.0
    design[2 * n_f:3 * n_f, n_w + 1] = 1.0

    labels = [f"Z{j + 1}" for j in range(n_w)] + ["b(1,0)", "b(0,1)"]
    centered = design - design.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centered, axis=0)
    pairs = []
    for i in range(design.shape[1]):
        for j in range(i + 1, design.shape[1]):
            if norms[i] > 0 and norms[j] > 0:
                c = abs(float(centered[:, i] @ centered[:, j]) / (norms[i] * norms[j]))
                pairs.append((c, labels[i], labels[j]))
    pairs.sort(reverse=True)
    print(f"  阶段 2 设计矩阵 cond = {np.linalg.cond(design):.3e}；最强共线列对：", end="")
    print("、".join(f"{a}~{b} |cos|={c:.3f}" for c, a, b in pairs[:3]))


def main() -> int:
    sc = load_scenario(SCENARIO)
    print("# Q-10 汇总：七种估计器的效果一览\n")
    print("配置：seed=%d，%d field × %d mark，实测 188 点；数据**无噪声**；"
          "场内基 order=2 去常数项（D-09）\n" % (sc.seed, sc.sampling.fields, sc.field.marks_per_field))

    # T0 是对照：如果全片指纹在单个场里**精确**是线性的，聚合两级应当没有截断误差
    for label, wafer_truth in (("T0 纯低阶（全片指纹在场内精确线性）", LOW_WAFER),
                              ("T1 无低阶", HIGH_WAFER),
                              ("T2 含低阶", {**HIGH_WAFER, **LOW_WAFER})):
        setup, a_true, b_true = build_setup(sc, wafer_truth)
        print(f"## {label}（wafer 真值 {wafer_truth}）\n")
        print("| 估计器 | 低阶 wafer (nm) | 高阶 wafer (nm) | 低阶场内 (nm) | 高阶场内 (nm) | 残差 RMS (nm) |")
        print("|---|---|---|---|---|---|")
        for fn in ALL_ESTIMATORS:
            res = fn(setup)
            e = errors(res, setup, a_true, b_true)
            print(f"| {res.name} | {e['低阶 wafer']:.4f} | {e['高阶 wafer']:.4f} | "
                  f"{e['低阶场内']:.4f} | {e['高阶场内']:.4f} | {e['残差 RMS']:.4f} |")
        stage2_diagnostics(setup)
        print()
    print("读法：低阶项看「归属」（能不能从混叠里分出来），高阶项看「统计效率」（观测数够不够）。")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
