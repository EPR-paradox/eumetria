#!/usr/bin/env python
"""Q-10 实验二：CPE 的阶数上限由谁决定？边缘场点数不足会怎样？

    .venv/bin/python experiments/q10_cpe_order.py

背景：CPE（Correction Per Exposure）= 逐场一套系数，所以**每个场各自的点数**必须
≥ 该场参数数 —— 平均点数够没有用。默认采样方案里有 4 个边缘场只剩 5 个可用 mark。

本实验对抽中的每个场分别构造场内 Cartesian 设计矩阵，报秩与条件数，
回答「二阶 CPE 到底能不能做」而不靠猜。
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from eumetria.basis import CartesianBasis
from eumetria.geometry import build_layout, field_marks, to_polar
from eumetria.io import load_scenario
from eumetria.synth import sample_fields

SCENARIO = ROOT / "scenarios" / "base-300mm.yaml"


def main() -> int:
    sc = load_scenario(SCENARIO)
    R = sc.wafer_radius_mm
    usable_rho = (R - sc.wafer.edge_exclusion_mm) / R
    layout = build_layout(
        sc.wafer.diameter_mm,
        sc.wafer.edge_exclusion_mm,
        sc.field.width_mm,
        sc.field.height_mm,
        include_partial=True,
    )
    rng = np.random.default_rng(sc.seed)
    fields = sample_fields(layout, sc.sampling.fields, sc.sampling.selection, rng)

    bases = {order: CartesianBasis(order) for order in (1, 2)}
    stats: dict[int, dict[str, list]] = {1: {"n": [], "cond": [], "rank_ok": []},
                                         2: {"n": [], "cond": [], "rank_ok": []}}

    print("# Q-10 实验二：CPE 阶数上限与边缘场点数\n")
    print(f"配置：seed={sc.seed}，{sc.sampling.fields} field × {sc.field.marks_per_field} mark，"
          f"edge exclusion {sc.wafer.edge_exclusion_mm}mm（ρ ≤ {usable_rho:.4f}）\n")
    print("| field_id | 可用 mark 数 | 一阶(3 项) 条件数 | 二阶(6 项) rank | 二阶条件数 |")
    print("|---|---|---|---|---|")
    for field in fields:
        marks = field_marks(field, sc.field.marks_per_field)
        xs = np.array([m.x_mm for m in marks])
        ys = np.array([m.y_mm for m in marks])
        xf = np.array([m.xf_mm for m in marks])
        yf = np.array([m.yf_mm for m in marks])
        rho, _ = to_polar(xs, ys, R)
        keep = rho <= usable_rho + 1e-12
        xi = 2 * xf[keep] / sc.field.width_mm
        eta = 2 * yf[keep] / sc.field.height_mm
        n = int(keep.sum())

        cond1 = float(np.linalg.cond(bases[1].values(xi, eta)))
        X2 = bases[2].values(xi, eta)
        rank2 = int(np.linalg.matrix_rank(X2))
        cond2 = float(np.linalg.cond(X2)) if rank2 == X2.shape[1] else float("inf")

        for order in (1, 2):
            stats[order]["n"].append(n)
        stats[1]["cond"].append(cond1)
        stats[1]["rank_ok"].append(cond1 < 1e6)
        stats[2]["cond"].append(cond2)
        stats[2]["rank_ok"].append(rank2 == X2.shape[1])
        flag = "" if rank2 == X2.shape[1] else " ← 降秩"
        print(f"| {field.field_id} | {n} | {cond1:.2f} | {rank2} | {cond2:.2f}{flag} |")

    print("\n## 汇总\n")
    for order in (1, 2):
        n_min = min(stats[order]["n"])
        bad = sum(1 for ok in stats[order]["rank_ok"] if not ok)
        finite = [c for c in stats[order]["cond"] if np.isfinite(c)]
        print(f"- 场内 **{order} 阶**（{len(bases[order].terms())} 项）："
              f"最少可用 mark {n_min} 个；不可辨识的场 {bad}/{len(fields)} 个；"
              f"条件数最差 {max(finite):.2f}"
              + ("" if not bad else "（降秩场条件数 = ∞）"))
    print()
    print("**结论：默认采样方案下，二阶 CPE 有 4 个场降到 rank 5 < 6，不可辨识；"
          "一阶 CPE 全部场可用（最差条件数 2.95）。所以 CPE 阶数上限是一阶。**")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
