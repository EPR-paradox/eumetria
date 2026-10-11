"""流水线编排：synth → fit → score → report。

**这一层是"怎么串"，不是"怎么算"，也不是"怎么解析参数"。**
- 算法在各自模块里（synth / fit / metrics / report）
- 参数解析与退出码在 cli.py
- 零参数入口在仓库根的 main.py

分层的原因是实打实的：编排原本埋在 `cli._cmd_run` 里，于是 main.py、notebook、
集成测试想跑一遍完整流水线都必须先构造一个 argparse.Namespace —— 前端一多就会各写一份。
现在三个前端都调 `run_all()`，编排只有一份。
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .fit import FitConfig, fit_field, write_fit_result
from .io import load_measurements, load_scenario, load_truth
from .io.schemas import Scenario
from .io.writers import write_measurements, write_truth
from .metrics import compare
from .synth import generate

MEASUREMENTS_NAME = "measurements.csv"
TRUTH_NAME = "truth.json"
FIT_RESULT_NAME = "fit_result.json"
METRICS_NAME = "metrics.json"
REPORT_DIRNAME = "report"


def scenario_with_seed(path: str | Path, seed: int | None = None) -> Scenario:
    scenario = load_scenario(path)
    if seed is not None:
        scenario = scenario.model_copy(update={"seed": int(seed)})
    return scenario


def read_measurements(scenario: Scenario, path: str | Path):
    """按 §6.2.5 的要求读 measurements（重算 ρ/θ 需要硅片半径）。"""
    return load_measurements(
        path,
        wafer_radius_mm=scenario.wafer_radius_mm,
        field_width_mm=scenario.field.width_mm,
        field_height_mm=scenario.field.height_mm,
    )


def fit_config_for(scenario: Scenario) -> FitConfig:
    return FitConfig(
        basis=scenario.fingerprint.basis,
        max_order=scenario.fingerprint.max_order,
        field_order=scenario.fingerprint.field_order,
    )


# ------------------------------------------------------------------ 四步


def synth(scenario: Scenario, out_dir: str | Path) -> dict[str, Any]:
    run = generate(scenario)
    out = Path(out_dir)
    ms_path = write_measurements(out / MEASUREMENTS_NAME, run.measurements)
    truth_path = write_truth(out / TRUTH_NAME, run.truth)
    return {
        "command": "synth",
        "seed": scenario.seed,
        "n_points": len(run.measurements),
        "n_fields": len({int(v) for v in run.measurements.field_id}),
        "measurements": str(ms_path),
        "truth": str(truth_path),
        "tis_estimate_nm": [round(v, 4) for v in run.tis_estimate_nm],
    }


def fit(
    scenario: Scenario, measurements_path: str | Path, out_dir: str | Path
) -> dict[str, Any]:
    ms = read_measurements(scenario, measurements_path)
    result = fit_field(ms, fit_config_for(scenario))
    path = write_fit_result(Path(out_dir) / FIT_RESULT_NAME, result)
    return {
        "command": "fit",
        "n_points": result.n_points,
        "n_params": result.n_params,
        "condition_number": round(result.condition_number, 4),
        "residual_rms_nm": round(result.residual_rms_nm, 6),
        "cpe_order_ceiling": result.identifiability.cpe_order_ceiling,
        "warnings": result.identifiability.warnings,
        "fit_result": str(path),
    }


def score(
    scenario: Scenario,
    measurements_path: str | Path,
    fit_path: str | Path,
    truth_path: str | Path,
    out_dir: str | Path,
) -> dict[str, Any]:
    ms = read_measurements(scenario, measurements_path)
    truth = load_truth(truth_path)
    fit_payload = json.loads(Path(fit_path).read_text(encoding="utf-8"))
    metrics = compare(fit_payload, truth, ms)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    path = out / METRICS_NAME
    path.write_text(
        json.dumps(metrics.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary = metrics.to_dict()
    summary["command"] = "score"
    summary["metrics"] = str(path)
    return summary


def report(
    scenario: Scenario, artifacts_dir: str | Path, report_dir: str | Path
) -> dict[str, Any]:
    from .report import fitted_values, plot_residual_map, plot_wafer_map, render_html

    in_dir = Path(artifacts_dir)
    ms = read_measurements(scenario, in_dir / MEASUREMENTS_NAME)
    fit_payload = json.loads((in_dir / FIT_RESULT_NAME).read_text(encoding="utf-8"))
    metrics_path = in_dir / METRICS_NAME
    metrics = json.loads(metrics_path.read_text(encoding="utf-8")) if metrics_path.exists() else {}

    out = Path(report_dir)
    dx_fit, dy_fit = fitted_values(fit_payload, ms)
    resid_x = ms.dx_nm - dx_fit
    resid_y = ms.dy_nm - dy_fit
    r = scenario.wafer_radius_mm
    e = scenario.wafer.edge_exclusion_mm

    figures = [
        plot_wafer_map(
            ms.dx_nm, ms.x_mm, ms.y_mm, out_path=out / "measured_dx.png",
            title="实测 dx 场", label="dx (nm)", radius_mm=r, edge_exclusion_mm=e,
        ),
        plot_wafer_map(
            dx_fit, ms.x_mm, ms.y_mm, out_path=out / "fitted_dx.png",
            title="拟合模型 dx 场（wafer + 场内）", label="dx (nm)",
            radius_mm=r, edge_exclusion_mm=e,
        ),
        plot_wafer_map(
            np.hypot(ms.dx_nm, ms.dy_nm), ms.x_mm, ms.y_mm,
            out_path=out / "overlay_quiver.png", title="overlay 矢量（长度 = 幅值）",
            radius_mm=r, edge_exclusion_mm=e, kind="quiver", dx=ms.dx_nm, dy=ms.dy_nm,
        ),
        plot_residual_map(
            resid_x, resid_y, ms.x_mm, ms.y_mm,
            out_path=out / "residual_map.png", radius_mm=r, edge_exclusion_mm=e,
        ),
    ]
    html = render_html(figures, metrics, out / "report.html", fit=fit_payload)
    return {
        "command": "report",
        "figures": [str(p) for p in figures],
        "html": str(html),
        "residual_rms_nm": round(float(np.sqrt(np.mean(resid_x**2))), 6),
    }


# ------------------------------------------------------------------ 串联


def run_all(
    config_path: str | Path, out_dir: str | Path, seed: int | None = None
) -> dict[str, Any]:
    """跑完整条流水线。三步共用同一个输出目录，中间产物路径由它自己推出。"""
    scenario = scenario_with_seed(config_path, seed)
    out = Path(out_dir)
    synth_summary = synth(scenario, out)
    fit_summary = fit(scenario, out / MEASUREMENTS_NAME, out)
    score_summary = score(
        scenario, out / MEASUREMENTS_NAME, out / FIT_RESULT_NAME, out / TRUTH_NAME, out
    )
    report_summary = report(scenario, out, out / REPORT_DIRNAME)
    return {
        "command": "run",
        "synth": synth_summary,
        "fit": fit_summary,
        "score": score_summary,
        "report": report_summary,
    }
