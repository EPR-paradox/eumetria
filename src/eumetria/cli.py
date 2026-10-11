"""命令行入口（§8）。

    eumetria synth -c scenarios/base-300mm.yaml -o out/base
    eumetria fit   -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
    eumetria score --fit out/base/fit_result.json --truth out/base/truth.json \
                   -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
    eumetria run   -c scenarios/base-300mm.yaml -o out/base          # 上面四步串联

约定：
- 退出码：0 成功；1 输入/schema 错误；2 计算错误（点数不足、降秩）；3 内部错误
- `--seed` 覆盖 scenario 里的 seed
- `--json` 把机器可读摘要打到 **stdout**；日志一律走 **stderr**（stdout 不许被污染）
- 生成物写进 `-o` 目录（`out/` 已在 .gitignore）

与详设 §8 的一处修订：`score` 增加了 `-i/--input` 与 `-c/--config`。
原因是 A-04（拟合残差与真值残差之差 < 0.01nm）**必须重建真值场**才能判定，
而它需要 measurements（点坐标）与 scenario（读 CSV 时的硅片半径，§6.2.5 强制要求）。
原契约里 score 只给 fit/truth 两个路径，A-04 无从计算。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .errors import (
    EumetriaError,
    InsufficientPointsError,
    InvalidMeasurementError,
    SchemaValidationError,
    SchemaVersionError,
)
from .fit import FitConfig, fit_field, write_fit_result
from .io import load_measurements, load_scenario, load_truth
from .io.writers import write_measurements, write_truth
from .io.schemas import Scenario
from .metrics import compare
from .synth import generate

EXIT_OK = 0
EXIT_INPUT = 1
EXIT_COMPUTE = 2
EXIT_INTERNAL = 3


def _log(message: str) -> None:
    print(message, file=sys.stderr)


def _emit(summary: dict[str, Any], as_json: bool) -> None:
    if as_json:
        print(json.dumps(summary, ensure_ascii=False, indent=2))
    else:
        for key, value in summary.items():
            _log(f"{key}: {value}")


def _scenario(args: argparse.Namespace) -> Scenario:
    scenario = load_scenario(args.config)
    seed = getattr(args, "seed", None)
    if seed is not None:
        scenario = scenario.model_copy(update={"seed": int(seed)})
    return scenario


def _measurements(scenario: Scenario, path: Path):
    return load_measurements(
        path,
        wafer_radius_mm=scenario.wafer_radius_mm,
        field_width_mm=scenario.field.width_mm,
        field_height_mm=scenario.field.height_mm,
    )


def _fit_config(scenario: Scenario, args: argparse.Namespace) -> FitConfig:
    return FitConfig(
        basis=scenario.fingerprint.basis,
        max_order=scenario.fingerprint.max_order,
        field_order=scenario.fingerprint.field_order,
    )


# ------------------------------------------------------------------ 子命令


def _cmd_synth(args: argparse.Namespace) -> dict[str, Any]:
    scenario = _scenario(args)
    run = generate(scenario)
    out = Path(args.output)
    ms_path = write_measurements(out / "measurements.csv", run.measurements)
    truth_path = write_truth(out / "truth.json", run.truth)
    return {
        "command": "synth",
        "seed": scenario.seed,
        "n_points": len(run.measurements),
        "n_fields": int(len(set(int(v) for v in run.measurements.field_id))),
        "measurements": str(ms_path),
        "truth": str(truth_path),
        "tis_estimate_nm": [round(v, 4) for v in run.tis_estimate_nm],
    }


def _cmd_fit(args: argparse.Namespace) -> dict[str, Any]:
    scenario = _scenario(args)
    ms = _measurements(scenario, Path(args.input))
    result = fit_field(ms, _fit_config(scenario, args))
    path = write_fit_result(Path(args.output) / "fit_result.json", result)
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


def _cmd_score(args: argparse.Namespace) -> dict[str, Any]:
    scenario = _scenario(args)
    ms = _measurements(scenario, Path(args.input))
    truth = load_truth(args.truth)
    fit_payload = json.loads(Path(args.fit).read_text(encoding="utf-8"))
    metrics = compare(fit_payload, truth, ms)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "metrics.json"
    path.write_text(
        json.dumps(metrics.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    summary = metrics.to_dict()
    summary["command"] = "score"
    summary["metrics"] = str(path)
    return summary


def _cmd_run(args: argparse.Namespace) -> dict[str, Any]:
    """synth → fit → score 串联（report 属 Phase 4）。

    四步共用同一个 `-o` 目录，所以 fit/score 需要的 `-i/--fit/--truth` 由 `run` 自己
    按固定文件名推出来 —— 不要求用户再写一遍，也不让两个信息源有机会不一致。
    """
    out = Path(args.output)
    step = argparse.Namespace(**vars(args))
    step.input = str(out / "measurements.csv")
    step.fit = str(out / "fit_result.json")
    step.truth = str(out / "truth.json")
    synth = _cmd_synth(step)
    fit = _cmd_fit(step)
    score = _cmd_score(step)
    return {
        "command": "run",
        "synth": synth,
        "fit": fit,
        "score": score,
    }


# ------------------------------------------------------------------ 解析


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="eumetria", description="光刻量测场分析")
    parser.add_argument("--json", action="store_true", help="机器可读摘要打到 stdout")
    sub = parser.add_subparsers(dest="command", required=True)

    p_synth = sub.add_parser("synth", help="正向生成 measurements.csv + truth.json")
    p_synth.add_argument("-c", "--config", required=True)
    p_synth.add_argument("-o", "--output", required=True)
    p_synth.add_argument("--seed", type=int, default=None)
    p_synth.set_defaults(func=_cmd_synth)

    p_fit = sub.add_parser("fit", help="联合最小二乘拟合")
    p_fit.add_argument("-i", "--input", required=True)
    p_fit.add_argument("-c", "--config", required=True)
    p_fit.add_argument("-o", "--output", required=True)
    p_fit.add_argument("--seed", type=int, default=None)
    p_fit.set_defaults(func=_cmd_fit)

    p_score = sub.add_parser("score", help="与真值对比（唯一读 truth.json 的命令）")
    p_score.add_argument("--fit", required=True)
    p_score.add_argument("--truth", required=True)
    p_score.add_argument("-i", "--input", required=True)
    p_score.add_argument("-c", "--config", required=True)
    p_score.add_argument("-o", "--output", required=True)
    p_score.add_argument("--seed", type=int, default=None)
    p_score.set_defaults(func=_cmd_score)

    p_run = sub.add_parser("run", help="synth → fit → score")
    p_run.add_argument("-c", "--config", required=True)
    p_run.add_argument("-o", "--output", required=True)
    p_run.add_argument("--seed", type=int, default=None)
    p_run.set_defaults(func=_cmd_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        summary = args.func(args)
    except (SchemaValidationError, SchemaVersionError, InvalidMeasurementError, FileNotFoundError) as exc:
        _log(f"输入/schema 错误：{exc}")
        return EXIT_INPUT
    except InsufficientPointsError as exc:
        _log(f"计算错误：{exc}")
        return EXIT_COMPUTE
    except EumetriaError as exc:  # pragma: no cover - 兜底
        _log(f"错误：{exc}")
        return EXIT_INTERNAL
    except Exception as exc:  # noqa: BLE001  pragma: no cover - 兜底
        _log(f"内部错误：{type(exc).__name__}: {exc}")
        return EXIT_INTERNAL
    _emit(summary, bool(getattr(args, "json", False)))
    return EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
