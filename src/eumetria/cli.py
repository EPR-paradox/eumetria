"""命令行入口（§8）。

    eumetria synth  -c scenarios/base-300mm.yaml -o out/base
    eumetria fit    -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
    eumetria score  --fit out/base/fit_result.json --truth out/base/truth.json \\
                    -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
    eumetria report -i out/base -c scenarios/base-300mm.yaml -o out/base/report
    eumetria run    -c scenarios/base-300mm.yaml -o out/base          # 上面四步串联

本层只做三件事：解析参数、调 `pipeline`、把异常映射成退出码。
**编排与计算都不在这里** —— 想让别的入口（main.py / notebook / 测试）复用同一条流水线，
调 `eumetria.pipeline.run_all()`，不要构造 argparse.Namespace 再喂进来。

约定：
- 退出码：0 成功；1 输入/schema 错误；2 计算错误（点数不足、降秩）；3 内部错误
- `--seed` 覆盖 scenario 里的 seed
- `--json` 把机器可读摘要打到 **stdout**；日志一律走 **stderr**（stdout 不许被污染）
- 生成物写进 `-o` 目录（`out/` 已在 .gitignore）

与详设 §8 的两处修订：`score` 与 `report` 增加了 `-c/--config`。
原因是 A-04（拟合残差与真值残差之差 < 0.01nm）**必须重建真值场**才能判定，而重建与读
CSV 都需要场景参数（硅片半径，§6.2.5 强制要求）。原契约里只给路径不给场景，A-04 无从计算。
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from . import pipeline
from .errors import (
    EumetriaError,
    InsufficientPointsError,
    InvalidMeasurementError,
    SchemaValidationError,
    SchemaVersionError,
)

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


# ------------------------------------------------------------------ 子命令


def _cmd_synth(args: argparse.Namespace) -> dict[str, Any]:
    return pipeline.synth(pipeline.scenario_with_seed(args.config, args.seed), args.output)


def _cmd_fit(args: argparse.Namespace) -> dict[str, Any]:
    scenario = pipeline.scenario_with_seed(args.config, args.seed)
    return pipeline.fit(scenario, args.input, args.output)


def _cmd_score(args: argparse.Namespace) -> dict[str, Any]:
    scenario = pipeline.scenario_with_seed(args.config, args.seed)
    return pipeline.score(scenario, args.input, args.fit, args.truth, args.output)


def _cmd_report(args: argparse.Namespace) -> dict[str, Any]:
    # 注意：这里**不**传 seed —— report 只读磁盘上的产物，seed 影响不到任何东西，
    # 所以 report 子命令也就没有 --seed。别为了"参数整齐"给它加一个空转的开关。
    scenario = pipeline.scenario_with_seed(args.config)
    return pipeline.report(scenario, args.input, args.output)


def _cmd_run(args: argparse.Namespace) -> dict[str, Any]:
    """四步串联。中间产物路径由 `-o` 目录推出，不让两个信息源有机会不一致。"""
    return pipeline.run_all(args.config, args.output, seed=args.seed)


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

    p_report = sub.add_parser("report", help="基础图 + 单文件 HTML")
    p_report.add_argument("-i", "--input", required=True, help="产物目录（含 measurements/fit_result/metrics）")
    p_report.add_argument("-c", "--config", required=True)
    p_report.add_argument("-o", "--output", required=True)
    p_report.set_defaults(func=_cmd_report)

    p_run = sub.add_parser("run", help="synth → fit → score → report")
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
