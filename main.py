#!/usr/bin/env python
"""一条命令跑完整条流水线 —— 零参数即可。

PyCharm 里右键点绿三角、或者在终端敲 `python main.py` 都直接跑，
默认用 `scenarios/base-300mm.yaml`，产物写进 `out/base/`。

只用标准库 + 本包。刻意**不**引入 argparse 子命令那一套：想要完整的命令行界面，
用 `eumetria` 命令（或 `python -m eumetria.cli`）；这里只图一件事 —— 打开就能跑。
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))  # 免安装也能直接跑（PyCharm 新克隆时很有用）

from eumetria import pipeline  # noqa: E402  (必须在 sys.path 调整之后)
from eumetria.errors import (  # noqa: E402
    InsufficientPointsError,
    InvalidMeasurementError,
    SchemaValidationError,
    SchemaVersionError,
)

DEFAULT_CONFIG = ROOT / "scenarios" / "base-300mm.yaml"
DEFAULT_OUT = ROOT / "out" / "base"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="eumetria 全流程（synth → fit → score → report）",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-c", "--config", default=str(DEFAULT_CONFIG), help="场景 YAML")
    parser.add_argument("-o", "--output", default=str(DEFAULT_OUT), help="输出目录")
    parser.add_argument("--seed", type=int, default=None, help="覆盖场景里的随机种子")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    print(f"场景   {args.config}")
    print(f"输出   {args.output}")
    try:
        summary = pipeline.run_all(args.config, args.output, seed=args.seed)
    except (SchemaValidationError, SchemaVersionError, InvalidMeasurementError, FileNotFoundError) as exc:
        print(f"输入/schema 错误：{exc}", file=sys.stderr)
        return 1
    except InsufficientPointsError as exc:
        print(f"计算错误：{exc}", file=sys.stderr)
        return 2
    except Exception as exc:  # noqa: BLE001
        print(f"内部错误：{type(exc).__name__}: {exc}", file=sys.stderr)
        return 3

    synth, fit, score, report = (
        summary["synth"], summary["fit"], summary["score"], summary["report"],
    )
    print()
    print(f"1 synth   {synth['n_points']} 点 / {synth['n_fields']} 场 → {synth['measurements']}")
    print(f"2 fit     {fit['n_params']} 参数 / 条件数 {fit['condition_number']} / "
          f"残差 RMS {fit['residual_rms_nm']} nm")
    for warning in fit["warnings"]:
        print(f"          告警 {warning}")
    print(f"3 score   A-02 {'不适用' if not score['a02_applicable'] else ('通过' if score['pass_A02'] else '未过')} / "
          f"A-03 {'通过' if score['pass_A03'] else '未过'} / "
          f"残差一致性 {score['residual_agreement_nm']} nm")
    print(f"4 report  {len(report['figures'])} 张图 + {report['html']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
