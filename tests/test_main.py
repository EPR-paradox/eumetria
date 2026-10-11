"""main.py 入口测试 —— 零参数直接跑是它存在的唯一理由，所以要真的跑一遍。"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:  # 让 `import main` 不依赖调用方式（python -m pytest / pytest 都行）
    sys.path.insert(0, str(ROOT))

MAIN = ROOT / "main.py"
SCENARIO_PATH = ROOT / "scenarios" / "base-300mm.yaml"


def test_main_runs_with_no_arguments_in_a_subprocess(tmp_path):
    """PyCharm 绿三角的等价物：`python main.py` 不带任何参数。

    必须走子进程 —— 目标之一就是验证"没装过包、没设 PYTHONPATH 也能跑"，
    进程内 import 会把 `sys.path` 的处理绕过去。
    """
    proc = subprocess.run(
        [sys.executable, str(MAIN), "-o", str(tmp_path / "out")],
        capture_output=True, text=True, cwd=str(tmp_path), timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    assert "1 synth" in proc.stdout and "4 report" in proc.stdout
    for name in ("measurements.csv", "truth.json", "fit_result.json", "metrics.json"):
        assert (tmp_path / "out" / name).exists(), name
    assert (tmp_path / "out" / "report" / "report.html").exists()


def test_main_default_paths_point_into_the_repo(tmp_path):
    """默认值要对准仓库自身，而不是当前工作目录 —— 从别处调用也得对。"""
    from main import DEFAULT_CONFIG, DEFAULT_OUT

    assert DEFAULT_CONFIG == ROOT / "scenarios" / "base-300mm.yaml"
    assert DEFAULT_OUT == ROOT / "out" / "base"
    assert DEFAULT_CONFIG.exists()


def test_main_respects_seed_and_writes_elsewhere(tmp_path):
    from main import main

    out_a, out_b = tmp_path / "a", tmp_path / "b"
    assert main(["-c", str(SCENARIO_PATH), "-o", str(out_a), "--seed", "7"]) == 0
    assert main(["-c", str(SCENARIO_PATH), "-o", str(out_b), "--seed", "7"]) == 0
    assert (out_a / "measurements.csv").read_bytes() == (out_b / "measurements.csv").read_bytes()
    assert main(["-c", str(SCENARIO_PATH), "-o", str(tmp_path / "c"), "--seed", "8"]) == 0
    assert (out_a / "measurements.csv").read_bytes() != (
        tmp_path / "c" / "measurements.csv"
    ).read_bytes()


def test_main_reports_input_errors_with_exit_code_1(tmp_path):
    from main import main

    assert main(["-c", str(tmp_path / "missing.yaml"), "-o", str(tmp_path / "o")]) == 1


def test_pipeline_run_all_is_reusable_without_argparse():
    """main.py 与 cli 共用的那一层：不构造 argparse.Namespace 也能跑完整条流水线。"""
    from eumetria import pipeline

    assert callable(pipeline.run_all)
    params = list(pipeline.run_all.__code__.co_varnames[: pipeline.run_all.__code__.co_argcount])
    assert params == ["config_path", "out_dir", "seed"]


def test_cli_and_pipeline_agree_on_the_numbers(tmp_path):
    """两个前端调同一条流水线，数字必须一致 —— 否则说明编排又被复制了一份。"""
    from eumetria import pipeline
    from eumetria.cli import EXIT_OK, main

    assert main(["run", "-c", str(SCENARIO_PATH), "-o", str(tmp_path / "via_cli")]) == EXIT_OK
    via_lib = pipeline.run_all(SCENARIO_PATH, tmp_path / "via_lib")
    a = (tmp_path / "via_cli" / "fit_result.json").read_text(encoding="utf-8")
    b = (tmp_path / "via_lib" / "fit_result.json").read_text(encoding="utf-8")
    assert a == b
    assert via_lib["fit"]["n_params"] > 0
    for key in ("synth", "fit", "score", "report"):
        assert key in via_lib


@pytest.mark.parametrize("name", ["pipeline", "cli", "report"])
def test_public_modules_stay_importable_without_matplotlib(name):
    """报告是可选依赖：import 本包不该因为缺 matplotlib 就炸。"""
    import importlib

    module = importlib.import_module(f"eumetria.{name}")
    assert module is not None
