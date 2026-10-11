"""CLI 集成测试：A-01、A-05 与退出码映射（§8）。

直接调用 `main(argv)` 而不是起子进程 —— 快，且 capsys 能检查 stdout/stderr 的分工。
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from eumetria.cli import EXIT_COMPUTE, EXIT_INPUT, EXIT_OK, main

SCENARIO_PATH = Path("/home/j/myprojects/eumetria/scenarios/base-300mm.yaml")


@pytest.fixture()
def workspace(tmp_path) -> Path:
    return tmp_path / "out"


# ------------------------------------------------------------------ A-01


def test_synth_produces_both_files(tmp_path, workspace):
    """A-01：给定 scenario + seed，一条命令生成 measurements.csv 与 truth.json。"""
    code = main(["synth", "-c", str(SCENARIO_PATH), "-o", str(workspace)])
    assert code == EXIT_OK
    assert (workspace / "measurements.csv").exists()
    assert (workspace / "truth.json").exists()
    with (workspace / "measurements.csv").open(encoding="utf-8") as fh:
        assert len(fh.readlines()) == 1 + 188


# ------------------------------------------------------------------ A-05


def test_same_seed_is_bit_identical(tmp_path):
    """A-05：同一 seed 两次运行产物 bit 级一致。"""
    a, b = tmp_path / "a", tmp_path / "b"
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(a)]) == EXIT_OK
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(b)]) == EXIT_OK
    for name in ("measurements.csv",):
        assert (a / name).read_bytes() == (b / name).read_bytes()
    ta = json.loads((a / "truth.json").read_text(encoding="utf-8"))
    tb = json.loads((b / "truth.json").read_text(encoding="utf-8"))
    ta.pop("generated_at")  # 时间戳必然不同
    tb.pop("generated_at")
    assert ta == tb


def test_seed_override_changes_the_data(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(a)]) == EXIT_OK
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(b), "--seed", "42"]) == EXIT_OK
    assert (a / "measurements.csv").read_bytes() != (b / "measurements.csv").read_bytes()


# --------------------------------------------------------- fit / score


def test_fit_then_score_end_to_end(tmp_path, workspace):
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(workspace)]) == EXIT_OK
    assert main(
        ["fit", "-i", str(workspace / "measurements.csv"), "-c", str(SCENARIO_PATH),
         "-o", str(workspace)]
    ) == EXIT_OK
    fit_payload = json.loads((workspace / "fit_result.json").read_text(encoding="utf-8"))
    assert fit_payload["n_points"] == 188
    assert fit_payload["identifiability"]["cpe_order_ceiling"] == 1

    assert main(
        ["score", "--fit", str(workspace / "fit_result.json"),
         "--truth", str(workspace / "truth.json"),
         "-i", str(workspace / "measurements.csv"), "-c", str(SCENARIO_PATH),
         "-o", str(workspace)]
    ) == EXIT_OK
    metrics = json.loads((workspace / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["pass_A03"] is True
    assert metrics["pass_A02"] is False
    assert metrics["residual_agreement_nm"] < 0.01


def test_run_chains_everything(tmp_path, workspace):
    assert main(["run", "-c", str(SCENARIO_PATH), "-o", str(workspace)]) == EXIT_OK
    for name in ("measurements.csv", "truth.json", "fit_result.json", "metrics.json"):
        assert (workspace / name).exists(), name


# ------------------------------------------------------------ 退出码


def test_missing_config_is_input_error(tmp_path, workspace):
    assert main(["synth", "-c", str(tmp_path / "nope.yaml"), "-o", str(workspace)]) == EXIT_INPUT


def test_too_few_points_is_compute_error(tmp_path, workspace):
    """点数 < 参数数 → 退出码 2（不是崩溃、也不是静默给解）。"""
    text = SCENARIO_PATH.read_text(encoding="utf-8").replace("fields: 25", "fields: 1").replace(
        "marks_per_field: 8", "marks_per_field: 4"
    )
    cfg = tmp_path / "tiny.yaml"
    cfg.write_text(text, encoding="utf-8")
    assert main(["synth", "-c", str(cfg), "-o", str(workspace)]) == EXIT_OK
    assert main(
        ["fit", "-i", str(workspace / "measurements.csv"), "-c", str(cfg), "-o", str(workspace)]
    ) == EXIT_COMPUTE


def test_tampered_csv_is_input_error(tmp_path, workspace):
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(workspace)]) == EXIT_OK
    path = workspace / "measurements.csv"
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = lines[1].split(",")
    fields[8] = "0.5"  # 篡改 rho -> §6.2.5 第 1 条拦截
    lines[1] = ",".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert main(
        ["fit", "-i", str(path), "-c", str(SCENARIO_PATH), "-o", str(workspace)]
    ) == EXIT_INPUT


# ------------------------------------------------------ stdout/stderr


def test_json_goes_to_stdout_and_logs_to_stderr(capsys, workspace):
    code = main(["--json", "synth", "-c", str(SCENARIO_PATH), "-o", str(workspace)])
    assert code == EXIT_OK
    captured = capsys.readouterr()
    payload = json.loads(captured.out)  # stdout 必须是干净的 JSON
    assert payload["command"] == "synth"
    assert payload["n_points"] == 188
    assert captured.err == ""


def test_human_summary_goes_to_stderr(capsys, workspace):
    assert main(["synth", "-c", str(SCENARIO_PATH), "-o", str(workspace)]) == EXIT_OK
    captured = capsys.readouterr()
    assert captured.out == ""  # 没有 --json 时不许往 stdout 写
    assert "n_points: 188" in captured.err
