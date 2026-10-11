"""report 层测试：基础图与单文件 HTML。"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from eumetria.cli import EXIT_OK, main
from eumetria.io import load_scenario
from eumetria.report import fitted_values, plot_wafer_map, render_html

SCENARIO_PATH = Path("/home/j/myprojects/eumetria/scenarios/base-300mm.yaml")


@pytest.fixture(scope="module")
def artifacts(tmp_path_factory):
    out = tmp_path_factory.mktemp("out")
    assert main(["run", "-c", str(SCENARIO_PATH), "-o", str(out)]) == EXIT_OK
    return out


def test_run_produces_the_four_figures_and_html(artifacts):
    report = artifacts / "report"
    for name in (
        "measured_dx.png",
        "fitted_dx.png",
        "overlay_quiver.png",
        "residual_map.png",
        "report.html",
    ):
        assert (report / name).exists(), name
    for name in ("measured_dx.png", "residual_map.png"):
        assert (report / name).stat().st_size > 5000  # 不是空图


def test_html_is_self_contained(artifacts):
    """P-05 的「单文件 HTML」：图以 base64 内嵌，没有外部引用。"""
    text = (artifacts / "report" / "report.html").read_text(encoding="utf-8")
    assert text.count("data:image/png;base64,") == 4
    assert "<img" in text and "src=\"./" not in text
    assert "A-03" in text and "CPE 阶数上限" in text


def test_html_shows_the_verdicts_and_notes(artifacts):
    metrics = json.loads((artifacts / "metrics.json").read_text(encoding="utf-8"))
    text = (artifacts / "report" / "report.html").read_text(encoding="utf-8")
    if metrics["a02_applicable"] is False:
        assert "不适用：数据带噪" in text
    for note in metrics["notes"]:
        assert note[:20] in text


def test_fitted_values_match_the_fit_result(artifacts):
    """报告里的"拟合场"必须与 fit 写出的模型一致（否则残差图就是假的）。"""
    scenario = load_scenario(SCENARIO_PATH)
    from eumetria.pipeline import read_measurements

    ms = read_measurements(scenario, artifacts / "measurements.csv")
    fit = json.loads((artifacts / "fit_result.json").read_text(encoding="utf-8"))
    dx_fit, _dy_fit = fitted_values(fit, ms)
    resid_x = ms.dx_nm - dx_fit
    assert float(np.sqrt(np.mean(resid_x**2))) == pytest.approx(fit["residual_rms_nm_x"], rel=1e-3)
    assert dx_fit.shape == ms.dx_nm.shape


def test_plot_helpers_validate_input(tmp_path):
    with pytest.raises(ValueError, match="quiver"):
        plot_wafer_map(
            np.zeros(3), np.zeros(3), np.zeros(3),
            out_path=tmp_path / "x.png", kind="quiver",
        )
    with pytest.raises(ValueError, match="未知的 kind"):
        plot_wafer_map(
            np.zeros(3), np.zeros(3), np.zeros(3), out_path=tmp_path / "x.png", kind="labels"
        )


def test_render_html_without_fit_or_metrics(tmp_path):
    """缺 fit_result / metrics 时也要能出报告骨架，而不是崩。"""
    png = tmp_path / "a.png"
    png.write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 100)
    out = render_html([png], {}, tmp_path / "r.html")
    text = out.read_text(encoding="utf-8")
    assert "验收判定" in text and "data:image/png;base64," in text


def test_cli_report_alone_on_existing_artifacts(artifacts, tmp_path):
    out = tmp_path / "report2"
    code = main(
        ["report", "-i", str(artifacts), "-c", str(SCENARIO_PATH), "-o", str(out)]
    )
    assert code == EXIT_OK
    assert (out / "report.html").exists()
    assert (out / "residual_map.png").exists()


def test_cli_report_missing_artifacts_is_input_error(tmp_path):
    code = main(
        [
            "report",
            "-i", str(tmp_path / "nope"),
            "-c", str(SCENARIO_PATH),
            "-o", str(tmp_path / "out"),
        ]
    )
    assert code == 1
