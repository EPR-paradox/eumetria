"""metrics 层测试：系数恢复、残差一致性、A-02/A-03 判定。"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from eumetria.fit import FitConfig, fit_field
from eumetria.io import load_scenario
from eumetria.metrics import coefficient_recovery, compare, residual_rms
from eumetria.synth import generate

SCENARIO_PATH = Path("/home/j/myprojects/eumetria/scenarios/base-300mm.yaml")


def _scenario(tmu_nm: float, tis_amplitude: float = 0.0):
    sc = load_scenario(SCENARIO_PATH)
    return sc.model_copy(
        update={
            "sampling": sc.sampling.model_copy(update={"tmu_nm": tmu_nm}),
            "tis": sc.tis.model_copy(update={"amplitude_nm": tis_amplitude}),
        }
    )


@pytest.fixture(scope="module")
def clean():
    sc = _scenario(0.0)
    run = generate(sc)
    return run, fit_field(run.measurements, FitConfig(max_order=12, field_order=2))


@pytest.fixture(scope="module")
def noisy():
    sc = load_scenario(SCENARIO_PATH)
    run = generate(sc)
    return run, fit_field(run.measurements, FitConfig(max_order=12, field_order=2))


# ---------------------------------------------------------- 基础工具


def test_coefficient_recovery_reports_all_four_numbers():
    out = coefficient_recovery({"Z4": 2.4, "Z5": 0.0}, {"Z4": 2.5, "Z5": 0.0})
    assert out["Z4"] == {
        "truth": 2.5,
        "fit": 2.4,
        "abs_err": pytest.approx(0.1),
        "rel_err": pytest.approx(0.04),
    }
    assert np.isnan(out["Z5"]["rel_err"])  # 真值为 0 时不做相对误差（会除零）


def test_residual_rms_known_value():
    assert residual_rms(np.array([1.0, 2.0]), np.array([1.0, 4.0])) == pytest.approx(np.sqrt(2.0))


# ------------------------------------------------------ 无噪声：全过


def test_clean_run_passes_a02_and_a03(clean):
    run, fit = clean
    m = compare(fit, run.truth, run.measurements)
    assert m.pass_A02 is True
    assert m.pass_A03 is True
    assert m.max_rel_err < 1e-9
    assert m.residual_rms_nm < 1e-10
    assert m.residual_agreement_nm < 1e-10
    assert m.spurious_max_abs_nm < 1e-9


# ------------------------------------------------------ 有噪声：A-03 过、A-02 挂


def test_noisy_run_fails_a02_but_passes_a03(noisy):
    run, fit = noisy
    m = compare(fit, run.truth, run.measurements)
    assert m.a02_applicable is False
    assert m.pass_A02 is False, "带噪数据上 A-02 的前提不成立，不该判过"
    assert m.pass_A03 is True
    assert m.max_rel_err < 0.10
    assert any("A-02 不适用" in n for n in m.notes)
    # A-04 的核心：拟合残差 RMS 与真值残差 RMS 之差 < 0.01nm
    assert m.residual_agreement_nm < 0.01
    assert m.residual_rms_truth_nm == pytest.approx(0.1 / np.sqrt(2.0), rel=0.1)


def test_a03_only_judges_the_first_six_wafer_terms(noisy):
    """A-03 的口径：真值里 j ≤ 6 的非零项。场内高阶项差得多也不该拖累它。"""
    run, fit = noisy
    m = compare(fit, run.truth, run.measurements)
    assert m.max_rel_err_field > m.max_rel_err  # 场内 (2,0) 确实更差（见 Q-04）
    assert m.pass_A03 is True


# ------------------------------------------------------ 假阳性与边角


def test_spurious_terms_are_noted(noisy):
    run, fit = noisy
    payload = fit.to_dict()
    payload["wafer_coeffs_x"]["Z9"] = 0.42  # 真值里 Z9 = 0
    m = compare(payload, run.truth, run.measurements)
    assert m.spurious_max_abs_nm == pytest.approx(0.42)
    assert any("Z9" in n for n in m.notes)


def test_compare_accepts_plain_dict(noisy):
    """score 命令从 fit_result.json 读回来的就是 dict（系数被 round 到 6 位）。"""
    run, fit = noisy
    from_dict = compare(fit.to_dict(), run.truth, run.measurements)
    from_obj = compare(fit, run.truth, run.measurements)
    assert from_dict.max_rel_err == pytest.approx(from_obj.max_rel_err, rel=1e-4)
    assert from_dict.pass_A03 == from_obj.pass_A03


def test_metrics_json_shape(noisy):
    run, fit = noisy
    payload = compare(fit, run.truth, run.measurements).to_dict()
    for key in (
        "version",
        "coefficient_recovery",
        "max_rel_err",
        "max_rel_err_field",
        "residual_rms_nm",
        "residual_rms_truth_nm",
        "residual_agreement_nm",
        "pass_A02",
        "pass_A03",
        "a02_applicable",
        "notes",
    ):
        assert key in payload, key
    assert "Z4" in payload["coefficient_recovery"]
    assert "field 1,0" in payload["coefficient_recovery"]
