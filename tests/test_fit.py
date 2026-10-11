"""fit 层测试：联合求解、守卫、可辨识性诊断、A-02/A-03/A-04 的前置条件。"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pytest

from eumetria.errors import InsufficientPointsError, SchemaValidationError
from eumetria.fit import (
    FitConfig,
    build_design,
    field_terms_for,
    fit_field,
)
from eumetria.fit.linear_fit import _correlation_warnings
from eumetria.io import load_scenario
from eumetria.synth import generate

SCENARIO_PATH = Path("/home/j/myprojects/eumetria/scenarios/base-300mm.yaml")


def _clean_scenario():
    """无噪声、无 TIS 的 scenario —— 用来做「系数能否精确恢复」的强测试。"""
    sc = load_scenario(SCENARIO_PATH)
    return sc.model_copy(
        update={
            "sampling": sc.sampling.model_copy(update={"tmu_nm": 0.0}),
            "tis": sc.tis.model_copy(update={"amplitude_nm": 0.0}),
        }
    )


@pytest.fixture(scope="module")
def noisy_run():
    return generate(load_scenario(SCENARIO_PATH))


@pytest.fixture(scope="module")
def clean_run():
    return generate(_clean_scenario())


# --------------------------------------------------------- 联合求解（D-09）


def test_field_basis_excludes_constant_term():
    """D-09：场内基不含 (0,0)，常量归 wafer 的 Z1。"""
    assert "0,0" not in field_terms_for(2)
    assert field_terms_for(2) == ["1,0", "0,1", "2,0", "1,1", "0,2"]
    assert field_terms_for(1) == ["1,0", "0,1"]


def test_design_matrix_is_the_stacked_blocks(noisy_run, ):
    measurements = noisy_run.measurements
    config = FitConfig(max_order=12, field_order=2)
    X, w_terms, f_terms = build_design(measurements, config)
    assert X.shape == (len(measurements), 12 + 5)
    assert w_terms == [f"Z{i}" for i in range(1, 13)]
    assert f_terms == field_terms_for(2)


def test_clean_data_recovers_truth_exactly(clean_run):
    """无噪声时联合求解必须精确恢复到真值（联合解是线性模型的无偏解）。"""
    res = fit_field(clean_run.measurements, FitConfig(max_order=12, field_order=2))
    assert res.n_points == 188
    assert res.n_params == 17
    assert res.residual_rms_nm < 1e-10

    for term, value in clean_run.truth.wafer_terms.items():
        assert res.wafer_coeffs_x[term] == pytest.approx(value, abs=1e-9), term
    for term, value in clean_run.truth.field_terms_x.items():
        assert res.field_coeffs_x[term] == pytest.approx(value, abs=1e-9), term
    # 真值里没给的项必须是 0，不能凭空冒出来
    for term in ("Z1", "Z2", "Z3"):
        assert res.wafer_coeffs_x[term] == pytest.approx(0.0, abs=1e-9)


def test_two_channels_are_independent(clean_run):
    """x/y 是两套独立系数：dx 乘 3、dy 乘 −2，系数应各自跟着缩放。"""
    base = clean_run.measurements
    scaled = dataclasses.replace(
        base,
        dx_nm=base.dx_nm * 3.0,
        dx_nm_0=None,
        dx_nm_180=None,
        dy_nm=base.dy_nm * (-2.0),
        dy_nm_0=None,
        dy_nm_180=None,
    )
    res = fit_field(scaled, FitConfig(max_order=12, field_order=2))
    for term in res.wafer_terms:
        # coef_x = 3·orig 而 coef_y = −2·orig  ⇒  coef_x = −1.5·coef_y
        assert res.wafer_coeffs_x[term] == pytest.approx(-1.5 * res.wafer_coeffs_y[term], abs=1e-9)
    for term in res.field_terms:
        assert res.field_coeffs_x[term] == pytest.approx(-1.5 * res.field_coeffs_y[term], abs=1e-9)


# --------------------------------------------------------------- 守卫


def test_insufficient_points_raises_with_numbers(clean_run):
    """点数 < 参数数 → 显式抛错，且信息里带「需要多少 / 实际多少」。"""
    base = clean_run.measurements
    tiny = dataclasses.replace(
        base,
        point_id=base.point_id[:10],
        field_id=base.field_id[:10],
        x_mm=base.x_mm[:10],
        y_mm=base.y_mm[:10],
        xf_mm=base.xf_mm[:10],
        yf_mm=base.yf_mm[:10],
        rho=base.rho[:10],
        theta_rad=base.theta_rad[:10],
        dx_nm=base.dx_nm[:10],
        dy_nm=base.dy_nm[:10],
        is_edge=base.is_edge[:10],
        dx_nm_0=None,
        dx_nm_180=None,
        dy_nm_0=None,
        dy_nm_180=None,
    )
    with pytest.raises(InsufficientPointsError) as exc:
        fit_field(tiny, FitConfig(max_order=12, field_order=2))
    assert exc.value.required == 17
    assert exc.value.available == 10
    assert "17" in str(exc.value) and "10" in str(exc.value)


def test_rank_deficiency_raises_instead_of_min_norm_solution(clean_run):
    """降秩（点全在同一个"位置"）必须抛错 —— lstsq 在欠定系统上不会报错。"""
    base = clean_run.measurements
    n = 40
    degenerate = dataclasses.replace(
        base,
        point_id=np.arange(n, dtype=np.int64),
        field_id=np.arange(n, dtype=np.int64),  # 每点一个场，绕开「场内坐标不重复」检查
        x_mm=np.full(n, 100.0),
        y_mm=np.full(n, 0.0),
        xf_mm=np.full(n, 5.0),
        yf_mm=np.full(n, 5.0),
        rho=np.full(n, 100.0 / 150.0),
        theta_rad=np.zeros(n),
        dx_nm=np.linspace(0.0, 1.0, n),
        dy_nm=np.linspace(0.0, 1.0, n),
        is_edge=np.zeros(n, dtype=bool),
        dx_nm_0=None,
        dx_nm_180=None,
        dy_nm_0=None,
        dy_nm_180=None,
    )
    with pytest.raises(InsufficientPointsError, match="降秩"):
        fit_field(degenerate, FitConfig(max_order=12, field_order=2))


def test_config_validation():
    with pytest.raises(NotImplementedError, match="Phase 2"):
        FitConfig(selection="aic")
    with pytest.raises(SchemaValidationError, match="fixed_terms"):
        FitConfig(selection="fixed")
    with pytest.raises(SchemaValidationError, match="超出"):
        FitConfig(selection="fixed", fixed_terms=["Z4", "Z99"], max_order=12)
    with pytest.raises(SchemaValidationError, match="zernike"):
        FitConfig(basis="cartesian")
    with pytest.raises(SchemaValidationError, match="field_order"):
        FitConfig(field_order=0)


def test_fixed_selection_uses_the_given_terms(clean_run):
    config = FitConfig(selection="fixed", fixed_terms=["Z4", "Z11"], max_order=12, field_order=2)
    res = fit_field(clean_run.measurements, config)
    assert res.selected_terms == ["Z4", "Z11"]
    assert res.n_params == 2 + 5
    assert set(res.wafer_coeffs_x) == {"Z4", "Z11"}


# --------------------------------------------------- 可辨识性与警告


def test_correlation_warning_helper_flags_duplicate_columns():
    x = np.arange(10.0)
    X = np.column_stack([x, x, np.sin(x)])  # 前两列完全相同
    warnings = _correlation_warnings(X, ["A", "B", "C"])
    assert any("A 与 B" in w and "1.00" in w for w in warnings)
    # 只报真正相关的项对，不要把无关的也报进来
    assert not any("A 与 C" in w for w in warnings)


def test_noisy_fit_reports_sane_diagnostics(noisy_run):
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=2))
    assert res.identifiability.rank_deficient is False
    assert res.identifiability.max_identifiable_order == 12
    assert res.condition_number < 1e6
    assert not any("条件数" in w for w in res.identifiability.warnings)
    assert res.residual_rms_nm == pytest.approx(0.1 / np.sqrt(2.0), rel=0.2)


def test_field_diagnostics_implement_the_cpe_ceiling_rule(noisy_run):
    """§3.1.1：CPE 阶数上限由最差的场决定。默认场景实测 ceiling = 1。"""
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=2))
    diag = res.field_diagnostics
    assert len(diag) == 25

    partial = [d for d in diag.values() if d.max_identifiable_order == 1]
    full = [d for d in diag.values() if d.max_identifiable_order == 2]
    assert len(partial) == 4 and len(full) == 21
    assert {d.field_id for d in partial} == {6, 15, 60, 69}
    assert all(d.n_marks == 5 for d in partial)
    assert all(not d.is_identifiable for d in partial)
    assert all(d.is_identifiable for d in full)

    assert res.identifiability.cpe_order_ceiling == 1
    assert any("CPE 阶数上限" in w for w in res.identifiability.warnings)
    assert any("不可辨识" in w for w in res.identifiability.warnings)


def test_field_order_1_is_identifiable_everywhere(noisy_run):
    """一阶场内基应该全场可用（这正是"上限降到一阶"的落地验证）。"""
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=1))
    assert res.identifiability.cpe_order_ceiling == 1
    assert all(d.is_identifiable for d in res.field_diagnostics.values())
    assert not any("不可辨识" in w for w in res.identifiability.warnings)
    assert not any("CPE 阶数上限" in w for w in res.identifiability.warnings)


# ------------------------------------------------- A-03 / A-04 的前置条件


def test_a03_main_terms_recovered_within_10_percent(noisy_run):
    """A-03：TMU 0.1nm 下前 6 阶主要项相对误差 < 10%。"""
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=2))
    worst = 0.0
    for term, truth in noisy_run.truth.wafer_terms.items():
        fit = res.wafer_coeffs_x[term]
        rel = abs(fit - truth) / abs(truth)
        worst = max(worst, rel)
    assert worst < 0.10, f"最大相对误差 {worst:.4f}"


def test_a04_residual_rms_agrees_with_true_residual(noisy_run):
    """A-04：拟合残差 RMS 与真值残差 RMS 之差 < 0.01nm。

    注意等效噪声是 σ/√2（§3.4），且拟合吃掉 17 个自由度 → 残差略小于噪声。
    """
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=2))
    true_resid = noisy_run.measurements.dx_nm - noisy_run.truth_values[:, 0]
    true_rms = float(np.sqrt(np.mean(true_resid**2)))
    assert abs(res.residual_rms_nm_x - true_rms) < 0.01


def test_result_serialises_to_the_documented_json(noisy_run):
    res = fit_field(noisy_run.measurements, FitConfig(max_order=12, field_order=2))
    payload = json.loads(json.dumps(res.to_dict()))
    for key in (
        "version",
        "basis",
        "max_order",
        "selected_terms",
        "selection_method",
        "wafer_coeffs_x",
        "wafer_coeffs_y",
        "field_coeffs_x",
        "field_coeffs_y",
        "residual_rms_nm",
        "n_points",
        "n_params",
        "condition_number",
        "identifiability",
    ):
        assert key in payload, key
    assert payload["identifiability"]["cpe_order_ceiling"] == 1
    assert len(payload["field_diagnostics"]) == 25
