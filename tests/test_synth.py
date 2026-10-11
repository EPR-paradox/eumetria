"""synth 层测试：采样、点表、指纹、TIS/噪声（§6.2.4、§7.3）。"""

from __future__ import annotations

import numpy as np
import pytest

from eumetria.basis import CartesianBasis, ZernikeBasis
from eumetria.geometry import build_layout, field_marks, to_polar
from eumetria.synth import (
    add_noise,
    apply_tis,
    combine_tis,
    make_fingerprint,
    sample_fields,
    sample_points,
    tis_vector_from_polar,
)

WAFER_MM = 300.0
R_WAFER_MM = 150.0
EDGE_MM = 3.0
FIELD_W_MM = 26.0
FIELD_H_MM = 33.0
SEED = 20261003


def _layout(include_partial: bool = True):
    return build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=include_partial)


def _sampled(n_fields: int = 25, mode: str = "regular"):
    return sample_fields(_layout(), n_fields, mode, np.random.default_rng(SEED))


# ------------------------------------------------------------------ 采样


def test_regular_sampling_is_spatially_spread():
    """regular 采样不能把点全赶到某一条带上（等距索引抽取的典型失败模式）。"""
    fields = _sampled(25)
    radii = np.hypot([f.cx_mm for f in fields], [f.cy_mm for f in fields])
    rho = radii / R_WAFER_MM
    assert len(fields) == 25
    assert len({f.field_id for f in fields}) == 25
    assert rho.min() < 0.5, f"抽到的 field 全在外圈（min rho = {rho.min():.3f}）"
    assert rho.max() > 0.85, f"没抽到边缘场（max rho = {rho.max():.3f}）"
    # 四个象限都要有代表
    quadrants = {(np.sign(f.cx_mm), np.sign(f.cy_mm)) for f in fields}
    assert len(quadrants) == 4


def test_sampling_modes_and_determinism():
    a = _sampled(25, "regular")
    b = _sampled(25, "regular")
    assert [f.field_id for f in a] == [f.field_id for f in b]

    rng1 = np.random.default_rng(7)
    rng2 = np.random.default_rng(7)
    s1 = sample_fields(_layout(), 25, "stratified", rng1)
    s2 = sample_fields(_layout(), 25, "stratified", rng2)
    assert [f.field_id for f in s1] == [f.field_id for f in s2]

    # random 与 stratified 用不同的随机流，结果不应恰好相同
    assert [f.field_id for f in s1] != [f.field_id for f in _sampled(25, "random")]


def test_sampling_validation():
    layout = _layout()
    with pytest.raises(ValueError, match="n_fields"):
        sample_fields(layout, 0, "regular", np.random.default_rng(1))
    with pytest.raises(ValueError, match="n_fields"):
        sample_fields(layout, len(layout) + 1, "regular", np.random.default_rng(1))
    with pytest.raises(ValueError, match="未知的采样模式"):
        sample_fields(layout, 5, "spiral", np.random.default_rng(1))
    with pytest.raises(ValueError, match="layout 为空"):
        sample_fields([], 1, "regular", np.random.default_rng(1))


def test_all_fields_can_be_selected():
    layout = _layout()
    picked = sample_fields(layout, len(layout), "regular", np.random.default_rng(SEED))
    assert len(picked) == len(layout)
    assert len({f.field_id for f in picked}) == len(layout)


# --------------------------------------------------------------- 点表


def test_point_table_filters_edge_exclusion():
    fields = _sampled(25)
    points = sample_points(fields, 8, R_WAFER_MM, EDGE_MM)
    assert len(points) <= 25 * 8
    assert len(points) < 25 * 8, "边缘场的部分 mark 必须被剔除，否则过滤没生效"
    assert points.rho.max() <= (R_WAFER_MM - EDGE_MM) / R_WAFER_MM + 1e-12
    assert points.rho.min() >= 0.0
    assert np.all(points.theta_rad >= -np.pi) and np.all(points.theta_rad < np.pi)
    # point_id 连续且与行号一致
    assert list(points.point_id) == list(range(len(points)))
    # 场内坐标 = 全场坐标 − field 中心
    assert np.allclose(points.xi_eta()[0], 2 * points.xf_mm / FIELD_W_MM)


def test_point_table_matches_documented_golden_field():
    """详设 §6.2.3 的样例：field 中心 (-117, -82.5) 只留下 5 个 mark，ρ 与文档一致。

    这一条同时验证了 mark 生成（4 角 + 4 边中点、距边 3mm）、极坐标换算与边缘过滤。
    """
    fields = _layout()
    target = next(f for f in fields if (f.cx_mm, f.cy_mm) == (-117.0, -82.5))
    marks = field_marks(target, 8)
    assert len(marks) == 8
    xs = np.array([m.x_mm for m in marks])
    ys = np.array([m.y_mm for m in marks])
    rho, theta = to_polar(xs, ys, R_WAFER_MM)
    keep = rho <= (R_WAFER_MM - EDGE_MM) / R_WAFER_MM + 1e-12
    assert keep.sum() == 5, "文档记载该场只剩 5 个可用 mark"

    got = {(round(float(x), 3), round(float(y), 3)) for x, y in zip(xs[keep], ys[keep], strict=True)}
    assert got == {
        (-107.0, -96.0),
        (-127.0, -69.0),
        (-107.0, -69.0),
        (-117.0, -69.0),
        (-107.0, -82.5),
    }
    got_rho = {round(float(r), 4) for r in rho[keep]}
    assert got_rho == {0.9584, 0.9636, 0.8488, 0.9055, 0.9007}
    # θ 按具体 mark 核对：(-127, -69) 对应文档里的 -2.6439
    idx = int(np.flatnonzero(np.isclose(xs, -127.0) & np.isclose(ys, -69.0))[0])
    assert keep[idx]
    assert round(float(rho[idx]), 4) == 0.9636
    assert round(float(theta[idx]), 4) == -2.6439


def test_point_table_validation():
    with pytest.raises(ValueError, match="fields 为空"):
        sample_points([], 8, R_WAFER_MM, EDGE_MM)
    with pytest.raises(ValueError, match="可用半径非正"):
        sample_points(_layout(), 8, 3.0, 3.0)


# ------------------------------------------------------------- 指纹


def test_fingerprint_is_linear_in_coefficients():
    points = sample_points(_sampled(25), 8, R_WAFER_MM, EDGE_MM)
    wafer = ZernikeBasis(12)
    field_b = CartesianBasis(2)
    wtx, wty = {"Z4": 2.5, "Z5": -1.8}, {"Z7": 0.9}
    ftx, fty = {"1,0": 1.2}, {"0,1": -0.8}

    full = make_fingerprint(points, wafer, wtx, wty, field_b, ftx, fty)
    assert full.shape == (len(points), 2)

    only_wafer = make_fingerprint(points, wafer, wtx, wty, field_b, {}, {})
    only_field = make_fingerprint(points, wafer, {}, {}, field_b, ftx, fty)
    assert np.allclose(full, only_wafer + only_field)
    assert np.allclose(
        only_wafer[:, 0],
        make_fingerprint(points, wafer, {"Z4": 2.5}, {}, field_b, {}, {})[:, 0]
        + make_fingerprint(points, wafer, {"Z5": -1.8}, {}, field_b, {}, {})[:, 0],
    )


def test_fingerprint_empty_terms_give_zero():
    points = sample_points(_sampled(10), 4, R_WAFER_MM, EDGE_MM)
    out = make_fingerprint(
        points, ZernikeBasis(12), {}, {}, CartesianBasis(2), {}, {}
    )
    assert np.allclose(out, 0.0)


def test_fingerprint_piston_only_shifts_constant():
    """只给 x 的项 → dy 恒为 0（x/y 是两套独立系数，不是共享一套）。"""
    points = sample_points(_sampled(10), 4, R_WAFER_MM, EDGE_MM)
    out = make_fingerprint(points, ZernikeBasis(4), {"Z1": 3.0}, {}, CartesianBasis(2), {}, {})
    assert np.allclose(out[:, 0], 3.0)
    assert np.allclose(out[:, 1], 0.0)


def test_fingerprint_field_term_uses_normalized_field_coords():
    """field 项必须用归一化场内坐标 ξ/η（[-1,1]），不是 mm —— 用 ξ 项在 ±1 处检验。"""
    points = sample_points(_sampled(10), 4, R_WAFER_MM, EDGE_MM)
    out = make_fingerprint(
        points, ZernikeBasis(4), {}, {}, CartesianBasis(2), {"1,0": 1.0}, {}
    )
    xi, _ = points.xi_eta()
    assert np.allclose(out[:, 0], xi)
    assert np.abs(out[:, 0]).max() <= 1.0 + 1e-12


# -------------------------------------------------------- TIS 与噪声


def test_tis_roundtrip_without_noise():
    """σ=0 时，(m0−m180)/2 必须精确还原真值，(m0+m180)/2 精确等于 TIS。"""
    true = np.array([[3.0, -1.5], [-2.0, 4.0]])
    tis = (0.26, 0.15)
    meas = apply_tis(true, tis, calibrated=True, sigma_nm=0.0, rng=np.random.default_rng(1))
    combined, tis_hat = combine_tis(meas["m0"], meas["m180"])
    assert np.allclose(combined, true)
    assert np.allclose(tis_hat, np.array(tis))


def test_tis_uncalibrated_keeps_bias():
    true = np.array([[3.0, -1.5]])
    tis = (0.26, 0.15)
    meas = apply_tis(true, tis, calibrated=False, sigma_nm=0.0, rng=np.random.default_rng(1))
    assert set(meas) == {"m"}
    assert np.allclose(meas["m"] - true, np.array(tis))


def test_tis_combination_halves_the_noise_variance():
    """0°/180° 组合把噪声压到 σ/√2 ——「测两次」的收益就在这里，不是白测。

        m0   = +true + TIS + ε0
        m180 = −true + TIS + ε180
        (m0 − m180)/2 = true + (ε0 − ε180)/2     →  Var = (σ²+σ²)/4 = σ²/2

    所以 TIS 校准后的等效 TMU 是 σ_tmu/√2（默认 0.1nm → 0.0707nm）。
    这也是 A-03 里系数恢复精度的实际噪声水平。
    """
    true = np.zeros((4000, 2))
    tis = (0.0, 0.0)
    a = apply_tis(true, tis, calibrated=True, sigma_nm=0.1, rng=np.random.default_rng(SEED))
    b = apply_tis(true, tis, calibrated=True, sigma_nm=0.1, rng=np.random.default_rng(SEED))
    assert np.allclose(a["m0"], b["m0"]) and np.allclose(a["m180"], b["m180"])

    combined, _ = combine_tis(a["m0"], a["m180"])
    assert combined.std() == pytest.approx(0.1 / np.sqrt(2.0), rel=0.05)
    # 单次读数本身仍是 σ
    assert a["m0"].std() == pytest.approx(0.1, rel=0.05)


def test_tis_validation():
    with pytest.raises(ValueError, match=r"\(n_points, 2\)"):
        apply_tis(np.zeros(3), (0.0, 0.0), True, 0.1, np.random.default_rng(1))
    with pytest.raises(ValueError, match="形状不一致"):
        combine_tis(np.zeros((2, 2)), np.zeros((3, 2)))
    with pytest.raises(ValueError, match="TIS 幅值"):
        tis_vector_from_polar(-1.0, 0.0)


def test_tis_vector_from_polar_matches_documented_default():
    dx, dy = tis_vector_from_polar(0.3, 30.0)
    assert round(dx, 2) == 0.26 and round(dy, 2) == 0.15  # 详设 §6.3 的示例值


def test_add_noise_zero_sigma_and_reproducibility():
    true = np.full((10, 2), 5.0)
    assert np.allclose(add_noise(true, 0.0, np.random.default_rng(1)), true)
    a = add_noise(true, 0.5, np.random.default_rng(SEED))
    b = add_noise(true, 0.5, np.random.default_rng(SEED))
    assert np.allclose(a, b)
    assert not np.allclose(a, true)
    with pytest.raises(ValueError, match="sigma_nm"):
        add_noise(true, -0.1, np.random.default_rng(1))
