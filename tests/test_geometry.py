"""几何层测试：坐标换算 + field 布局（详设 §7.1、验收无直接对应项，属 Phase 1 前置）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from eumetria.geometry import (
    MARK_INSET_MM,
    build_layout,
    field_marks,
    normalize_field,
    to_polar,
)

# 默认场景（P-07）
WAFER_MM = 300.0
R_WAFER_MM = WAFER_MM / 2.0
EDGE_MM = 3.0
FIELD_W_MM = 26.0
FIELD_H_MM = 33.0


# ---------------------------------------------------------------- to_polar


def test_to_polar_normalization_uses_wafer_radius():
    rho, theta = to_polar([R_WAFER_MM], [0.0], R_WAFER_MM)
    assert rho[0] == pytest.approx(1.0)
    assert theta[0] == pytest.approx(0.0)

    # 147mm（可用区半径）不等于 1 —— 归一化基准是物理半径 150mm
    rho_usable, _ = to_polar([R_WAFER_MM - EDGE_MM], [0.0], R_WAFER_MM)
    assert rho_usable[0] == pytest.approx(147.0 / 150.0)


def test_to_polar_theta_is_half_open_at_pi():
    _rho, theta = to_polar([-1.0, 0.0, -1.0], [0.0, 1.0, -1e-300], R_WAFER_MM)
    assert theta[0] == pytest.approx(-math.pi)  # 负实轴 → -π，不是 +π
    assert theta[1] == pytest.approx(math.pi / 2)
    assert theta[2] == pytest.approx(-math.pi)
    assert np.all(theta >= -math.pi) and np.all(theta < math.pi)


def test_to_polar_allows_outside_wafer_without_clipping():
    rho, _ = to_polar([200.0], [0.0], R_WAFER_MM)
    assert rho[0] == pytest.approx(200.0 / 150.0)  # 不静默截断到 1


def test_to_polar_rejects_bad_input():
    with pytest.raises(ValueError, match="形状不一致"):
        to_polar([1.0, 2.0], [1.0], R_WAFER_MM)
    for bad_radius in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="wafer_radius_mm"):
            to_polar([0.0], [0.0], bad_radius)


def test_to_polar_scalar_input():
    rho, theta = to_polar(3.0, 4.0, 10.0)
    assert float(rho) == pytest.approx(0.5)
    assert float(theta) == pytest.approx(math.atan2(4.0, 3.0))


# ---------------------------------------------------------- normalize_field


def test_normalize_field_maps_corners_to_unit_square():
    xi, eta = normalize_field(
        [-FIELD_W_MM / 2, FIELD_W_MM / 2, 0.0],
        [-FIELD_H_MM / 2, FIELD_H_MM / 2, 0.0],
        FIELD_W_MM,
        FIELD_H_MM,
    )
    assert xi[0] == pytest.approx(-1.0) and eta[0] == pytest.approx(-1.0)
    assert xi[1] == pytest.approx(1.0) and eta[1] == pytest.approx(1.0)
    assert xi[2] == pytest.approx(0.0) and eta[2] == pytest.approx(0.0)


def test_normalize_field_rejects_bad_input():
    with pytest.raises(ValueError, match="形状不一致"):
        normalize_field([0.0], [0.0, 1.0], FIELD_W_MM, FIELD_H_MM)
    with pytest.raises(ValueError, match="field_width_mm"):
        normalize_field([0.0], [0.0], 0.0, FIELD_H_MM)
    with pytest.raises(ValueError, match="field_height_mm"):
        normalize_field([0.0], [0.0], FIELD_W_MM, -1.0)


# ------------------------------------------------------------- build_layout


def test_build_layout_default_scenario_counts():
    """默认场景 300mm / 26×33 / 3mm：口径 B = 76，口径 A = 64。"""
    partial = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=True)
    strict = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=False)
    assert len(partial) == 76
    assert len(strict) == 64
    # 口径 A 是口径 B 的子集
    strict_centers = {(f.cx_mm, f.cy_mm) for f in strict}
    partial_centers = {(f.cx_mm, f.cy_mm) for f in partial}
    assert strict_centers <= partial_centers


def test_build_layout_is_deterministic():
    a = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM)
    b = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM)
    assert a == b


def test_build_layout_field_centers_inside_usable_circle():
    usable = R_WAFER_MM - EDGE_MM
    for f in build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=True):
        assert math.hypot(f.cx_mm, f.cy_mm) <= usable


def test_build_layout_flags_are_consistent():
    partial = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=True)
    # 口径 B 下必然存在被裁切的 field；口径 A 下不存在
    assert any(f.is_partial for f in partial)
    strict = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=False)
    assert not any(f.is_partial for f in strict)
    # is_edge 比 is_partial 更严：被裁切不等于 mark 落进排除区
    assert all(f.is_edge for f in partial if f.is_edge)  # 平凡真，保证字段是 bool
    assert sum(f.is_edge for f in partial) <= sum(f.is_partial for f in partial)


def test_build_layout_field_ids_are_sequential():
    for fields in (
        build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM),
        build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM, include_partial=False),
    ):
        assert [f.field_id for f in fields] == list(range(len(fields)))


def test_build_layout_rejects_bad_input():
    with pytest.raises(ValueError, match="wafer_diameter_mm"):
        build_layout(0.0, EDGE_MM, FIELD_W_MM, FIELD_H_MM)
    with pytest.raises(ValueError, match="edge_exclusion_mm"):
        build_layout(WAFER_MM, R_WAFER_MM, FIELD_W_MM, FIELD_H_MM)
    with pytest.raises(ValueError, match="field_width_mm"):
        build_layout(WAFER_MM, EDGE_MM, 0.0, FIELD_H_MM)
    with pytest.raises(ValueError, match="field_height_mm"):
        build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, 0.0)


# --------------------------------------------------------------- field_marks


def test_field_marks_4_and_8_geometry():
    field = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM)[0]

    marks4 = field_marks(field, 4)
    assert [m.mark_id for m in marks4] == [0, 1, 2, 3]
    dx = FIELD_W_MM / 2 - MARK_INSET_MM
    dy = FIELD_H_MM / 2 - MARK_INSET_MM
    assert {(m.xf_mm, m.yf_mm) for m in marks4} == {
        (-dx, -dy),
        (-dx, dy),
        (dx, -dy),
        (dx, dy),
    }
    # 全场坐标 = field 中心 + 场内坐标
    for m in marks4:
        assert m.x_mm == pytest.approx(field.cx_mm + m.xf_mm)
        assert m.y_mm == pytest.approx(field.cy_mm + m.yf_mm)
    # 所有 mark 都在 field 内
    for m in marks4:
        assert abs(m.xf_mm) < FIELD_W_MM / 2
        assert abs(m.yf_mm) < FIELD_H_MM / 2

    marks8 = field_marks(field, 8)
    assert len(marks8) == 8
    assert {(m.xf_mm, m.yf_mm) for m in marks4} <= {(m.xf_mm, m.yf_mm) for m in marks8}
    assert (0.0, dy) in {(m.xf_mm, m.yf_mm) for m in marks8}
    assert (dx, 0.0) in {(m.xf_mm, m.yf_mm) for m in marks8}


def test_field_marks_rejects_undefined_density():
    field = build_layout(WAFER_MM, EDGE_MM, FIELD_W_MM, FIELD_H_MM)[0]
    with pytest.raises(ValueError, match="U-10"):
        field_marks(field, 19)
    with pytest.raises(ValueError, match="不支持的"):
        field_marks(field, 5)
