"""硅片 field 布局与场内 mark 位置。

设计约束（§3.6）：
- field 是规则网格、中心对齐；候选网格尺寸由硅片直径决定，不由可用区决定。
- 两种口径：B（默认，field 中心在可用圆内）与 A（四角都在可用圆内）。
  默认场景 300mm / 26×33mm / 3mm 排除 → B = 76，A = 64。
- 全部确定性，无随机。
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# 4 点 mark 布局距 field 边缘的内缩量（mm）。is_edge 的判定依赖它（见 build_layout）。
MARK_INSET_MM = 3.0
# 已定义且无歧义的 mark 密度。19 点布局在详设 §7.1 里的定义自相矛盾（见 U-10），不实现。
SUPPORTED_MARKS_PER_FIELD = (4, 8)


@dataclass(frozen=True)
class Field:
    """一个曝光场。cx/cy 为全场坐标下的 field 中心（mm）。"""

    field_id: int
    cx_mm: float
    cy_mm: float
    width_mm: float
    height_mm: float
    is_partial: bool
    is_edge: bool

    def corners_mm(self) -> tuple[tuple[float, float], ...]:
        """四个角的全场坐标。"""
        hw, hh = self.width_mm / 2.0, self.height_mm / 2.0
        return (
            (self.cx_mm - hw, self.cy_mm - hh),
            (self.cx_mm - hw, self.cy_mm + hh),
            (self.cx_mm + hw, self.cy_mm - hh),
            (self.cx_mm + hw, self.cy_mm + hh),
        )


@dataclass(frozen=True)
class Mark:
    """一个量测标记。

    `xf_mm` / `yf_mm` 是相对 field 中心的场内坐标；`x_mm` / `y_mm` 是对应的全场坐标。
    `mark_id` 只在所属 field 内唯一 —— 全局编号由采样层（synth.sampling）负责，
    因为 `field_marks` 的签名里没有全局上下文。
    """

    mark_id: int
    xf_mm: float
    yf_mm: float
    x_mm: float
    y_mm: float


def build_layout(
    wafer_diameter_mm: float,
    edge_exclusion_mm: float,
    field_width_mm: float,
    field_height_mm: float,
    include_partial: bool = True,
) -> list[Field]:
    """规则网格铺满候选区（网格 = ceil(2R/W) × ceil(2R/H)），中心对齐。

    include_partial=True （默认，口径 B）：field 中心在可用圆内即纳入，标 is_partial。
    include_partial=False（口径 A）：仅四角都在可用圆内的 field。

    `is_edge` 的含义是「该 field 的默认 4 点 mark 布局中，至少一个 mark 落在边缘
    排除区内（r > R - edge_exclusion）」—— 它比 is_partial 更严：一个被裁切的 field
    可能因为内缩 3mm 而让所有 mark 仍然落在可用区内。`is_partial` 描述几何，
    `is_edge` 描述量测点，两者不等价。
    """
    radius_mm = wafer_diameter_mm / 2.0
    if not math.isfinite(radius_mm) or radius_mm <= 0:
        raise ValueError(f"wafer_diameter_mm 必须为正的有限值，得到 {wafer_diameter_mm!r}")
    if not math.isfinite(edge_exclusion_mm) or not 0 <= edge_exclusion_mm < radius_mm:
        raise ValueError(
            f"edge_exclusion_mm 必须在 [0, R) 内，得到 {edge_exclusion_mm!r}（R={radius_mm!r}）"
        )
    if not math.isfinite(field_width_mm) or field_width_mm <= 0:
        raise ValueError(f"field_width_mm 必须为正的有限值，得到 {field_width_mm!r}")
    if not math.isfinite(field_height_mm) or field_height_mm <= 0:
        raise ValueError(f"field_height_mm 必须为正的有限值，得到 {field_height_mm!r}")
    if min(field_width_mm, field_height_mm) <= 2.0 * MARK_INSET_MM:
        raise ValueError(
            f"field 尺寸 {field_width_mm}×{field_height_mm} 小于 mark 内缩 {MARK_INSET_MM}mm 的要求"
        )

    usable_radius_mm = radius_mm - edge_exclusion_mm
    nx = math.ceil(2.0 * radius_mm / field_width_mm)
    ny = math.ceil(2.0 * radius_mm / field_height_mm)
    mark_dx = field_width_mm / 2.0 - MARK_INSET_MM
    mark_dy = field_height_mm / 2.0 - MARK_INSET_MM

    fields: list[Field] = []
    for j in range(ny):
        cy_mm = (j - (ny - 1) / 2.0) * field_height_mm
        for i in range(nx):
            cx_mm = (i - (nx - 1) / 2.0) * field_width_mm
            r_center_mm = math.hypot(cx_mm, cy_mm)
            corners_inside = all(
                math.hypot(cx, cy) <= usable_radius_mm
                for cx, cy in (
                    (cx_mm - field_width_mm / 2.0, cy_mm - field_height_mm / 2.0),
                    (cx_mm - field_width_mm / 2.0, cy_mm + field_height_mm / 2.0),
                    (cx_mm + field_width_mm / 2.0, cy_mm - field_height_mm / 2.0),
                    (cx_mm + field_width_mm / 2.0, cy_mm + field_height_mm / 2.0),
                )
            )

            if include_partial:
                if r_center_mm > usable_radius_mm:
                    continue
            elif not corners_inside:
                continue

            is_edge = any(
                math.hypot(cx_mm + sx * mark_dx, cy_mm + sy * mark_dy) > usable_radius_mm
                for sx in (-1.0, 1.0)
                for sy in (-1.0, 1.0)
            )
            fields.append(
                Field(
                    field_id=len(fields),
                    cx_mm=cx_mm,
                    cy_mm=cy_mm,
                    width_mm=field_width_mm,
                    height_mm=field_height_mm,
                    is_partial=not corners_inside,
                    is_edge=is_edge,
                )
            )
    return fields


def field_marks(field: Field, marks_per_field: int) -> list[Mark]:
    """返回场内的 mark 位置（场内坐标 + 全场坐标）。

    4  → 四角内侧固定偏移（距边 MARK_INSET_MM）
    8  → 4 角 + 4 边中点（同一内缩量）
    19 → 详设 §7.1 的定义自相矛盾（25 - 6 既非 90° 对称也非 180° 对称），
         列为待定项 U-10，此处显式拒绝，不静默降级。
    """
    if marks_per_field not in SUPPORTED_MARKS_PER_FIELD:
        raise ValueError(
            f"不支持的 marks_per_field={marks_per_field!r}；"
            f"已定义的是 {SUPPORTED_MARKS_PER_FIELD}（19 点布局待定，见详设 U-10）"
        )

    hw = field.width_mm / 2.0 - MARK_INSET_MM
    hh = field.height_mm / 2.0 - MARK_INSET_MM
    if hw <= 0 or hh <= 0:
        raise ValueError(
            f"field {field.width_mm}×{field.height_mm}mm 去掉 {MARK_INSET_MM}mm 内缩后没有可用 mark 区"
        )

    local: list[tuple[float, float]] = [
        (-hw, -hh),
        (-hw, hh),
        (hw, -hh),
        (hw, hh),
    ]
    if marks_per_field == 8:
        local += [(-hw, 0.0), (hw, 0.0), (0.0, -hh), (0.0, hh)]

    return [
        Mark(
            mark_id=k,
            xf_mm=xf,
            yf_mm=yf,
            x_mm=field.cx_mm + xf,
            y_mm=field.cy_mm + yf,
        )
        for k, (xf, yf) in enumerate(local)
    ]
