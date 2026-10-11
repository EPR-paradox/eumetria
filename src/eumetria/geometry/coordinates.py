"""坐标换算：全场坐标 → wafer 归一化极坐标；场内坐标 → 单位方形坐标。

约定（§9.1）：坐标 mm，角度 rad。
`to_polar` 的分母是**物理硅片半径** R_wafer，不是排除后的可用区半径（§6.2）——
否则 rho 的归一化基准会跟着采样方案变，Zernike 系数失去跨片可比性。
"""

from __future__ import annotations

from typing import Union

import numpy as np
from numpy.typing import NDArray

ArrLike = Union[float, NDArray[np.float64], list, tuple]


def to_polar(
    x_mm: ArrLike,
    y_mm: ArrLike,
    wafer_radius_mm: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """返回 (rho, theta)。

    rho = r / R_wafer ∈ [0, 1]（点若在片外则 > 1，本函数不截断、不报错，
    是否允许越界由调用方决定）；theta ∈ [-π, π)。
    """
    x = np.asarray(x_mm, dtype=np.float64)
    y = np.asarray(y_mm, dtype=np.float64)
    if x.shape != y.shape:
        raise ValueError(f"x_mm 与 y_mm 形状不一致：{x.shape} vs {y.shape}")
    if not np.isfinite(wafer_radius_mm) or wafer_radius_mm <= 0:
        raise ValueError(f"wafer_radius_mm 必须为正的有限值，得到 {wafer_radius_mm!r}")

    rho = np.hypot(x, y) / float(wafer_radius_mm)

    theta = np.arctan2(y, x)
    # arctan2 在负实轴上返回 +π；按契约映射到 -π，保证半开区间 [-π, π)
    theta = np.where(theta >= np.pi, -np.pi, theta)
    return rho, theta


def normalize_field(
    xf_mm: ArrLike,
    yf_mm: ArrLike,
    field_width_mm: float,
    field_height_mm: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """返回 (xi, eta) ∈ [-1, 1]^2：场内坐标相对于 field 半宽/半高归一化。

    field-level 基函数（Cartesian 幂基 / Legendre）的自变量就是这个坐标；
    用 mm 直接当自变量会让不同 field 尺寸的场景不可比。
    """
    xf = np.asarray(xf_mm, dtype=np.float64)
    yf = np.asarray(yf_mm, dtype=np.float64)
    if xf.shape != yf.shape:
        raise ValueError(f"xf_mm 与 yf_mm 形状不一致：{xf.shape} vs {yf.shape}")
    if not np.isfinite(field_width_mm) or field_width_mm <= 0:
        raise ValueError(f"field_width_mm 必须为正的有限值，得到 {field_width_mm!r}")
    if not np.isfinite(field_height_mm) or field_height_mm <= 0:
        raise ValueError(f"field_height_mm 必须为正的有限值，得到 {field_height_mm!r}")

    xi = 2.0 * xf / float(field_width_mm)
    eta = 2.0 * yf / float(field_height_mm)
    return xi, eta
