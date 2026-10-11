"""基函数接口。

不同基的自变量不同：
- wafer-level Zernike：`values(rho, theta, terms=...)`，自变量是归一化极坐标（单位圆）
- field-level Cartesian 幂基：`values(xi, eta, terms=...)`，自变量是单位方形坐标 [-1,1]^2

因此协议用 `*coords` 表达「基自己的自然坐标」，而不是强行统一成一个数组。
"""

from __future__ import annotations

from typing import Any, Protocol, runtime_checkable

import numpy as np
from numpy.typing import NDArray


@runtime_checkable
class Basis(Protocol):
    """设计矩阵的提供者。

    `terms()` 返回该项基的项名列表；`values(...)` 返回形状 (n_points, n_terms)
    的设计矩阵，列顺序与传入的 `terms` 严格一致。`terms=None` 表示取该项基的全部项。
    """

    name: str

    def terms(self) -> list[str]: ...

    def values(
        self, coord1: Any, coord2: Any, terms: list[str] | None = None
    ) -> NDArray[np.float64]: ...
