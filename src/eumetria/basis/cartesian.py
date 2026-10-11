"""Cartesian 幂基（field-level）：ξ^p·η^q，(p+q) ≤ order，定义域 [-1,1]^2。

项名形如 "p,q"。枚举顺序按总阶数升序、同阶内 p 降序，于是 order=2 时是

    (0,0) (1,0) (0,1) (2,0) (1,1) (0,2)      —— 共 6 项

与详设 §6.2.2「field 2 阶 = 2×6 = 12 项」的口径一致。

为什么 field-level 用幂基而不用 Zernike（§3.1）：
扫描机的场内校正能力本身就是按 ξ/η 的多项式定义的，用同一套基可以直接和
扫描机系数对话；而 wafer-level 用 Zernike 是因为它在圆上正交（系数稳定）。
两层的自变量不同（单位方形 vs 单位圆），不可互换。
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

ArrLike = float | NDArray[np.float64] | list[float] | tuple[float, ...]


class CartesianBasis:
    """幂基 ξ^p η^q，(p+q) ≤ order。"""

    name = "cartesian"

    def __init__(self, order: int = 2) -> None:
        if not isinstance(order, (int, np.integer)) or isinstance(order, bool) or order < 0:
            raise ValueError(f"order 必须是非负整数，得到 {order!r}")
        self.order = int(order)
        # 总阶数升序，同阶内 p 降序
        self._terms = [
            f"{p},{total - p}" for total in range(self.order + 1) for p in range(total, -1, -1)
        ]

    def terms(self) -> list[str]:
        return list(self._terms)

    def __repr__(self) -> str:  # pragma: no cover - 便于调试
        return f"CartesianBasis(order={self.order})"

    @staticmethod
    def parse_term(term: str) -> tuple[int, int]:
        parts = str(term).split(",")
        if len(parts) != 2:
            raise ValueError(f"Cartesian 项名必须是 'p,q' 形式，得到 {term!r}")
        try:
            p, q = int(parts[0]), int(parts[1])
        except ValueError as exc:
            raise ValueError(f"Cartesian 项名必须是 'p,q' 形式，得到 {term!r}") from exc
        if p < 0 or q < 0:
            raise ValueError(f"幂次必须非负，得到 {term!r}")
        return p, q

    def values(
        self,
        xi: ArrLike,
        eta: ArrLike,
        terms: list[str] | None = None,
    ) -> NDArray[np.float64]:
        """设计矩阵，形状 (n_points, n_terms)。列顺序与 terms 严格一致。"""
        xi_arr = np.asarray(xi, dtype=np.float64)
        eta_arr = np.asarray(eta, dtype=np.float64)
        if xi_arr.shape != eta_arr.shape:
            raise ValueError(f"xi 与 eta 形状不一致：{xi_arr.shape} vs {eta_arr.shape}")

        term_list = self._terms if terms is None else list(terms)
        if len(set(term_list)) != len(term_list):
            raise ValueError(f"terms 中存在重复项：{term_list}")

        out = np.empty((xi_arr.size, len(term_list)), dtype=np.float64)
        for col, term in enumerate(term_list):
            p, q = self.parse_term(term)
            if p + q > self.order:
                raise ValueError(f"项 {term!r} 的总阶数 {p + q} 超过已初始化的 order={self.order}")
            out[:, col] = np.power(xi_arr.ravel(), p) * np.power(eta_arr.ravel(), q)
        return out
