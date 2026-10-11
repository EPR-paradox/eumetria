"""Zernike 基（Noll 索引，单位方差 / 正交归一化约定）。

约定声明（§9.2）—— 引用本模块前必须先读这一段
==============================================

索引：**Noll**。下面这张表是唯一的真值来源，任何「Z9 是球差」的说法都是错的
（那是 ANSI / Fringe / Wyant 的 1-based 编号，不是 Noll）：

    j= 1 (n=0,m= 0) piston        j= 2 (1, 1) tilt-x      j= 3 (1,-1) tilt-y
    j= 4 (2, 0)      defocus      j= 5 (2,-2) astig 45°   j= 6 (2, 2) astig 0°
    j= 7 (3,-1)      coma-y       j= 8 (3, 1) coma-x      j= 9 (3,-3) trefoil-y
    j=10 (3, 3)      trefoil-x    j=11 (4, 0) **spherical**（不是 Z9！）
    j=12 (4, 2) / j=13 (4,-2) / j=14 (4,4) / j=15 (4,-4)

Noll 与 ANSI 的错位正是从 j=9 开始的：Noll 的球差是 j=11，ANSI/Fringe 的球差是 9。
本项目全程用 Noll，并在 `terms()` 里用 "Z<Noll>" 命名，避免歧义。

归一化：单位方差，
    N_n^m = sqrt(2(n+1) / (1 + δ_{m0}))
    角向：m > 0 → cos(mθ)，m < 0 → sin(|m|θ)，m = 0 → 1
    于是 (1/π) ∬ Z_j Z_k ρ dρ dθ = δ_jk（单位圆上正交且单位方差）。
    `zernpy` 用的是同一套归一化，本仓库的性质测试直接对它数值对照。

实现：`noll_to_nm` 按 Noll 的壳层规则程序化生成（不抄索引表），径向多项式用闭式
求和（不抄系数表），θ 属于 [-π, π)，定义域 ρ ∈ [0, 1]。
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
from numpy.typing import NDArray

ArrLike = float | NDArray[np.float64] | list[float] | tuple[float, ...]

# 浮点容差：仅用于「ρ 是否越界」的判定，不用于数值截断
_RHO_TOL = 1e-12


class NollConventionError(ValueError):
    """Noll 索引非法（j < 1 或非整数）。"""


def _shell_orders(n: int) -> list[int]:
    """第 n 壳层内按 Noll 顺序排列的 m 序列（不含符号表的硬编码，按规则生成）。

    规则：壳层内的起始符号 s = (-1)^(n // 2)，每对 (m, -m) 共享该符号；
    n 为偶数时 |m| 取 0,2,4,...,n；n 为奇数时取 1,3,5,...,n。
    """
    s = -1 if (n // 2) % 2 else 1
    if n % 2 == 0:
        abs_ms = [0] + list(range(2, n + 1, 2))
    else:
        abs_ms = list(range(1, n + 1, 2))

    out: list[int] = []
    for m_abs in abs_ms:
        if m_abs == 0:
            out.append(0)
        else:
            out.append(s * m_abs)
            out.append(-s * m_abs)
    return out


def noll_to_nm(j: int) -> tuple[int, int]:
    """Noll 索引 j（从 1 开始）→ (n, m)。"""
    if not isinstance(j, (int, np.integer)) or isinstance(j, bool):
        raise NollConventionError(f"Noll 索引必须是整数，得到 {j!r}")
    j = int(j)
    if j < 1:
        raise NollConventionError(f"Noll 索引从 1 开始，得到 {j}")

    count = 0
    n = 0
    while True:
        shell = _shell_orders(n)
        if count + len(shell) >= j:
            return n, shell[j - count - 1]
        count += len(shell)
        n += 1


def noll_count(max_order: int) -> int:
    """j ≤ max_order 的项数（= max_order，保留函数是为了让调用方不必假设恒等）。"""
    _check_max_order(max_order)
    return int(max_order)


def noll_terms(max_order: int) -> list[str]:
    """返回 ["Z1", "Z2", ..., "Z<max_order>"]。"""
    return [f"Z{j}" for j in range(1, noll_count(max_order) + 1)]


def _check_max_order(max_order: int) -> None:
    if not isinstance(max_order, (int, np.integer)) or isinstance(max_order, bool):
        raise NollConventionError(f"max_order 必须是整数，得到 {max_order!r}")
    if int(max_order) < 1:
        raise NollConventionError(f"max_order 必须 ≥ 1，得到 {max_order}")
    if int(max_order) > 1539:
        raise NollConventionError(f"max_order={max_order} 超出可计算范围（≤1539）")


def _as_float_array(value: ArrLike) -> NDArray[np.float64]:
    return np.asarray(value, dtype=np.float64)


class ZernikeBasis:
    """Noll 索引的 Zernike 基，单位方差归一化。

    `max_order` 是**最大 Noll 索引 j**（§9.4：max_order=12 → Z1..Z12），
    不是最大径向阶数 n。
    """

    name = "zernike"

    def __init__(self, max_order: int = 12) -> None:
        _check_max_order(max_order)
        self.max_order = int(max_order)
        self._terms = noll_terms(self.max_order)

    # ------------------------------------------------------------- 元信息

    def terms(self) -> list[str]:
        return list(self._terms)

    def __repr__(self) -> str:  # pragma: no cover - 便于调试
        return f"ZernikeBasis(max_order={self.max_order}, convention='noll')"

    # --------------------------------------------------------- 索引与径向

    @staticmethod
    def noll_to_nm(j: int) -> tuple[int, int]:
        return noll_to_nm(j)

    @staticmethod
    def radial(n: int, m: int, rho: ArrLike) -> NDArray[np.float64]:
        """径向多项式 R_n^{|m|}(ρ)，闭式求和（整数系数，无硬编码系数表）。"""
        m_abs = abs(int(m))
        n = int(n)
        if n < 0 or m_abs > n or (n - m_abs) % 2 != 0:
            raise ValueError(f"非法 (n, m) = ({n}, {m})：需要 n ≥ |m| 且 (n - |m|) 为偶数")

        rho_arr = _as_float_array(rho)
        kmax = (n - m_abs) // 2
        result = np.zeros_like(rho_arr, dtype=np.float64)
        for k in range(kmax + 1):
            num = math.factorial(n - k)
            den = (
                math.factorial(k)
                * math.factorial((n + m_abs) // 2 - k)
                * math.factorial((n - m_abs) // 2 - k)
            )
            coeff = num // den if (num % den == 0) else num / den
            if k % 2:
                result -= coeff * np.power(rho_arr, n - 2 * k)
            else:
                result += coeff * np.power(rho_arr, n - 2 * k)
        return result

    @staticmethod
    def normalization(n: int, m: int) -> float:
        """N_n^m = sqrt(2(n+1) / (1 + δ_{m0}))。"""
        return math.sqrt(2.0 * (int(n) + 1) / (1.0 + (1.0 if int(m) == 0 else 0.0)))

    # --------------------------------------------------------- 设计矩阵

    def values(
        self,
        rho: ArrLike,
        theta: ArrLike,
        terms: list[str] | None = None,
    ) -> NDArray[np.float64]:
        """设计矩阵，形状 (n_points, n_terms)。

        ρ 必须落在 [0, 1]（单位圆）；越界抛 ValueError，不静默截断（§9.2）。
        θ 属于 [-π, π)（任何实数都可，周期函数），单位为弧度。
        """
        rho_arr = _as_float_array(rho)
        theta_arr = _as_float_array(theta)
        if rho_arr.shape != theta_arr.shape:
            raise ValueError(f"rho 与 theta 形状不一致：{rho_arr.shape} vs {theta_arr.shape}")

        self._check_domain(rho_arr)
        term_list = self._terms if terms is None else list(terms)
        if len(set(term_list)) != len(term_list):
            raise ValueError(f"terms 中存在重复项：{term_list}")

        out = np.empty((rho_arr.size, len(term_list)), dtype=np.float64)
        for col, term in enumerate(term_list):
            j = self._parse_term(term)
            n, m = noll_to_nm(j)
            radial = self.radial(n, m, rho_arr.ravel())
            if m > 0:
                angular = np.cos(m * theta_arr.ravel())
            elif m < 0:
                angular = np.sin(abs(m) * theta_arr.ravel())
            else:
                angular = np.ones(rho_arr.size, dtype=np.float64)
            out[:, col] = self.normalization(n, m) * radial * angular
        return out

    # ------------------------------------------------------------- 内部

    @staticmethod
    def _check_domain(rho_arr: NDArray[np.float64]) -> None:
        if not np.all(np.isfinite(rho_arr)):
            raise ValueError("rho 含非有限值")
        if np.any(rho_arr < -_RHO_TOL) or np.any(rho_arr > 1.0 + _RHO_TOL):
            worst = float(np.nanmax(np.abs(rho_arr)))
            raise ValueError(
                f"rho 超出单位圆定义域 [0,1]（最大 |rho| = {worst:.6f}）；"
                "越界点应在采样阶段排除，拟合不静默截断"
            )

    def _parse_term(self, term: str) -> int:
        if not isinstance(term, str) or not term.startswith("Z"):
            raise ValueError(f"Zernike 项名必须是 'Z<j>' 形式，得到 {term!r}")
        try:
            j = int(term[1:])
        except ValueError as exc:
            raise ValueError(f"Zernike 项名必须是 'Z<j>' 形式，得到 {term!r}") from exc
        if not 1 <= j <= self.max_order:
            raise ValueError(f"项 {term!r} 超出已初始化的范围 Z1..Z{self.max_order}")
        return j

    def subset(self, terms: Iterable[str]) -> list[str]:
        """校验一组项名是否合法并原样返回（顺序保留）。"""
        return [t for t in terms if self._parse_term(t) is not None]
