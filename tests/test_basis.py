"""basis 单元测试：Noll 索引、归一化、已知取值、设计矩阵与最小二乘里程碑（A-02）。"""

from __future__ import annotations

import math

import numpy as np
import pytest

from eumetria.basis import ZernikeBasis, noll_terms, noll_to_nm
from eumetria.geometry import build_layout, field_marks, to_polar

# Noll 索引的真值表（j → (n, m)）。这张表就是 §9.2 约定，测试是它的执行版本。
NOLL_TABLE = {
    1: (0, 0),
    2: (1, 1),
    3: (1, -1),
    4: (2, 0),
    5: (2, -2),
    6: (2, 2),
    7: (3, -1),
    8: (3, 1),
    9: (3, -3),
    10: (3, 3),
    11: (4, 0),
    12: (4, 2),
    13: (4, -2),
    14: (4, 4),
    15: (4, -4),
    16: (5, 1),
    17: (5, -1),
    18: (5, 3),
    19: (5, -3),
    20: (5, 5),
    21: (5, -5),
    22: (6, 0),
    23: (6, -2),
    24: (6, 2),
    25: (6, -4),
    26: (6, 4),
    27: (6, -6),
    28: (6, 6),
    29: (7, -1),
    30: (7, 1),
    31: (7, -3),
    32: (7, 3),
    33: (7, -5),
    34: (7, 5),
    35: (7, -7),
    36: (7, 7),
}


# ------------------------------------------------------------- Noll 索引


@pytest.mark.parametrize(("j", "nm"), sorted(NOLL_TABLE.items()))
def test_noll_to_nm_matches_convention(j, nm):
    assert noll_to_nm(j) == nm


def test_noll_convention_traps_the_ansi_confusion():
    """Noll 的球差是 j=11；j=9 是三叶草像差。写成 Z9=spherical 是 ANSI/Fringe 的编号。"""
    assert noll_to_nm(11) == (4, 0)  # spherical
    assert noll_to_nm(9) == (3, -3)  # trefoil, 不是 spherical
    # 常用像差在 Noll 下的位置（与文献核对过，可放心引用）
    assert noll_to_nm(4) == (2, 0)  # defocus
    assert {noll_to_nm(5), noll_to_nm(6)} == {(2, -2), (2, 2)}  # astigmatism
    assert {noll_to_nm(7), noll_to_nm(8)} == {(3, -1), (3, 1)}  # coma


def test_noll_to_nm_full_shells_are_permutations():
    """每个壳层恰好覆盖该 n 的全部合法 m（无遗漏、无重复）。"""
    seen: dict[int, set[int]] = {}
    for j in range(1, 37):
        n, m = noll_to_nm(j)
        assert abs(m) <= n and (n - abs(m)) % 2 == 0
        seen.setdefault(n, set()).add(m)
    for n, ms in seen.items():
        expected = {m for m in range(-n, n + 1) if (n - abs(m)) % 2 == 0}
        assert ms == expected, f"n={n} 的 m 集合不完整"


def test_noll_invalid_indices():
    for bad in (0, -3):
        with pytest.raises(ValueError, match="从 1 开始"):
            noll_to_nm(bad)
    with pytest.raises(ValueError, match="整数"):
        noll_to_nm(2.5)  # type: ignore[arg-type]


def test_terms_naming():
    assert noll_terms(4) == ["Z1", "Z2", "Z3", "Z4"]
    assert ZernikeBasis(6).terms() == noll_terms(6)


# -------------------------------------------------------- 径向与归一化


@pytest.mark.parametrize(
    ("n", "m"),
    [(0, 0), (2, 0), (3, 1), (4, 0), (4, 4), (6, 0), (7, 7), (8, 0)],
)
def test_radial_polynomial_is_one_at_rho_one(n, m):
    rho = np.array([1.0])
    assert ZernikeBasis.radial(n, m, rho)[0] == pytest.approx(1.0)


def test_radial_polynomial_known_closed_forms():
    rho = np.linspace(0, 1, 11)
    assert np.allclose(ZernikeBasis.radial(2, 0, rho), 2 * rho**2 - 1)
    assert np.allclose(ZernikeBasis.radial(3, 1, rho), 3 * rho**3 - 2 * rho)
    assert np.allclose(ZernikeBasis.radial(4, 0, rho), 6 * rho**4 - 6 * rho**2 + 1)
    assert np.allclose(ZernikeBasis.radial(2, 2, rho), rho**2)


def test_radial_rejects_illegal_nm():
    with pytest.raises(ValueError, match="非法"):
        ZernikeBasis.radial(2, 3, 0.5)  # |m| > n
    with pytest.raises(ValueError, match="非法"):
        ZernikeBasis.radial(3, 2, 0.5)  # n - |m| 为奇数


def test_normalization_matches_unit_variance_convention():
    assert ZernikeBasis.normalization(0, 0) == pytest.approx(1.0)
    assert ZernikeBasis.normalization(1, 1) == pytest.approx(2.0)
    assert ZernikeBasis.normalization(2, 0) == pytest.approx(math.sqrt(3.0))
    assert ZernikeBasis.normalization(4, 0) == pytest.approx(math.sqrt(5.0))
    assert ZernikeBasis.normalization(4, 2) == pytest.approx(math.sqrt(10.0))


def test_known_zernike_values():
    rho = np.array([0.0, 0.5, 1.0])
    theta = np.array([0.3, 1.2, 2.7])
    basis = ZernikeBasis(11)

    assert np.allclose(basis.values(rho, theta, ["Z1"])[:, 0], 1.0)
    assert np.allclose(basis.values(rho, theta, ["Z2"])[:, 0], 2 * rho * np.cos(theta))
    assert np.allclose(basis.values(rho, theta, ["Z3"])[:, 0], 2 * rho * np.sin(theta))
    assert np.allclose(
        basis.values(rho, theta, ["Z4"])[:, 0], math.sqrt(3) * (2 * rho**2 - 1)
    )
    assert np.allclose(
        basis.values(rho, theta, ["Z11"])[:, 0],
        math.sqrt(5) * (6 * rho**4 - 6 * rho**2 + 1),
    )


# ------------------------------------------------------------- 设计矩阵


def test_design_matrix_shape_and_column_order():
    basis = ZernikeBasis(12)
    rho = np.linspace(0.05, 1.0, 20)
    theta = np.linspace(-np.pi, np.pi, 20, endpoint=False)
    full = basis.values(rho, theta)
    assert full.shape == (20, 12)

    picked = ["Z4", "Z1", "Z11"]
    sub = basis.values(rho, theta, picked)
    assert sub.shape == (20, 3)
    for col, term in enumerate(picked):
        assert np.allclose(sub[:, col], basis.values(rho, theta, [term])[:, 0])


def test_design_matrix_input_validation():
    basis = ZernikeBasis(6)
    with pytest.raises(ValueError, match="形状不一致"):
        basis.values(np.zeros(3), np.zeros(4))
    with pytest.raises(ValueError, match="重复"):
        basis.values(np.array([0.5]), np.array([0.0]), ["Z1", "Z1"])
    with pytest.raises(ValueError, match="超出已初始化的范围"):
        basis.values(np.array([0.5]), np.array([0.0]), ["Z9"])
    with pytest.raises(ValueError, match="形式"):
        basis.values(np.array([0.5]), np.array([0.0]), ["P4"])


# --------------------------------------------------- 里程碑：系数可恢复性


def _disk_points(n_points: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    rho = np.sqrt(rng.uniform(0.01**2, 0.98**2, n_points))
    theta = rng.uniform(-np.pi, np.pi, n_points)
    return rho, theta


def test_milestone_lstsq_recovers_known_coefficients():
    """第一个真里程碑（检查表第 3 步）+ 验收条件 A-02。

    无噪声、点数 ≥ 8× 系数个数 → 全部系数恢复相对误差 < 1%。
    """
    basis = ZernikeBasis(max_order=12)
    terms = basis.terms()
    rho, theta = _disk_points(8 * len(terms), seed=20261003)
    rng = np.random.default_rng(7)
    truth = rng.normal(0.0, 3.0, size=len(terms))  # 量级 nm

    design = basis.values(rho, theta, terms)
    y = design @ truth  # 无噪声
    recovered, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)

    assert rank == len(terms)
    rel_err = np.abs(recovered - truth) / np.maximum(np.abs(truth), 1e-12)
    assert rel_err.max() < 0.01, f"最大相对误差 {rel_err.max():.4g} 超过 1%"

    cond = float(np.linalg.cond(design))
    assert cond < 1e6, f"条件数 {cond:.3g} 超过 §9.3 的告警阈值"


def test_milestone_on_real_layout_marks():
    """用真实布局的 mark 点集跑一遍：geometry → basis → 最小二乘。

    顺带钉住一个必须记住的事实：口径 B 纳入的 partial field 会让一部分 mark 落到
    **片外**（本轮实测 rho max = 1.061，即 r ≈ 159mm > 150mm）。所以「按 edge exclusion
    过滤 mark」不是可选优化，而是采样层（§7.3）的必需步骤 —— 不过滤就直接喂进 Zernike
    会触发定义域报错，这正是 §9.2 要求 ρ>1 抛错而不是静默截断的意义。
    """
    fields = build_layout(300.0, 3.0, 26.0, 33.0, include_partial=True)
    xs: list[float] = []
    ys: list[float] = []
    for field in fields:
        for mark in field_marks(field, 4):
            xs.append(mark.x_mm)
            ys.append(mark.y_mm)
    rho, theta = to_polar(np.array(xs), np.array(ys), 150.0)

    assert rho.max() > 1.0, "口径 B 应当产生片外 mark；若不再成立，采样层的必要性论证要重写"
    basis = ZernikeBasis(max_order=8)
    terms = basis.terms()

    usable_rho = (150.0 - 3.0) / 150.0
    keep = rho <= usable_rho
    assert keep.sum() >= 8 * len(terms)
    with pytest.raises(ValueError, match="单位圆"):
        basis.values(rho, theta, terms)  # 未过滤 → 定义域报错

    rng = np.random.default_rng(11)
    truth = rng.normal(0.0, 2.0, size=len(terms))
    design = basis.values(rho[keep], theta[keep], terms)
    recovered, _, _, _ = np.linalg.lstsq(design, design @ truth, rcond=None)
    assert np.allclose(recovered, truth, rtol=1e-6, atol=1e-8)
