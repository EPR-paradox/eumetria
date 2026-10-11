"""Zernike 基的性质测试（§13.2）：正交性、单位方差、定义域处理。"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from eumetria.basis import ZernikeBasis

M_ORDER = 7  # 只测到 n ≤ 7 的壳层，够覆盖 astigmatism / coma / spherical 以上
N_RADIAL = 80
N_ANGULAR = 256


def _quadrature_weights() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """单位圆的乘积型求积：(rho_i, rho_w_i, theta_j, theta_w_j)。

    ∫∫ f ρ dρ dθ ≈ Σ_i Σ_j ρ_w_i θ_w_j f(ρ_i, θ_j)；归一化时除以 π。
    """
    x, w = leggauss(N_RADIAL)
    rho = 0.5 * (x + 1.0)
    rho_w = 0.5 * w * rho  # 含 Jacobian ρ
    theta = np.linspace(0.0, 2.0 * np.pi, N_ANGULAR, endpoint=False)
    theta_w = np.full(N_ANGULAR, 2.0 * np.pi / N_ANGULAR)
    return rho, rho_w, theta, theta_w


def _inner_product(basis: ZernikeBasis, ti: str, tj: str) -> float:
    rho, rho_w, theta, theta_w = _quadrature_weights()
    rho_grid, theta_grid = np.meshgrid(rho, theta, indexing="ij")
    vi = basis.values(rho_grid.ravel(), theta_grid.ravel(), [ti])[:, 0]
    vj = basis.values(rho_grid.ravel(), theta_grid.ravel(), [tj])[:, 0]
    grid = (vi * vj).reshape(rho.size, theta.size)
    return float(np.einsum("i,j,ij->", rho_w, theta_w, grid) / np.pi)


def test_zernike_basis_is_orthonormal_on_unit_disk():
    """(1/π) ∬ Z_i Z_j ρ dρ dθ = δ_ij，相对容差 1e-10（§13.2）。"""
    basis = ZernikeBasis(max_order=28)  # 覆盖 n ≤ 6 的全部壳层
    terms = basis.terms()
    worst = 0.0
    for i, ti in enumerate(terms):
        for tj in terms[i:]:
            value = _inner_product(basis, ti, tj)
            expected = 1.0 if ti == tj else 0.0
            worst = max(worst, abs(value - expected))
    assert worst < 1e-10, f"最大偏离 {worst:.3e} 超过 1e-10"


def test_zernike_unit_variance_of_individual_terms():
    basis = ZernikeBasis(max_order=15)
    for term in basis.terms():
        assert _inner_product(basis, term, term) == pytest.approx(1.0, abs=1e-10)


def test_domain_outside_unit_circle_raises():
    """ρ > 1 必须报错，不静默截断（§9.2）。"""
    basis = ZernikeBasis(4)
    with pytest.raises(ValueError, match="单位圆"):
        basis.values(np.array([1.0001]), np.array([0.3]))
    with pytest.raises(ValueError, match="单位圆"):
        basis.values(np.array([0.5, 2.0]), np.array([0.3, 0.3]))
    with pytest.raises(ValueError, match="非有限"):
        basis.values(np.array([np.nan]), np.array([0.3]))


def test_domain_boundary_is_inclusive():
    basis = ZernikeBasis(4)
    values = basis.values(np.array([0.0, 1.0]), np.array([0.0, 0.0]))
    assert values.shape == (2, 4)
    assert np.all(np.isfinite(values))


def test_angular_symmetry_under_theta_shift_by_2pi():
    """θ 是周期量：θ 与 θ+2π 必须给出相同取值（zernpy 只接受 [0,2π] 时的桥接依据）。"""
    basis = ZernikeBasis(11)
    rho = np.linspace(0.1, 0.9, 9)
    theta = np.linspace(-np.pi, np.pi, 9, endpoint=False)
    assert np.allclose(
        basis.values(rho, theta),
        basis.values(rho, theta + 2.0 * np.pi),
    )
