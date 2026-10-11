"""与 zernpy 的数值一致性对照（§13.2，可选依赖，未安装则 skip）。

zernpy 只是**对照后端**（P-03）：本项目的 Zernike 由自己实现，这个测试用来证明
自己实现没抄错，而不是把计算结果建立在第三方之上。
"""

from __future__ import annotations

import numpy as np
import pytest

from eumetria.basis import ZernikeBasis, noll_to_nm

zernpy = pytest.importorskip("zernpy", reason="可选依赖 zernpy 未安装")

TOL = 1e-8  # §13.2 规定的容差
MAX_NOLL = 36


@pytest.mark.optional
def test_noll_index_mapping_matches_zernpy():
    for j in range(1, MAX_NOLL + 1):
        polym = zernpy.ZernPol(noll_index=j)
        (m, n), _osa, noll, _fringe = polym.get_indices()
        assert noll == j
        assert noll_to_nm(j) == (n, m), f"j={j} 的 (n, m) 与 zernpy 不一致"


@pytest.mark.optional
def test_values_match_zernpy_on_random_points():
    rng = np.random.default_rng(20261003)
    rho = np.sqrt(rng.uniform(0.0, 1.0, 64))
    theta = rng.uniform(0.0, 2.0 * np.pi, 64)  # zernpy 要求 [0, 2π]

    basis = ZernikeBasis(max_order=MAX_NOLL)
    mine = basis.values(rho, theta)
    theirs = np.column_stack(
        [
            np.asarray(zernpy.ZernPol(noll_index=j).polynomial_value(rho, theta), dtype=float)
            for j in range(1, MAX_NOLL + 1)
        ]
    )
    deviation = np.abs(mine - theirs)
    assert deviation.max() < TOL, f"与 zernpy 的最大偏差 {deviation.max():.3e} 超过 {TOL}"


@pytest.mark.optional
def test_values_match_zernpy_at_boundary_points():
    rho = np.array([0.0, 1.0, 1.0, 0.5])
    theta = np.array([0.0, 0.0, np.pi, 3.0 * np.pi / 2.0])
    basis = ZernikeBasis(max_order=15)
    for col, j in enumerate(range(1, 16)):
        expected = np.asarray(
            zernpy.ZernPol(noll_index=j).polynomial_value(rho, theta), dtype=float
        )
        assert np.allclose(basis.values(rho, theta)[:, col], expected, atol=TOL)
