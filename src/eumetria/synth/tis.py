"""TIS（Tool Induced Shift）模型与 0°/180° 双次测量的组合（§3.4）。

量测工具自身带系统偏移 TIS。真实产线对同一个 mark 做两次测量：一次正常，一次把
硅片（或 mark）转 180°。转 180° 让**真实值反号**，而 TIS 跟着工具走、**不反号**：

    m_0   = +true + TIS + ε_0
    m_180 = −true + TIS + ε_180

于是

    true_hat = (m_0 − m_180) / 2        组合掉 TIS
    TIS_hat  = (m_0 + m_180) / 2        被丢弃的那一半，可用于校准验证

噪声在这里加（两次测量各有独立噪声），不要在外面再加一次 —— 否则 σ 会被放大 √2 倍。
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

Arr: TypeAlias = NDArray[np.float64]


def tis_vector_from_polar(amplitude_nm: float, angle_deg: float) -> tuple[float, float]:
    """由幅值（nm）与方向（deg）得到 TIS 向量 (dx_nm, dy_nm)。方向固定、不随机制造。"""
    if amplitude_nm < 0:
        raise ValueError(f"TIS 幅值必须非负，得到 {amplitude_nm}")
    angle_rad = np.deg2rad(angle_deg)
    return float(amplitude_nm * np.cos(angle_rad)), float(amplitude_nm * np.sin(angle_rad))


def apply_tis(
    true: Arr,
    tis_vector_nm: tuple[float, float],
    calibrated: bool,
    sigma_nm: float,
    rng: np.random.Generator,
) -> dict[str, Arr]:
    """对真值施加 TIS 与噪声。

    calibrated=True  → 返回 {"m0": …, "m180": …}，两组测量，各自独立噪声
    calibrated=False → 返回 {"m": …}，只有一组带 TIS 的读数（Q-07 的对照臂）

    返回的数组形状与 `true` 相同（(n_points, 2)）。
    """
    true = np.asarray(true, dtype=np.float64)
    if true.ndim != 2 or true.shape[1] != 2:
        raise ValueError(f"true 必须是 (n_points, 2) 的 (dx, dy)，得到 {true.shape}")
    if sigma_nm < 0 or not np.isfinite(sigma_nm):
        raise ValueError(f"sigma_nm 必须是非负有限值，得到 {sigma_nm}")

    tis = np.asarray(tis_vector_nm, dtype=np.float64).reshape(1, 2)

    if not calibrated:
        eps = rng.normal(0.0, sigma_nm, size=true.shape) if sigma_nm else np.zeros_like(true)
        return {"m": true + tis + eps}

    eps_0 = rng.normal(0.0, sigma_nm, size=true.shape) if sigma_nm else np.zeros_like(true)
    eps_180 = rng.normal(0.0, sigma_nm, size=true.shape) if sigma_nm else np.zeros_like(true)
    return {
        "m0": true + tis + eps_0,
        "m180": -true + tis + eps_180,
    }


def combine_tis(m0: Arr, m180: Arr) -> tuple[Arr, Arr]:
    """返回 (measured, tis_hat)：前者进拟合，后者用于校准验证。"""
    m0 = np.asarray(m0, dtype=np.float64)
    m180 = np.asarray(m180, dtype=np.float64)
    if m0.shape != m180.shape:
        raise ValueError(f"m0 与 m180 形状不一致：{m0.shape} vs {m180.shape}")
    return (m0 - m180) / 2.0, (m0 + m180) / 2.0
