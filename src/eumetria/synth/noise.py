"""量测噪声模型（§3.4）。

    measured = true + TIS + N(0, σ_tmu)

所有随机数走调用方传入的 `np.random.Generator`（§9.5：生成器函数不接受全局随机状态），
所以调用顺序决定数值结果 —— 流水线里的顺序在 §6.2.4 固定，不要改动。
"""

from __future__ import annotations

from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

Arr: TypeAlias = NDArray[np.float64]


def add_noise(true: Arr, sigma_nm: float, rng: np.random.Generator) -> Arr:
    """加零均值高斯噪声，形状与 `true` 相同。σ = 0 时原样返回（便于无噪声验收 A-02）。"""
    if sigma_nm < 0:
        raise ValueError(f"sigma_nm 必须非负，得到 {sigma_nm}")
    if not np.isfinite(sigma_nm):
        raise ValueError(f"sigma_nm 必须是有限值，得到 {sigma_nm}")
    true = np.asarray(true, dtype=np.float64)
    if sigma_nm == 0.0:
        return true.copy()
    return true + rng.normal(0.0, sigma_nm, size=true.shape)
