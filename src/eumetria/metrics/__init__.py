"""评分层：系数恢复误差与残差一致性（§7.7）。

**唯一允许读 `truth.json` 的模块**（§6.3）。真值永远是显式入参，
所以 `fit` 无法通过任何路径碰到它。
"""

from .recovery import MetricsResult, coefficient_recovery, compare, residual_rms

__all__ = [
    "MetricsResult",
    "coefficient_recovery",
    "compare",
    "residual_rms",
]
