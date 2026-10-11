"""拟合层：两级模型的联合最小二乘（D-09）与可辨识性诊断。

依赖：basis + geometry + io.schemas + errors。**不依赖 synth** —— 拟合函数拿不到真值，
真值只在 metrics 的显式入参里出现（§1.2 的隔离规则）。
"""

from .linear_fit import FitConfig, build_design, field_terms_for, fit_field
from .result import FieldDiagnostics, FitResult, IdentifiabilityInfo
from .writers import write_fit_result

__all__ = [
    "FitConfig",
    "FitResult",
    "FieldDiagnostics",
    "IdentifiabilityInfo",
    "build_design",
    "field_terms_for",
    "fit_field",
    "write_fit_result",
]
