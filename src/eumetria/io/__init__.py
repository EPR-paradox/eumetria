"""io 层：数据契约的校验、读写。

依赖：pydantic + pyyaml + numpy（§4.2）。**不依赖 synth / fit** —— 反向依赖会让
「只有生成器知道真值」这条隔离失效（io 被 fit 用到，若 io 反过来 import synth，
fit 就能顺着 io 摸到真值）。
"""

from .readers import load_measurements, load_scenario, load_truth
from .schemas import (
    MEASUREMENT_COLUMNS,
    SCHEMA_VERSION,
    MeasurementSet,
    Scenario,
    Truth,
)
from .writers import now_cst_isoformat, write_measurements, write_truth

__all__ = [
    "MEASUREMENT_COLUMNS",
    "SCHEMA_VERSION",
    "MeasurementSet",
    "Scenario",
    "Truth",
    "load_measurements",
    "load_scenario",
    "load_truth",
    "now_cst_isoformat",
    "write_measurements",
    "write_truth",
]
