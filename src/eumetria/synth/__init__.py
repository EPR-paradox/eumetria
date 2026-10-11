"""正向生成层：布局 → 采样 → 指纹 → TIS/噪声 → 测量值。

这是整条流水线的唯一真值来源。只有本层的生成函数能看到真值（§1.2 的硬规则）。
"""

from .fingerprint import SamplePoints, make_fingerprint, sample_fields, sample_points
from .noise import add_noise
from .pipeline import GeneratedRun, generate
from .tis import apply_tis, combine_tis, tis_vector_from_polar

__all__ = [
    "GeneratedRun",
    "SamplePoints",
    "add_noise",
    "apply_tis",
    "combine_tis",
    "generate",
    "make_fingerprint",
    "sample_fields",
    "sample_points",
    "tis_vector_from_polar",
]
