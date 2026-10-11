"""报告层：wafer 场图与单文件 HTML（§7.8）。

只读入参与生成物（measurements / fit_result / metrics），**不做任何计算判定** ——
判定在 `metrics`，报告只负责把它画出来、指出来。
"""

from .html import render_html
from .plots import fitted_values, plot_residual_map, plot_wafer_map

__all__ = ["render_html", "fitted_values", "plot_residual_map", "plot_wafer_map"]
