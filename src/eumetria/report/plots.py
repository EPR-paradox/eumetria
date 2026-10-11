"""基础图（§7.8）：wafer map / 残差图 / 矢量图。

约定：
- matplotlib 走 `Agg` 后端（本机无显示环境），且**惰性导入** —— 它只是 pyproject 里的
  `[plot]` 额外依赖，缺了要给一句能读的错，而不是 ImportError 堆栈。
- 图上必须画出硅片边界与 edge exclusion 圆：量的空间含义靠这两条线读。
- 单位一律写在坐标轴标签里（mm / nm），不靠标题暗示。
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, TypeAlias

import numpy as np
from numpy.typing import NDArray

from ..basis import CartesianBasis, ZernikeBasis
from ..io.schemas import MeasurementSet

Arr: TypeAlias = NDArray[np.float64]


_CJK_PREFERENCE = (
    "Noto Sans CJK SC",
    "Noto Sans CJK JP",  # 覆盖率与 SC 相同，个别字形是日式变体，可读
    "Source Han Sans SC",
    "Droid Sans Fallback",
    "AR PL UMing CN",
    "WenQuanYi Zen Hei",
)


def _plt():  # pragma: no cover - 取决于环境
    try:
        import matplotlib
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "绘图需要 matplotlib：uv pip install --python .venv/bin/python matplotlib"
        ) from exc
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager

    # 图里的标题/标签是中文，DejaVu Sans 没有 CJK 字形，会画成一排豆腐块。
    # 按可用性挑一个中文字体排在最前；一个都没有时不报错，只是字会缺 —— 值得在
    # 日志里说一句，否则"图上的字是方块"会被当成绘图 bug 查很久。
    available = {f.name for f in font_manager.fontManager.ttflist}
    chosen = [name for name in _CJK_PREFERENCE if name in available]
    plt.rcParams["font.sans-serif"] = [*chosen, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def fitted_values(fit: Any, measurements: MeasurementSet) -> tuple[Arr, Arr]:
    """由拟合系数重建每个点上的模型值（用于残差图）。"""
    def _get(key: str):
        if isinstance(fit, Mapping):
            return fit[key]
        return getattr(fit, key)

    w_map = {str(k): float(v) for k, v in dict(_get("wafer_coeffs_x")).items()}
    f_map = {str(k): float(v) for k, v in dict(_get("field_coeffs_x")).items()}
    w_terms = list(w_map)
    f_terms = list(f_map)
    zw = ZernikeBasis(int(_get("max_order"))).values(
        measurements.rho, measurements.theta_rad, w_terms
    )
    xi, eta = measurements.xi_eta()
    order = max((int(t.split(",")[0]) + int(t.split(",")[1]) for t in f_terms), default=0)
    cf = CartesianBasis(order).values(xi, eta, f_terms)
    dx = zw @ np.array([w_map[t] for t in w_terms]) + cf @ np.array([f_map[t] for t in f_terms])

    w_map_y = {str(k): float(v) for k, v in dict(_get("wafer_coeffs_y")).items()}
    f_map_y = {str(k): float(v) for k, v in dict(_get("field_coeffs_y")).items()}
    dx2 = dx
    dy = (
        ZernikeBasis(int(_get("max_order"))).values(
            measurements.rho, measurements.theta_rad, list(w_map_y)
        )
        @ np.array([w_map_y[t] for t in w_map_y])
        + CartesianBasis(order).values(xi, eta, list(f_map_y))
        @ np.array([f_map_y[t] for t in f_map_y])
    )
    return np.asarray(dx2, dtype=np.float64), np.asarray(dy, dtype=np.float64)


def _draw_wafer_outline(ax, radius_mm: float, edge_exclusion_mm: float) -> None:
    theta = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.plot(radius_mm * np.cos(theta), radius_mm * np.sin(theta), color="0.25", lw=1.0)
    usable = radius_mm - edge_exclusion_mm
    ax.plot(
        usable * np.cos(theta),
        usable * np.sin(theta),
        color="0.6",
        lw=0.8,
        ls="--",
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")


def plot_wafer_map(
    values: Arr,
    x_mm: Arr,
    y_mm: Arr,
    *,
    out_path: str | Path,
    title: str = "",
    label: str = "value (nm)",
    radius_mm: float = 150.0,
    edge_exclusion_mm: float = 3.0,
    kind: str = "scalar",
    dx: Arr | None = None,
    dy: Arr | None = None,
) -> Path:
    """画一张 wafer 场图。`kind='scalar'` 用 `values` 上色，`kind='quiver'` 用 (dx, dy) 画箭头。"""
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7.2, 7.2), dpi=140)

    if kind == "scalar":
        scatter = ax.scatter(x_mm, y_mm, c=np.asarray(values), s=26, cmap="coolwarm", zorder=3)
        fig.colorbar(scatter, ax=ax, label=label, shrink=0.8)
    elif kind == "quiver":
        if dx is None or dy is None:
            raise ValueError("kind='quiver' 需要同时给 dx 与 dy")
        ax.quiver(x_mm, y_mm, dx, dy, np.hypot(dx, dy), cmap="coolwarm", scale=200, zorder=3)
    else:
        raise ValueError(f"未知的 kind={kind!r}，只支持 'scalar' / 'quiver'")

    _draw_wafer_outline(ax, radius_mm, edge_exclusion_mm)
    if title:
        ax.set_title(title)
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def plot_residual_map(
    resid_x: Arr,
    resid_y: Arr,
    x_mm: Arr,
    y_mm: Arr,
    *,
    out_path: str | Path,
    radius_mm: float = 150.0,
    edge_exclusion_mm: float = 3.0,
) -> Path:
    """残差图：左 X 残差、右 Y 残差（同一色标，便于直接比）。"""
    plt = _plt()
    limit = float(max(np.abs(resid_x).max(), np.abs(resid_y).max()))
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 6.6), dpi=140)
    for ax, resid, name in zip(axes, (resid_x, resid_y), ("X", "Y")):
        scatter = ax.scatter(
            x_mm, y_mm, c=np.asarray(resid), s=26, cmap="coolwarm",
            vmin=-limit, vmax=limit, zorder=3,
        )
        _draw_wafer_outline(ax, radius_mm, edge_exclusion_mm)
        ax.set_title(f"{name} 残差")
        fig.colorbar(scatter, ax=ax, label="residual (nm)", shrink=0.8)
    fig.suptitle("拟合残差场（已扣除 wafer + 场内模型）")
    fig.tight_layout()
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
