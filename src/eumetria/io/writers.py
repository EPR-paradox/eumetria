"""写出（§6.2 / §6.3）。

**精度：全量写出（12 位有效数字），不四舍五入。** §6.2.3 里那些 4 位小数只是给人看的
展示格式；文件里要是只写 4 位，`readers` 的 §6.2.5 第 1 条校验（rho/theta 重算比对，
容差 1e-9）就会失败。写和读用同一套精度约定。
"""

from __future__ import annotations

import csv
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

from ..errors import SchemaValidationError
from .schemas import MEASUREMENT_COLUMNS, MeasurementSet, Truth

_FLOAT_FMT = ""  # 浮点走 repr()：最短且精确往返（见 _fmt）
_CST = timezone(timedelta(hours=8))


def _fmt(value: float | int | bool | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not np.isfinite(number):
            raise SchemaValidationError(f"待写出的数值不是有限值：{value}")
        # repr(float) 是 Python 保证精确往返的最短表示（17 位有效数字以内按需给出）。
        # 用 .12g 之类的固定精度会在 |v| ~ 2.5 这种量级上把误差推到 5e-11，
        # 直接顶穿 readers 里 1e-9 的 rho/theta 重算校验。
        return repr(number)
    return str(value)


def write_measurements(path: str | Path, ms: MeasurementSet) -> Path:
    """写 `measurements.csv`（19 列，列序固定，生成物不入 git）。"""
    ms.check_consistency()
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    n = len(ms)

    columns: dict[str, list[str]] = {
        "wafer_id": [ms.wafer_id] * n,
        "lot_id": [ms.lot_id or ""] * n,
        "point_id": [_fmt(v) for v in ms.point_id],
        "field_id": [_fmt(v) for v in ms.field_id],
        "x_mm": [_fmt(v) for v in ms.x_mm],
        "y_mm": [_fmt(v) for v in ms.y_mm],
        "xf_mm": [_fmt(v) for v in ms.xf_mm],
        "yf_mm": [_fmt(v) for v in ms.yf_mm],
        "rho": [_fmt(v) for v in ms.rho],
        "theta_rad": [_fmt(v) for v in ms.theta_rad],
        "dx_nm": [_fmt(v) for v in ms.dx_nm],
        "dy_nm": [_fmt(v) for v in ms.dy_nm],
        "is_edge": [_fmt(bool(v)) for v in ms.is_edge],
    }
    for name in ("dx_nm_0", "dy_nm_0", "dx_nm_180", "dy_nm_180", "dz_nm", "cd_nm"):
        arr = getattr(ms, name)
        columns[name] = [""] * n if arr is None else [_fmt(v) for v in arr]

    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(MEASUREMENT_COLUMNS)
        for i in range(n):
            writer.writerow([columns[c][i] for c in MEASUREMENT_COLUMNS])
    return path


def write_truth(path: str | Path, truth: Truth) -> Path:
    """写 `truth.json`（§6.3）。生成物不入 git。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = truth.model_dump(mode="json")
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return path


def now_cst_isoformat() -> str:
    """生成本地时间戳（CST），供 truth.json 的 `generated_at` 使用。"""
    return datetime.now(_CST).replace(microsecond=0).isoformat()
