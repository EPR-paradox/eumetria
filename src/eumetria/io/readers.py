"""读取与强制校验（§6.2.5）。

`measurements.csv` 载入后必须做四件事，一件都不能省：

1. 重算 rho/theta 与文件里的值比对（容差 1e-9）——不让文件伪造坐标系
2. 断言 rho ≤ 1.0，越界抛 `InvalidMeasurementError`
3. 若四列 `_0`/`_180` 都在，断言 `dx_nm == (dx_nm_0 − dx_nm_180)/2`
4. 断言同一 `field_id` 内 `(xf_mm, yf_mm)` 不重复

为什么第 1 条重要：拟合只看得见 `measurements.csv`，如果 (ρ, θ) 可以由文件随便写，
那么「真值的唯一来源是生成器」这条隔离（§1.2）就被绕过了 —— 谁都能塞一份手改的
极坐标进来，让拟合"拟合出"想要的结果。
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import yaml

from ..errors import (
    InvalidMeasurementError,
    SchemaValidationError,
    SchemaVersionError,
)
from ..geometry import to_polar
from .schemas import (
    MEASUREMENT_COLUMNS,
    NULLABLE_COLUMNS,
    SCHEMA_VERSION,
    MeasurementSet,
    Scenario,
    Truth,
)

RHO_TOL = 1e-9


def load_scenario(path: str | Path) -> Scenario:
    """读 `scenario.yaml`。未知字段、非法项名、越界取值一律报错（extra='forbid'）。"""
    raw = Path(path).read_text(encoding="utf-8")
    try:
        data = yaml.safe_load(raw)
    except yaml.YAMLError as exc:  # pragma: no cover - 取决于文件内容
        raise SchemaValidationError(f"scenario YAML 解析失败：{exc}") from exc
    if not isinstance(data, dict):
        raise SchemaValidationError(f"scenario 必须是映射，得到 {type(data).__name__}")
    _check_version(data.get("version"), "scenario")
    return Scenario.model_validate(data)


def load_truth(path: str | Path) -> Truth:
    """读 `truth.json`。**调用方必须是 metrics 模块**（§6.3 的隔离规则）。"""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise SchemaValidationError(f"{path} 必须是 JSON 对象")
    _check_version(data.get("version"), "truth")
    return Truth.model_validate(data)


def _check_version(value: object, kind: str) -> None:
    """`version` 必须在进 pydantic 之前检查 —— 否则 SchemaVersionError 会被
    pydantic 包成 ValidationError，调用方拿到的异常类型就跟 §6.6 的约定不符了。"""
    if value != SCHEMA_VERSION:
        raise SchemaVersionError(
            f"{kind} 的 version={value!r} 与当前代码的 {SCHEMA_VERSION} 不兼容"
        )


def _require_columns(fieldnames: Sequence[str] | None, path: Path) -> None:
    if not fieldnames:
        raise SchemaValidationError(f"{path} 没有表头")
    missing = [c for c in MEASUREMENT_COLUMNS if c not in fieldnames]
    if missing:
        raise SchemaValidationError(f"{path} 缺少列：{missing}")
    unexpected = [c for c in fieldnames if c not in MEASUREMENT_COLUMNS]
    if unexpected:  # §6.6：增列允许（向后兼容），但不允许悄悄出现拼错的列
        raise SchemaValidationError(
            f"{path} 出现未定义的列：{unexpected}（新增列必须同时更新 MEASUREMENT_COLUMNS）"
        )


def _floats(rows: list[dict[str, str]], column: str, path: Path) -> np.ndarray | None:
    if column in NULLABLE_COLUMNS:
        values = [r[column] for r in rows]
        if all(v == "" for v in values):
            return None
        try:
            return np.array([np.nan if v == "" else float(v) for v in values], dtype=np.float64)
        except ValueError as exc:
            raise SchemaValidationError(f"{path} 的列 {column} 含非数值：{exc}") from exc
    try:
        return np.array([float(r[column]) for r in rows], dtype=np.float64)
    except ValueError as exc:
        raise SchemaValidationError(f"{path} 的列 {column} 含非数值：{exc}") from exc


def _ints(rows: list[dict[str, str]], column: str, path: Path) -> np.ndarray:
    try:
        return np.array([int(r[column]) for r in rows], dtype=np.int64)
    except ValueError as exc:
        raise SchemaValidationError(f"{path} 的列 {column} 含非整数：{exc}") from exc


def _bools(rows: list[dict[str, str]], column: str, path: Path) -> np.ndarray:
    allowed = {"true": True, "false": False, "1": True, "0": False}
    out: list[bool] = []
    for r in rows:
        token = r[column].strip().lower()
        if token not in allowed:
            raise SchemaValidationError(
                f"{path} 的列 {column} 只能是 true/false，得到 {r[column]!r}"
            )
        out.append(allowed[token])
    return np.array(out, dtype=np.bool_)


def load_measurements(
    path: str | Path,
    *,
    wafer_radius_mm: float,
    field_width_mm: float,
    field_height_mm: float,
) -> MeasurementSet:
    """读 `measurements.csv` 并做 §6.2.5 的四条强制校验。

    `wafer_radius_mm` 是必需的：rho = r / R_wafer，没有 R 就无法重算 rho，
    也就无法完成第 1 条校验。调用方通常从 scenario 里取。
    """
    path = Path(path)
    if not np.isfinite(wafer_radius_mm) or wafer_radius_mm <= 0:
        raise SchemaValidationError(f"wafer_radius_mm 必须为正的有限值，得到 {wafer_radius_mm}")

    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        _require_columns(reader.fieldnames, path)
        rows = list(reader)
    if not rows:
        raise InvalidMeasurementError(f"{path} 没有数据行")

    wafer_ids = {r["wafer_id"] for r in rows}
    if len(wafer_ids) != 1:
        raise InvalidMeasurementError(f"{path} 含多个 wafer_id：{sorted(wafer_ids)}")
    lot_ids = {r["lot_id"] for r in rows}
    if len(lot_ids) != 1:
        raise InvalidMeasurementError(f"{path} 含多个 lot_id：{sorted(lot_ids)}")

    ms = MeasurementSet(
        wafer_id=rows[0]["wafer_id"],
        lot_id=rows[0]["lot_id"] or None,
        point_id=_ints(rows, "point_id", path),
        field_id=_ints(rows, "field_id", path),
        x_mm=_floats(rows, "x_mm", path),
        y_mm=_floats(rows, "y_mm", path),
        xf_mm=_floats(rows, "xf_mm", path),
        yf_mm=_floats(rows, "yf_mm", path),
        rho=_floats(rows, "rho", path),
        theta_rad=_floats(rows, "theta_rad", path),
        dx_nm=_floats(rows, "dx_nm", path),
        dy_nm=_floats(rows, "dy_nm", path),
        is_edge=_bools(rows, "is_edge", path),
        dx_nm_0=_floats(rows, "dx_nm_0", path),
        dy_nm_0=_floats(rows, "dy_nm_0", path),
        dx_nm_180=_floats(rows, "dx_nm_180", path),
        dy_nm_180=_floats(rows, "dy_nm_180", path),
        dz_nm=_floats(rows, "dz_nm", path),
        cd_nm=_floats(rows, "cd_nm", path),
        wafer_radius_mm=float(wafer_radius_mm),
        field_width_mm=float(field_width_mm),
        field_height_mm=float(field_height_mm),
    )
    _validate(ms, path)
    return ms


def _validate(ms: MeasurementSet, path: Path) -> None:
    # --- 第 1 条：重算 rho/theta 并比对 ---
    rho_calc, theta_calc = to_polar(ms.x_mm, ms.y_mm, ms.wafer_radius_mm)
    if not np.allclose(ms.rho, rho_calc, rtol=0.0, atol=RHO_TOL):
        worst = float(np.max(np.abs(ms.rho - rho_calc)))
        raise InvalidMeasurementError(
            f"{path} 的 rho 与重算值不一致（最大偏差 {worst:.3e} > {RHO_TOL:.0e}）"
        )
    if not np.allclose(ms.theta_rad, theta_calc, rtol=0.0, atol=RHO_TOL):
        worst = float(np.max(np.abs(ms.theta_rad - theta_calc)))
        raise InvalidMeasurementError(
            f"{path} 的 theta_rad 与重算值不一致（最大偏差 {worst:.3e} > {RHO_TOL:.0e}）"
        )

    # --- 第 2 条：rho 不得越界 ---
    if np.any(ms.rho > 1.0 + RHO_TOL):
        worst = float(np.max(ms.rho))
        raise InvalidMeasurementError(
            f"{path} 有 rho > 1 的点（最大 {worst:.6f}）；越界点应在采样阶段排除"
        )

    # --- 第 3 条：四列都在则校验 TIS 组合 ---
    have_all = all(
        arr is not None
        for arr in (ms.dx_nm_0, ms.dy_nm_0, ms.dx_nm_180, ms.dy_nm_180)
    )
    if have_all:
        for a0, a180, combined, name in (
            (ms.dx_nm_0, ms.dx_nm_180, ms.dx_nm, "dx"),
            (ms.dy_nm_0, ms.dy_nm_180, ms.dy_nm, "dy"),
        ):
            if a0 is None or a180 is None:  # pragma: no cover - have_all 已保证
                raise InvalidMeasurementError(f"{path} 的 {name} 双测列不完整")
            expected = (a0 - a180) / 2.0
            if not np.allclose(combined, expected, rtol=0.0, atol=RHO_TOL):
                worst = float(np.max(np.abs(combined - expected)))
                raise InvalidMeasurementError(
                    f"{path} 的 {name}_nm 与 ({name}_nm_0 − {name}_nm_180)/2 不一致"
                    f"（最大偏差 {worst:.3e}）"
                )

    # --- 第 4 条 + 长度/索引一致性 ---
    ms.check_consistency()
