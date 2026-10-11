"""数据契约（§6）：scenario / measurements / truth 的 pydantic 模型与内存表示。

规则（§6）：读写失败必须报错，不静默补默认值。所以所有模型 `extra="forbid"` ——
YAML/JSON 里多一个拼错的键就报错，而不是被悄悄忽略。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from ..errors import (
    InvalidMeasurementError,
    SchemaValidationError,
    SchemaVersionError,
)
from ..geometry import MARK_INSET_MM

SCHEMA_VERSION = 1

# measurements.csv 的列顺序（§6.2，19 列）
MEASUREMENT_COLUMNS: tuple[str, ...] = (
    "wafer_id",
    "lot_id",
    "point_id",
    "field_id",
    "x_mm",
    "y_mm",
    "xf_mm",
    "yf_mm",
    "rho",
    "theta_rad",
    "dx_nm_0",
    "dy_nm_0",
    "dx_nm_180",
    "dy_nm_180",
    "dx_nm",
    "dy_nm",
    "dz_nm",
    "cd_nm",
    "is_edge",
)

# 可为空的列（§6.2 的「空值」列）
NULLABLE_COLUMNS: frozenset[str] = frozenset(
    {"lot_id", "dx_nm_0", "dy_nm_0", "dx_nm_180", "dy_nm_180", "dz_nm", "cd_nm"}
)


class _Spec(BaseModel):
    """所有 scenario 子模型的基类：拒绝未知字段。"""

    model_config = ConfigDict(extra="forbid")


class WaferSpec(_Spec):
    diameter_mm: float = 300.0
    edge_exclusion_mm: float = 3.0

    @model_validator(mode="after")
    def _check(self) -> WaferSpec:
        if not np.isfinite(self.diameter_mm) or self.diameter_mm <= 0:
            raise ValueError(f"wafer.diameter_mm 必须为正的有限值，得到 {self.diameter_mm}")
        radius = self.diameter_mm / 2.0
        if not 0 <= self.edge_exclusion_mm < radius:
            raise ValueError(
                f"wafer.edge_exclusion_mm 必须在 [0, {radius}) 内，得到 {self.edge_exclusion_mm}"
            )
        return self


class FieldSpec(_Spec):
    width_mm: float = 26.0
    height_mm: float = 33.0
    marks_per_field: int = 8

    @model_validator(mode="after")
    def _check(self) -> FieldSpec:
        if min(self.width_mm, self.height_mm) <= 2.0 * MARK_INSET_MM:
            raise ValueError(
                f"field 尺寸 {self.width_mm}×{self.height_mm}mm 放不下距边 {MARK_INSET_MM}mm 的 mark"
            )
        if self.marks_per_field not in (4, 8):
            raise ValueError(
                f"field.marks_per_field={self.marks_per_field} 未实现；"
                "当前只支持 4 与 8（19 点布局的定义自相矛盾，见详设 U-10）"
            )
        return self


class SamplingSpec(_Spec):
    mode: Literal["static", "adaptive"] = "static"
    fields: int = 25
    selection: Literal["regular", "random", "stratified"] = "regular"
    tmu_nm: float = 0.1

    @model_validator(mode="after")
    def _check(self) -> SamplingSpec:
        if self.fields < 1:
            raise ValueError(f"sampling.fields 必须 ≥ 1，得到 {self.fields}")
        if self.tmu_nm < 0 or not np.isfinite(self.tmu_nm):
            raise ValueError(f"sampling.tmu_nm 必须是非负有限值，得到 {self.tmu_nm}")
        return self


class TisSpec(_Spec):
    enabled: bool = True
    amplitude_nm: float = 0.3
    angle_deg: float = 30.0
    calibrated: bool = True

    @model_validator(mode="after")
    def _check(self) -> TisSpec:
        if self.amplitude_nm < 0 or not np.isfinite(self.amplitude_nm):
            raise ValueError(f"tis.amplitude_nm 必须是非负有限值，得到 {self.amplitude_nm}")
        return self


class FingerprintSpec(_Spec):
    basis: Literal["zernike", "cartesian"] = "zernike"
    max_order: int = 12
    terms: dict[str, float] = Field(default_factory=dict)
    field_order: int = 2
    field_terms: dict[str, float] = Field(default_factory=dict)

    @model_validator(mode="after")
    def _check(self) -> FingerprintSpec:
        if self.max_order < 1:
            raise ValueError(f"fingerprint.max_order 必须 ≥ 1，得到 {self.max_order}")
        for term in self.terms:
            if not term.startswith("Z") or not term[1:].isdigit() or int(term[1:]) < 1:
                raise ValueError(
                    f"fingerprint.terms 的键必须是合法 Noll 索引 'Z<n>'（n ≥ 1），得到 {term!r}"
                )
            if int(term[1:]) > self.max_order:
                raise ValueError(
                    f"项 {term!r} 超出 fingerprint.max_order={self.max_order}"
                    "（约束：max_order 必须 ≥ 最大指定项）"
                )
        if self.field_order < 0:
            raise ValueError(f"fingerprint.field_order 必须 ≥ 0，得到 {self.field_order}")
        for term in self.field_terms:
            parts = term.split(",")
            if len(parts) != 2 or not all(p.isdigit() for p in parts):
                raise ValueError(
                    f"fingerprint.field_terms 的键必须是 'p,q'（非负整数），得到 {term!r}"
                )
            if int(parts[0]) + int(parts[1]) > self.field_order:
                raise ValueError(
                    f"项 {term!r} 的总阶数超过 fingerprint.field_order={self.field_order}"
                )
        return self


class FocusSpec(_Spec):
    enabled: bool = False
    terms: dict[str, float] = Field(default_factory=dict)


class BossungSpec(_Spec):
    cd0_nm: float = 32.0
    dose0_mj: float = 20.0
    focus0_nm: float = 0.0
    k_dose: float = 0.8
    k_focus: float = -0.05


class CduSpec(_Spec):
    enabled: bool = False
    bossung: BossungSpec = Field(default_factory=BossungSpec)
    dose_terms: dict[str, float] = Field(default_factory=dict)
    lwr_nm: float = 1.2


class Scenario(_Spec):
    """`scenario.yaml`（§6.1）。

    `wafer_id` / `lot_id` 是相对 §6.1 的补充字段：measurements.csv 需要它们，
    而按 §6.2 又必须有值（lot_id 允许空）。"""""

    version: int = SCHEMA_VERSION
    seed: int = 20261003
    wafer_id: str = "W0001"
    lot_id: str | None = None
    wafer: WaferSpec = Field(default_factory=WaferSpec)
    field: FieldSpec = Field(default_factory=FieldSpec)
    sampling: SamplingSpec = Field(default_factory=SamplingSpec)
    tis: TisSpec = Field(default_factory=TisSpec)
    fingerprint: FingerprintSpec = Field(default_factory=FingerprintSpec)
    focus: FocusSpec = Field(default_factory=FocusSpec)
    cdu: CduSpec = Field(default_factory=CduSpec)

    @field_validator("version")
    @classmethod
    def _check_version(cls, value: int) -> int:
        if value != SCHEMA_VERSION:
            raise SchemaVersionError(
                f"scenario version={value} 与当前代码的 {SCHEMA_VERSION} 不兼容"
            )
        return value

    @property
    def wafer_radius_mm(self) -> float:
        return self.wafer.diameter_mm / 2.0


class Truth(BaseModel):
    """`truth.json`（§6.3）。**只有 metrics 模块可以读它。**"""

    model_config = ConfigDict(extra="forbid")

    version: int = SCHEMA_VERSION
    seed: int
    wafer_id: str
    basis: str
    max_order: int
    wafer_terms: dict[str, float]
    wafer_terms_y: dict[str, float]
    field_terms_x: dict[str, float]
    field_terms_y: dict[str, float]
    tis_vector_nm: tuple[float, float]
    noise_sigma_nm: float
    generated_at: str
    eumetria_version: str


@dataclass(frozen=True)
class MeasurementSet:
    """measurements.csv 的内存表示（19 列，numpy 数组）。

    `wafer_radius_mm` / `field_width_mm` / `field_height_mm` 不是文件里的列，但读入后
    必须携带：§6.2.5 要求重算 rho/theta 校验，而 rho 的分母是硅片半径；
    fit 算 ξ/η 也要用 field 尺寸。不引 pandas —— 188 行的表不需要它（P-02）。
    """

    wafer_id: str
    point_id: NDArray[np.int64]
    field_id: NDArray[np.int64]
    x_mm: NDArray[np.float64]
    y_mm: NDArray[np.float64]
    xf_mm: NDArray[np.float64]
    yf_mm: NDArray[np.float64]
    rho: NDArray[np.float64]
    theta_rad: NDArray[np.float64]
    dx_nm: NDArray[np.float64]
    dy_nm: NDArray[np.float64]
    is_edge: NDArray[np.bool_]
    wafer_radius_mm: float
    field_width_mm: float
    field_height_mm: float
    lot_id: str | None = None
    dx_nm_0: NDArray[np.float64] | None = None
    dy_nm_0: NDArray[np.float64] | None = None
    dx_nm_180: NDArray[np.float64] | None = None
    dy_nm_180: NDArray[np.float64] | None = None
    dz_nm: NDArray[np.float64] | None = None
    cd_nm: NDArray[np.float64] | None = None

    def __len__(self) -> int:
        return int(self.x_mm.size)

    def xi_eta(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """场内归一化坐标 (ξ, η) ∈ [-1,1]^2。"""
        return (
            2.0 * self.xf_mm / self.field_width_mm,
            2.0 * self.yf_mm / self.field_height_mm,
        )

    def check_consistency(self) -> None:
        """内部一致性：长度对齐、索引唯一（§6.2.5 第 4 条的通用版本）。"""
        n = len(self)
        for name in (
            "point_id",
            "field_id",
            "x_mm",
            "y_mm",
            "xf_mm",
            "yf_mm",
            "rho",
            "theta_rad",
            "dx_nm",
            "dy_nm",
            "is_edge",
        ):
            arr = getattr(self, name)
            if arr.shape != (n,):
                raise SchemaValidationError(f"{name} 长度为 {arr.shape}，期望 ({n},)")
        for name in ("dx_nm_0", "dy_nm_0", "dx_nm_180", "dy_nm_180", "dz_nm", "cd_nm"):
            arr = getattr(self, name)
            if arr is not None and arr.shape != (n,):
                raise SchemaValidationError(f"{name} 长度为 {arr.shape}，期望 ({n},)")
        if len(np.unique(self.point_id)) != n:
            raise InvalidMeasurementError("point_id 必须全局唯一")
        seen: set[tuple[int, float, float]] = set()
        # strict=True：三列长度不一致时宁可报错，也不静默截断把重复点放过去
        for fid, xf, yf in zip(self.field_id, self.xf_mm, self.yf_mm, strict=True):
            key = (int(fid), round(float(xf), 9), round(float(yf), 9))
            if key in seen:
                raise InvalidMeasurementError(
                    f"同一 field {int(fid)} 内 (xf_mm, yf_mm)=({xf:.6f}, {yf:.6f}) 重复"
                )
            seen.add(key)
