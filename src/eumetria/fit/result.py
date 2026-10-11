"""单个场 / 整体可辨识性的诊断结构（对应 §6.4 的 `identifiability` 与逐场证据）。

`FieldDiagnostics.max_identifiable_order` 是 §3.1.1「CPE 阶数上限由最差场决定」的直接证据：
它逐场回答"这个场最多能支撑几阶场内基"，整体上限取所有场的最小值。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from ..io.schemas import SCHEMA_VERSION


@dataclass(frozen=True)
class FieldDiagnostics:
    """一个场的可辨识性诊断。"""

    field_id: int
    n_marks: int
    configured_order: int
    rank: int
    n_params: int
    condition_number: float
    max_identifiable_order: int

    @property
    def is_identifiable(self) -> bool:
        """配置的那一阶能不能辨识（rank 满）。"""
        return self.rank == self.n_params

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_marks": self.n_marks,
            "configured_order": self.configured_order,
            "rank": self.rank,
            "n_params": self.n_params,
            "condition_number": (
                None if self.condition_number == float("inf") else round(self.condition_number, 6)
            ),
            "identifiable": self.is_identifiable,
            "max_identifiable_order": self.max_identifiable_order,
        }


@dataclass(frozen=True)
class IdentifiabilityInfo:
    """可辨识性汇总。`warnings` 非空不等于失败，`rank_deficient` 才是失败。"""

    max_identifiable_order: int
    cpe_order_ceiling: int
    warnings: list[str] = field(default_factory=list)
    rank_deficient: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "max_identifiable_order": self.max_identifiable_order,
            "cpe_order_ceiling": self.cpe_order_ceiling,
            "warnings": list(self.warnings),
            "rank_deficient": self.rank_deficient,
        }


@dataclass(frozen=True)
class FitResult:
    """一次拟合的全部输出。字段与 §6.4 对齐，可直接写 JSON。

    `linear_terms` / `correctable` / `cpe_terms` / `uncorrectable_residual_rms_nm` 由 Phase 3 的
    `decompose` 填写，这里不出现（不用 0 冒充已知）。
    """

    basis: str
    max_order: int
    selected_terms: list[str]
    selection_method: str
    wafer_terms: list[str]
    field_terms: list[str]
    wafer_coeffs_x: dict[str, float]
    wafer_coeffs_y: dict[str, float]
    field_coeffs_x: dict[str, float]
    field_coeffs_y: dict[str, float]
    residual_rms_nm: float
    residual_rms_nm_x: float
    residual_rms_nm_y: float
    n_points: int
    n_params: int
    condition_number: float
    identifiability: IdentifiabilityInfo
    field_diagnostics: dict[int, FieldDiagnostics]
    version: int = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "basis": self.basis,
            "max_order": self.max_order,
            "selected_terms": list(self.selected_terms),
            "selection_method": self.selection_method,
            "wafer_coeffs_x": _round_dict(self.wafer_coeffs_x),
            "wafer_coeffs_y": _round_dict(self.wafer_coeffs_y),
            "field_coeffs_x": _round_dict(self.field_coeffs_x),
            "field_coeffs_y": _round_dict(self.field_coeffs_y),
            "residual_rms_nm": round(self.residual_rms_nm, 6),
            "residual_rms_nm_x": round(self.residual_rms_nm_x, 6),
            "residual_rms_nm_y": round(self.residual_rms_nm_y, 6),
            "n_points": self.n_points,
            "n_params": self.n_params,
            "condition_number": round(self.condition_number, 6),
            "identifiability": self.identifiability.to_dict(),
            "field_diagnostics": {
                str(k): v.to_dict() for k, v in sorted(self.field_diagnostics.items())
            },
        }


def _round_dict(values: dict[str, float]) -> dict[str, float]:
    return {k: round(float(v), 6) for k, v in values.items()}
