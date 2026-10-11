"""评分：把拟合结果与真值对比（§7.7、§6.5）。

**只有本模块可以读 `truth.json`**（§6.3）。真值只作为显式入参出现 —— `fit` 的函数签名里
不存在 truth，这是"评分器不知道真值"的代码级保证。

与详设的两点差异（有意）：

1. `sampling_budget_curve`（Q-03）不在本期实现，属 Phase 2。
2. A-02 / A-03 的判定需要知道「哪些项是真值里非零的」，因此从 `truth` 的系数字典里取，
   而不是从 `fit` 的项集里取 —— 否则"漏掉一整项"会被误判成满分。
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from ..basis import CartesianBasis, ZernikeBasis
from ..io.schemas import SCHEMA_VERSION, MeasurementSet, Truth

# A-02：无噪声 + 点数 ≥ 8× 系数个数时，全部系数相对误差 < 1%
A02_REL_TOL = 0.01
A02_POINT_RATIO = 8
# A-03：TMU 下前 6 阶主要项相对误差 < 10%
A03_REL_TOL = 0.10
A03_MAX_ORDER = 6
# 真值里没给、但拟合却给出显著值的项（假阳性）的提示阈值（nm）
SPURIOUS_TOL_NM = 0.05


def coefficient_recovery(
    fit: Mapping[str, float], truth: Mapping[str, float]
) -> dict[str, dict[str, float]]:
    """逐项对照。分母用 |truth|，truth 为 0 的项不做相对误差（会除零）。"""
    out: dict[str, dict[str, float]] = {}
    for term, truth_value in truth.items():
        fit_value = float(fit.get(term, 0.0))
        abs_err = abs(fit_value - truth_value)
        rel_err = abs_err / abs(truth_value) if truth_value != 0 else float("nan")
        out[term] = {
            "truth": float(truth_value),
            "fit": fit_value,
            "abs_err": abs_err,
            "rel_err": rel_err,
        }
    return out


def residual_rms(predicted: NDArray[np.float64], measured: NDArray[np.float64]) -> float:
    """残差 RMS（nm）。"""
    diff = np.asarray(predicted, dtype=np.float64) - np.asarray(measured, dtype=np.float64)
    return float(np.sqrt(np.mean(diff**2)))


@dataclass(frozen=True)
class MetricsResult:
    coefficient_recovery: dict[str, dict[str, float]]
    max_rel_err: float
    max_rel_err_field: float
    residual_rms_nm: float
    residual_rms_truth_nm: float
    residual_agreement_nm: float
    pass_A02: bool
    pass_A03: bool
    a02_applicable: bool
    spurious_max_abs_nm: float
    notes: list[str] = field(default_factory=list)
    version: int = SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "version": self.version,
            "coefficient_recovery": {
                k: {kk: (None if np.isnan(vv) else round(vv, 8)) for kk, vv in v.items()}
                for k, v in self.coefficient_recovery.items()
            },
            "max_rel_err": None if np.isnan(self.max_rel_err) else round(self.max_rel_err, 8),
            "max_rel_err_field": (
                None if np.isnan(self.max_rel_err_field) else round(self.max_rel_err_field, 8)
            ),
            "residual_rms_nm": round(self.residual_rms_nm, 8),
            "residual_rms_truth_nm": round(self.residual_rms_truth_nm, 8),
            "residual_agreement_nm": round(self.residual_agreement_nm, 8),
            "spurious_max_abs_nm": round(self.spurious_max_abs_nm, 8),
            "pass_A02": self.pass_A02,
            "pass_A03": self.pass_A03,
            "a02_applicable": self.a02_applicable,
            "notes": list(self.notes),
        }


def _as_mapping(fit: Any, key: str) -> dict[str, float]:
    """兼容 `FitResult` 与 `fit_result.json` 读出的 dict。"""
    raw = getattr(fit, key, None)
    if raw is None:
        if not isinstance(fit, Mapping):
            raise TypeError(f"fit 既不是 FitResult 也不是映射：{type(fit).__name__}")
        raw = fit[key]
    return {str(k): float(v) for k, v in dict(raw).items()}


def _as_scalar(fit: Any, key: str) -> float:
    value = getattr(fit, key, None)
    if value is None:
        if not isinstance(fit, Mapping):
            raise TypeError(f"fit 既不是 FitResult 也不是映射：{type(fit).__name__}")
        value = fit[key]
    return float(value)


def _truth_values(
    measurements: MeasurementSet, truth: Truth, channel: str
) -> NDArray[np.float64]:
    """由真值系数重建真值场，用于算"真值残差"（A-04 需要它）。"""
    if channel == "x":
        w_map, f_map = truth.wafer_terms, truth.field_terms_x
    else:
        w_map, f_map = truth.wafer_terms_y, truth.field_terms_y
    wafer_terms = list(w_map)
    field_terms = list(f_map)
    zw = ZernikeBasis(truth.max_order).values(
        measurements.rho, measurements.theta_rad, wafer_terms
    )
    xi, eta = measurements.xi_eta()
    cf = CartesianBasis(_field_order_of(field_terms)).values(xi, eta, field_terms)
    a = np.array([w_map[t] for t in wafer_terms], dtype=np.float64)
    b = np.array([f_map[t] for t in field_terms], dtype=np.float64)
    return np.asarray(zw @ a + cf @ b, dtype=np.float64)


def _field_order_of(field_terms: list[str]) -> int:
    if not field_terms:
        return 0
    return max(int(t.split(",")[0]) + int(t.split(",")[1]) for t in field_terms)


def compare(fit: Any, truth: Truth, measurements: MeasurementSet) -> MetricsResult:
    """对比拟合结果与真值，输出 A-02 / A-03 的判定。

    `fit` 可以是 `FitResult` 也可以是 `fit_result.json` 读出的 dict。
    """
    fit_x = _as_mapping(fit, "wafer_coeffs_x")
    fit_y = _as_mapping(fit, "wafer_coeffs_y")
    fit_fx = _as_mapping(fit, "field_coeffs_x")
    fit_fy = _as_mapping(fit, "field_coeffs_y")
    del fit_fy  # Phase 1 的 A-02/A-03 只看 x 通道（y 通道同构，系数一并写进 metrics）

    wafer_rec = coefficient_recovery(fit_x, truth.wafer_terms)
    field_rec = coefficient_recovery(fit_fx, truth.field_terms_x)
    recovery = dict(wafer_rec)
    recovery.update({f"field {k}": v for k, v in field_rec.items()})
    if truth.wafer_terms_y:
        recovery.update(
            {f"y {k}": v for k, v in coefficient_recovery(fit_y, truth.wafer_terms_y).items()}
        )

    notes: list[str] = []

    def _max_rel(rec: dict[str, dict[str, float]], pairs: list[tuple[str, float]]) -> float:
        values = [rec[t]["rel_err"] for t, v in pairs if v != 0 and t in rec]
        return float(np.max(values)) if values else float("nan")

    wafer_pairs = list(truth.wafer_terms.items())
    max_rel = _max_rel(wafer_rec, wafer_pairs)
    max_rel_field = _max_rel(field_rec, list(truth.field_terms_x.items()))

    # A-03：只看真值里 j ≤ 6 的非零 wafer 项（"前 6 阶主要项"）
    main_pairs = [(t, v) for t, v in wafer_pairs if int(t[1:]) <= A03_MAX_ORDER]
    max_rel_main = _max_rel(wafer_rec, main_pairs)
    pass_A03 = bool(np.isfinite(max_rel_main) and max_rel_main < A03_REL_TOL)
    if not main_pairs:
        notes.append("真值里没有 Z1..Z6 的非零项，A-03 无从判定")

    n_points = int(_as_scalar(fit, "n_points"))
    n_params = int(_as_scalar(fit, "n_params"))
    noise_sigma = float(truth.noise_sigma_nm)
    a02_applicable = noise_sigma == 0.0
    pass_A02 = bool(
        a02_applicable
        and np.isfinite(max_rel)
        and max_rel < A02_REL_TOL
        and n_points >= A02_POINT_RATIO * n_params
    )
    if not a02_applicable:
        notes.append(
            f"A-02 不适用：数据带噪（σ={noise_sigma}nm），而该验收的前提是「无噪声」。"
            "带噪数据的判定看 A-03 与残差一致性 —— 否则一条 1% 的阈值会被噪声蒙对"
        )
    elif n_points < A02_POINT_RATIO * n_params:
        notes.append(
            f"A-02 前置不满足：点数 {n_points} < {A02_POINT_RATIO}×参数数 {n_params}"
        )
    elif not pass_A02:
        notes.append(f"A-02 未过：最大相对误差 {max_rel:.4f} ≥ {A02_REL_TOL}")

    # 真值里没给、拟合却给出显著值的项
    truth_terms = set(truth.wafer_terms)
    spurious = {t: abs(v) for t, v in fit_x.items() if t not in truth_terms and abs(v) > SPURIOUS_TOL_NM}
    spurious_max = max(spurious.values()) if spurious else 0.0
    if spurious:
        notes.append(
            "以下项真值为 0 但拟合给出显著值："
            + ", ".join(f"{t}={v:.3f}nm" for t, v in sorted(spurious.items()))
        )

    resid_fit = float(_as_scalar(fit, "residual_rms_nm"))
    truth_x = _truth_values(measurements, truth, "x")
    truth_y = _truth_values(measurements, truth, "y")
    rms_truth = float(
        np.sqrt(
            (
                residual_rms(truth_x, measurements.dx_nm) ** 2
                + residual_rms(truth_y, measurements.dy_nm) ** 2
            )
            / 2.0
        )
    )

    return MetricsResult(
        coefficient_recovery=recovery,
        max_rel_err=max_rel,
        max_rel_err_field=max_rel_field,
        residual_rms_nm=resid_fit,
        residual_rms_truth_nm=rms_truth,
        residual_agreement_nm=abs(resid_fit - rms_truth),
        pass_A02=pass_A02,
        pass_A03=pass_A03,
        a02_applicable=a02_applicable,
        spurious_max_abs_nm=float(spurious_max),
        notes=notes,
    )
