"""两级模型的联合最小二乘（§7.4、D-09）。

模型（§3.1）：

    dx = Σ_j a_j^x·Z_j(ρ,θ)  +  Σ_{p,q} b_pq^x·ξ^p η^q  +  ε
    dy = 同构，x/y 各有独立系数集

**做法：一次联合求解** —— 把两块基拼成一个设计矩阵

    X = [ Z(ρ,θ) | C(ξ,η) ]     形状 (n_points, n_wafer_terms + n_field_terms)

每个通道各做一次 `np.linalg.lstsq(X, y)`，再按列切开成 `wafer_coeffs_*` / `field_coeffs_*`。
分层是可解释性的需求（wafer 系数给 APC、场内系数给 i-HOPC/CPE），切开系数向量就能分层，
**不需要用有偏的顺序解去换**（D-09；证据见 docs/experiments.md → Q-10）。

守卫（§3.1.1、§10.1）：
- `n_points < n_params` 或设计矩阵降秩 → `InsufficientPointsError`，**不给解**
- 这两条必须显式检查：`np.linalg.lstsq` 在欠定系统上不报错，会静默返回最小范数解
- 条件数 > 1e6、项间相关系数 > 0.9 → 记入 `identifiability.warnings`（软警告，继续跑）
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from ..basis import CartesianBasis, ZernikeBasis
from ..errors import InsufficientPointsError, SchemaValidationError
from ..io.schemas import MeasurementSet
from .result import FieldDiagnostics, FitResult, IdentifiabilityInfo

# §9.3 的阈值
COND_WARN = 1e6
CORR_WARN = 0.9
# 场内基不含常数项（D-09）：常量归 wafer-level 的 Tx = Z1
FIELD_TERM_SKIP = "0,0"
# 逐场可辨识性最多检查到几阶（CPE 阶数上限的搜索范围）
MAX_FIELD_ORDER_CHECK = 3


@dataclass(frozen=True)
class FitConfig:
    """拟合配置（§7.4）。`aic` 属 Phase 2，本版只实现 `none` / `fixed`。"""

    basis: str = "zernike"
    max_order: int = 12
    selection: Literal["none", "fixed", "aic"] = "none"
    fixed_terms: list[str] | None = None
    field_order: int = 2

    def __post_init__(self) -> None:
        if self.basis != "zernike":
            raise SchemaValidationError(f"本期只实现 zernike 基，得到 {self.basis!r}")
        if self.selection not in ("none", "fixed", "aic"):
            raise SchemaValidationError(f"未知的 selection={self.selection!r}")
        if self.selection == "aic":
            raise NotImplementedError(
                "selection='aic' 属 Phase 2（前向逐步 + AIC，见 Q-02）；"
                "Phase 1 用 'none' 或 'fixed'"
            )
        if self.selection == "fixed" and not self.fixed_terms:
            raise SchemaValidationError("selection='fixed' 必须给 fixed_terms")
        if self.max_order < 1:
            raise SchemaValidationError(f"max_order 必须 ≥ 1，得到 {self.max_order}")
        if self.field_order < 1:
            raise SchemaValidationError(f"field_order 必须 ≥ 1，得到 {self.field_order}")
        if self.selection == "fixed" and self.fixed_terms:
            # 配置错误在构造时就报，不要拖到拟合时才炸
            all_terms = set(ZernikeBasis(self.max_order).terms())
            unknown = [t for t in self.fixed_terms if t not in all_terms]
            if unknown:
                raise SchemaValidationError(
                    f"fixed_terms 含超出 max_order={self.max_order} 的项：{unknown}"
                )


def field_terms_for(order: int) -> list[str]:
    """场内项集：order 阶幂基，去掉常数项（D-09）。"""
    return [t for t in CartesianBasis(order).terms() if t != FIELD_TERM_SKIP]


def _select_wafer_terms(config: FitConfig) -> list[str]:
    """按 selection 决定 wafer 项集。'none' = 全部 Z1..Z<max_order>。"""
    all_terms = ZernikeBasis(config.max_order).terms()
    if config.selection == "fixed":
        assert config.fixed_terms is not None  # __post_init__ 已保证
        unknown = [t for t in config.fixed_terms if t not in all_terms]
        if unknown:
            raise SchemaValidationError(
                f"fixed_terms 含超出 max_order={config.max_order} 的项：{unknown}"
            )
        return list(config.fixed_terms)
    return all_terms


def build_design(
    measurements: MeasurementSet, config: FitConfig
) -> tuple[NDArray[np.float64], list[str], list[str]]:
    """返回 (X, wafer_terms, field_terms)。X 的两块按 wafer → field 顺序拼接。"""
    wafer_terms = _select_wafer_terms(config)
    field_terms = field_terms_for(config.field_order)
    zw = ZernikeBasis(config.max_order).values(
        measurements.rho, measurements.theta_rad, wafer_terms
    )
    xi, eta = measurements.xi_eta()
    cf = CartesianBasis(config.field_order).values(xi, eta, field_terms)
    return np.hstack([zw, cf]), wafer_terms, field_terms


def _solve_channel(
    X: NDArray[np.float64], y: NDArray[np.float64], channel: str
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """一次最小二乘 + 显式守卫，返回 (系数, 残差)。"""
    n_points, n_params = X.shape
    if n_points < n_params:
        raise InsufficientPointsError(
            f"{channel} 通道点数不足：方程数 < 未知数", required=n_params, available=n_points
        )
    rank = int(np.linalg.matrix_rank(X))
    if rank < n_params:
        raise InsufficientPointsError(
            f"{channel} 通道设计矩阵降秩（rank={rank}）：存在不可辨识方向，拒绝给出最小范数解",
            required=n_params,
            available=rank,
        )
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return coef, y - X @ coef


def _max_identifiable_order(measurements: MeasurementSet, config: FitConfig) -> int:
    """最大的 j，使 [Z1..Zj] 单独作设计矩阵时满秩且条件数可接受。"""
    basis = ZernikeBasis(config.max_order)
    order = 0
    for j in range(1, config.max_order + 1):
        terms = [f"Z{k}" for k in range(1, j + 1)]
        design = basis.values(measurements.rho, measurements.theta_rad, terms)
        if int(np.linalg.matrix_rank(design)) < j or float(np.linalg.cond(design)) > COND_WARN:
            break
        order = j
    return order


def _correlation_warnings(X: NDArray[np.float64], terms: list[str]) -> list[str]:
    """项间相关系数 > 0.9 的项对（§9.3）。"""
    warnings: list[str] = []
    centred = X - X.mean(axis=0, keepdims=True)
    norms = np.linalg.norm(centred, axis=0)
    for i in range(X.shape[1]):
        if norms[i] == 0:
            continue
        for j in range(i + 1, X.shape[1]):
            if norms[j] == 0:
                continue
            corr = abs(float(centred[:, i] @ centred[:, j]) / (norms[i] * norms[j]))
            if corr > CORR_WARN:
                warnings.append(
                    f"{terms[i]} 与 {terms[j]} 在本采样布局下相关系数 {corr:.2f}，"
                    "建议增加点数或去掉其一"
                )
    return warnings


def _field_diagnostics(
    measurements: MeasurementSet, config: FitConfig
) -> dict[int, FieldDiagnostics]:
    """逐场的秩与条件数，以及每个场能支撑的最高场内阶数。

    逐场检查时**用含常数项的完整幂基**：在逐场拟合里常数项不与任何东西混叠
    （混叠只在跨场共享系数时才出现，见 D-09）。
    """
    xi, eta = measurements.xi_eta()
    configured = CartesianBasis(config.field_order)
    configured_terms = configured.terms()
    order_bases = {
        order: CartesianBasis(order) for order in range(1, MAX_FIELD_ORDER_CHECK + 1)
    }

    out: dict[int, FieldDiagnostics] = {}
    for field_id in sorted({int(v) for v in measurements.field_id}):
        mask = measurements.field_id == field_id
        design = configured.values(xi[mask], eta[mask], configured_terms)
        rank = int(np.linalg.matrix_rank(design))
        cond = float(np.linalg.cond(design)) if rank == len(configured_terms) else float("inf")

        # 这个场能支撑到几阶：逐阶检查秩是否满
        best = 0
        for order, basis in order_bases.items():
            terms = basis.terms()
            sub = basis.values(xi[mask], eta[mask], terms)
            if int(np.linalg.matrix_rank(sub)) == len(terms):
                best = order
            else:
                break

        out[field_id] = FieldDiagnostics(
            field_id=field_id,
            n_marks=int(mask.sum()),
            configured_order=config.field_order,
            rank=rank,
            n_params=len(configured_terms),
            condition_number=cond,
            max_identifiable_order=best,
        )
    return out


def fit_field(measurements: MeasurementSet, config: FitConfig | None = None) -> FitResult:
    """两级模型的联合最小二乘拟合。

    与详设的差异（有意）：签名里去掉 `layout` —— `MeasurementSet` 自带场中心、场尺寸
    与场内坐标，(ρ,θ) 也已经在读取时重算校验过；再传一份 layout 只会多一个能和它矛盾
    的信息源。
    """
    config = config or FitConfig()
    measurements.check_consistency()

    X, wafer_terms, field_terms = build_design(measurements, config)
    n_params = X.shape[1]

    coef_x, resid_x = _solve_channel(X, measurements.dx_nm, "dx")
    coef_y, resid_y = _solve_channel(X, measurements.dy_nm, "dy")

    rms_x = float(np.sqrt(np.mean(resid_x**2)))
    rms_y = float(np.sqrt(np.mean(resid_y**2)))
    rms = float(np.sqrt((rms_x**2 + rms_y**2) / 2.0))
    cond = float(np.linalg.cond(X))

    diagnostics = _field_diagnostics(measurements, config)
    ceiling = min(d.max_identifiable_order for d in diagnostics.values())

    warnings = _correlation_warnings(X, wafer_terms + field_terms)
    if cond > COND_WARN:
        warnings.append(
            f"设计矩阵条件数 {cond:.3g} 超过 {COND_WARN:.0e}：解对输入扰动敏感，"
            "检查是否有重复列或线性相关的项"
        )
    short = [d for d in diagnostics.values() if not d.is_identifiable]
    if short:
        ids = ", ".join(str(d.field_id) for d in short)
        warnings.append(
            f"以下场在 {config.field_order} 阶场内基下不可辨识（mark 数或秩不足）：{ids}"
            " —— 逐场校正（CPE）不得使用该阶数（§3.1.1），也不许逐场降阶"
        )
    if ceiling < config.field_order:
        warnings.append(
            f"CPE 阶数上限被最差的场限制为 {ceiling} 阶（配置要的是 {config.field_order} 阶）："
            "全局统一降档，或在采样层解决，不逐场降阶"
        )

    info = IdentifiabilityInfo(
        max_identifiable_order=_max_identifiable_order(measurements, config),
        cpe_order_ceiling=ceiling,
        warnings=warnings,
        rank_deficient=False,  # 降秩在 _solve_channel 里已抛错，能走到这里就不是降秩
    )
    n_w = len(wafer_terms)
    return FitResult(
        basis=config.basis,
        max_order=config.max_order,
        selected_terms=list(wafer_terms),
        selection_method=config.selection,
        wafer_terms=list(wafer_terms),
        field_terms=list(field_terms),
        wafer_coeffs_x=dict(zip(wafer_terms, coef_x[:n_w])),
        field_coeffs_x=dict(zip(field_terms, coef_x[n_w:])),
        wafer_coeffs_y=dict(zip(wafer_terms, coef_y[:n_w])),
        field_coeffs_y=dict(zip(field_terms, coef_y[n_w:])),
        residual_rms_nm=rms,
        residual_rms_nm_x=rms_x,
        residual_rms_nm_y=rms_y,
        n_points=len(measurements),
        n_params=n_params,
        condition_number=cond,
        identifiability=info,
        field_diagnostics=diagnostics,
    )
