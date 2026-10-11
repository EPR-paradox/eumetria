"""正向生成的完整流水线（§6.2.4 的 10 步），由一个 `Scenario` 驱动。

这是 **唯一** 能看到真值的地方：`generate()` 读 scenario、造指纹、加 TIS/噪声、组装
`MeasurementSet` 与 `Truth`。下游的 fit 只拿到 `MeasurementSet`。

噪声的随机流顺序是固定的（field → mark → TIS→ 0°/180°），所以同一个 seed
必然产生同一份文件（A-05 的 bit 级可复现）。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

from .. import __version__
from ..basis import CartesianBasis, ZernikeBasis
from ..geometry import build_layout
from ..io.schemas import MeasurementSet, Scenario, Truth
from ..io.writers import now_cst_isoformat
from .fingerprint import SamplePoints, make_fingerprint, sample_fields, sample_points
from .tis import apply_tis, combine_tis, tis_vector_from_polar

Arr: TypeAlias = NDArray[np.float64]


@dataclass(frozen=True)
class GeneratedRun:
    """一次正向生成的全部产物。`truth_values` 与 `truth` 只允许 metrics 读。"""

    points: SamplePoints
    truth_values: Arr
    measurements: MeasurementSet
    truth: Truth
    tis_estimate_nm: tuple[float, float]


def generate(scenario: Scenario, rng: np.random.Generator | None = None) -> GeneratedRun:
    """按 §6.2.4 的顺序生成一份完整数据集。"""
    if rng is None:
        rng = np.random.default_rng(scenario.seed)

    # 1. 布局（口径 B）
    layout = build_layout(
        scenario.wafer.diameter_mm,
        scenario.wafer.edge_exclusion_mm,
        scenario.field.width_mm,
        scenario.field.height_mm,
        include_partial=True,
    )

    # 2. 抽 field（regular / stratified 会消耗 rng，位置与 §6.2.4 一致）
    fields = sample_fields(
        layout, scenario.sampling.fields, scenario.sampling.selection, rng
    )

    # 3~6. 展开 mark、算全场坐标、转极坐标、过滤 edge exclusion
    points = sample_points(
        fields,
        scenario.field.marks_per_field,
        scenario.wafer_radius_mm,
        scenario.wafer.edge_exclusion_mm,
    )

    # 7. 真值指纹
    fp = scenario.fingerprint
    wafer_basis = ZernikeBasis(fp.max_order)
    field_basis = CartesianBasis(fp.field_order)
    true_values = make_fingerprint(
        points,
        wafer_basis,
        dict(fp.terms),
        dict(fp.terms),
        field_basis,
        dict(fp.field_terms),
        dict(fp.field_terms),
    )
    # 若 scenario 里显式给了 y 的独立系数则覆盖（§6.1 只列了一套 terms，
    # 这里保持 x/y 相同，与 §6.2.4 的「x/y 各有独立系数集」不冲突 ——
    # 只是本版 scenario 没有暴露 y 的入口，truth.json 里两组都写出以便对照）。
    tis_vec = tis_vector_from_polar(scenario.tis.amplitude_nm, scenario.tis.angle_deg)

    # 8~9. TIS 与噪声、0°/180° 组合
    if scenario.tis.enabled:
        meas = apply_tis(
            true_values, tis_vec, scenario.tis.calibrated, scenario.sampling.tmu_nm, rng
        )
    else:
        meas = apply_tis(
            true_values, (0.0, 0.0), False, scenario.sampling.tmu_nm, rng
        )
    if "m0" in meas:
        combined, tis_hat = combine_tis(meas["m0"], meas["m180"])
        m0, m180 = meas["m0"], meas["m180"]
    else:
        combined = meas["m"]
        tis_hat = np.zeros_like(combined)
        m0 = m180 = None

    n = len(points)
    measurements = MeasurementSet(
        wafer_id=scenario.wafer_id,
        lot_id=scenario.lot_id,
        point_id=points.point_id,
        field_id=points.field_id,
        x_mm=points.x_mm,
        y_mm=points.y_mm,
        xf_mm=points.xf_mm,
        yf_mm=points.yf_mm,
        rho=points.rho,
        theta_rad=points.theta_rad,
        dx_nm=np.ascontiguousarray(combined[:, 0]),
        dy_nm=np.ascontiguousarray(combined[:, 1]),
        is_edge=points.is_edge,
        dx_nm_0=None if m0 is None else np.ascontiguousarray(m0[:, 0]),
        dy_nm_0=None if m0 is None else np.ascontiguousarray(m0[:, 1]),
        dx_nm_180=None if m180 is None else np.ascontiguousarray(m180[:, 0]),
        dy_nm_180=None if m180 is None else np.ascontiguousarray(m180[:, 1]),
        dz_nm=None,
        cd_nm=None,
        wafer_radius_mm=scenario.wafer_radius_mm,
        field_width_mm=scenario.field.width_mm,
        field_height_mm=scenario.field.height_mm,
    )
    measurements.check_consistency()

    # 10. 真值记录（只有 metrics 会读）
    truth = Truth(
        version=scenario.version,
        seed=scenario.seed,
        wafer_id=scenario.wafer_id,
        basis=fp.basis,
        max_order=fp.max_order,
        wafer_terms=dict(fp.terms),
        wafer_terms_y=dict(fp.terms),
        field_terms_x=dict(fp.field_terms),
        field_terms_y=dict(fp.field_terms),
        tis_vector_nm=tis_vec,
        noise_sigma_nm=scenario.sampling.tmu_nm,
        generated_at=now_cst_isoformat(),
        eumetria_version=__version__,
    )
    return GeneratedRun(
        points=points,
        truth_values=true_values,
        measurements=measurements,
        truth=truth,
        tis_estimate_nm=(float(np.mean(tis_hat[:, 0])), float(np.mean(tis_hat[:, 1]))),
    )
