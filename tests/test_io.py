"""io 层测试：schema 校验、读写往返、§6.2.5 的四条强制校验。"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
from pydantic import ValidationError

from eumetria.errors import (
    InvalidMeasurementError,
    SchemaVersionError,
    SchemaValidationError,
)
from eumetria.io import (
    MEASUREMENT_COLUMNS,
    MeasurementSet,
    load_measurements,
    load_scenario,
    load_truth,
    write_measurements,
    write_truth,
)
from eumetria.synth import generate

SCENARIO_PATH = Path("/home/j/myprojects/eumetria/scenarios/base-300mm.yaml")


@pytest.fixture(scope="module")
def scenario():
    return load_scenario(SCENARIO_PATH)


@pytest.fixture(scope="module")
def run(scenario):
    return generate(scenario)


@pytest.fixture()
def csv_path(tmp_path, run) -> Path:
    return write_measurements(tmp_path / "out" / "measurements.csv", run.measurements)


def _load(path: Path, scenario) -> MeasurementSet:
    return load_measurements(
        path,
        wafer_radius_mm=scenario.wafer_radius_mm,
        field_width_mm=scenario.field.width_mm,
        field_height_mm=scenario.field.height_mm,
    )


def _rewrite(path: Path, mutate) -> None:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    mutate(rows)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(MEASUREMENT_COLUMNS))
        writer.writeheader()
        writer.writerows(rows)


# ------------------------------------------------------------- scenario


def test_scenario_defaults_match_document(scenario):
    assert scenario.wafer.diameter_mm == 300
    assert scenario.wafer.edge_exclusion_mm == 3.0
    assert scenario.field.marks_per_field == 8
    assert scenario.sampling.fields == 25
    assert scenario.sampling.selection == "regular"
    assert scenario.sampling.tmu_nm == 0.1
    assert scenario.tis.amplitude_nm == 0.3 and scenario.tis.angle_deg == 30.0
    assert scenario.fingerprint.max_order == 12
    assert scenario.fingerprint.terms == {
        "Z4": 2.5,
        "Z5": -1.8,
        "Z6": 1.1,
        "Z7": 0.9,
        "Z11": -2.2,
    }
    assert scenario.fingerprint.field_terms == {"1,0": 1.2, "0,1": -0.8, "2,0": 0.5}
    assert scenario.wafer_radius_mm == 150.0


def test_scenario_rejects_extra_key(tmp_path, scenario):
    text = SCENARIO_PATH.read_text(encoding="utf-8").replace(
        "wafer_id: W0001", "wafer_id: W0001\nseedx: 1"
    )
    p = tmp_path / "bad.yaml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ValidationError):
        load_scenario(p)


def test_scenario_rejects_term_beyond_max_order(tmp_path):
    text = SCENARIO_PATH.read_text(encoding="utf-8").replace("Z11: -2.2", "Z20: -2.2")
    p = tmp_path / "bad.yaml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ValidationError, match="max_order"):
        load_scenario(p)


def test_scenario_rejects_unsupported_mark_density(tmp_path):
    text = SCENARIO_PATH.read_text(encoding="utf-8").replace(
        "marks_per_field: 8", "marks_per_field: 19"
    )
    p = tmp_path / "bad.yaml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(ValidationError, match="U-10"):
        load_scenario(p)


def test_scenario_rejects_bad_version(tmp_path):
    text = SCENARIO_PATH.read_text(encoding="utf-8").replace("version: 1", "version: 2", 1)
    p = tmp_path / "bad.yaml"
    p.write_text(text, encoding="utf-8")
    with pytest.raises(SchemaVersionError):
        load_scenario(p)


# --------------------------------------------------------- 读写往返


def test_measurements_roundtrip(run, csv_path, scenario):
    back = _load(csv_path, scenario)
    assert len(back) == len(run.measurements)
    for name in ("x_mm", "y_mm", "xf_mm", "yf_mm", "rho", "theta_rad", "dx_nm", "dy_nm"):
        assert np.allclose(getattr(back, name), getattr(run.measurements, name), rtol=0, atol=1e-11)
    assert np.array_equal(back.field_id, run.measurements.field_id)
    assert np.array_equal(back.is_edge, run.measurements.is_edge)
    for name in ("dx_nm_0", "dy_nm_0", "dx_nm_180", "dy_nm_180"):
        got, want = getattr(back, name), getattr(run.measurements, name)
        assert (got is None) == (want is None)
        if want is not None:
            assert np.allclose(got, want, rtol=0, atol=1e-11)


def test_csv_header_is_the_documented_19_columns(csv_path):
    with csv_path.open(encoding="utf-8") as fh:
        header = next(csv.reader(fh))
    assert header == list(MEASUREMENT_COLUMNS)
    assert len(header) == 19


def test_is_edge_column_is_not_all_false(run):
    """边缘场的点必须被标出来，否则这一列没有信息量（§6.2 的语义修订）。"""
    assert run.measurements.is_edge.any()
    # 被标 edge 的场 = 有 mark 落进排除区的场 = 恰好是那 4 个各剩 5 个可用 mark 的场
    edge_fields = sorted(set(run.measurements.field_id[run.measurements.is_edge]))
    assert len(edge_fields) == 4
    counts = {int(f): int((run.measurements.field_id == f).sum()) for f in edge_fields}
    assert set(counts.values()) == {5}, f"边缘场各自的可用 mark 数：{counts}"
    assert int(run.measurements.is_edge.sum()) == 20


# --------------------------------------------- §6.2.5 的四条强制校验


def test_reader_rejects_tampered_rho(csv_path, scenario):
    def mutate(rows):
        rows[3]["rho"] = "0.5"
    _rewrite(csv_path, mutate)
    with pytest.raises(InvalidMeasurementError, match="rho 与重算值不一致"):
        _load(csv_path, scenario)


def test_reader_rejects_rho_out_of_range(csv_path, scenario):
    """把坐标推到片外（rho > 1）：第 2 条校验必须拦下。"""
    def mutate(rows):
        for r in rows[:1]:
            r["x_mm"] = "200.0"
            r["y_mm"] = "0.0"
            # 同时把 rho/theta 改成"自洽"的值，确保被拦是因为越界而不是第 1 条
            r["rho"] = f"{200.0 / 150.0:.12g}"
            r["theta_rad"] = "0"
    _rewrite(csv_path, mutate)
    with pytest.raises(InvalidMeasurementError, match=r"rho > 1"):
        _load(csv_path, scenario)


def test_reader_rejects_inconsistent_tis_combination(csv_path, scenario):
    def mutate(rows):
        rows[0]["dx_nm"] = "999.0"
    _rewrite(csv_path, mutate)
    with pytest.raises(InvalidMeasurementError, match="dx_nm 与"):
        _load(csv_path, scenario)


def test_reader_rejects_duplicate_mark_in_field(csv_path, scenario):
    def mutate(rows):
        # 让两行的场内坐标相同（同一 field 内）
        rows[1]["xf_mm"], rows[1]["yf_mm"] = rows[0]["xf_mm"], rows[0]["yf_mm"]
    _rewrite(csv_path, mutate)
    with pytest.raises(InvalidMeasurementError, match="重复"):
        _load(csv_path, scenario)


def test_reader_rejects_missing_or_unknown_column(csv_path, scenario):
    with csv_path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        del r["dz_nm"]
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[c for c in MEASUREMENT_COLUMNS if c != "dz_nm"])
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(SchemaValidationError, match="缺少列"):
        _load(csv_path, scenario)


def test_reader_rejects_bad_bool_token(csv_path, scenario):
    def mutate(rows):
        rows[0]["is_edge"] = "yes"
    _rewrite(csv_path, mutate)
    with pytest.raises(SchemaValidationError, match="is_edge"):
        _load(csv_path, scenario)


# ------------------------------------------------------------- truth


def test_truth_roundtrip(tmp_path, run):
    p = write_truth(tmp_path / "truth.json", run.truth)
    back = load_truth(p)
    assert back == run.truth
    assert back.tis_vector_nm == pytest.approx((0.2598076, 0.15), abs=1e-6)
    assert back.noise_sigma_nm == 0.1
    assert back.eumetria_version == run.truth.eumetria_version


def test_truth_records_independent_x_y_slots(run):
    """truth.json 必须同时写出 x/y 两组槽位，将来 scenario 支持 y 独立系数时不用改契约。"""
    payload = run.truth.model_dump()
    for key in ("wafer_terms", "wafer_terms_y", "field_terms_x", "field_terms_y"):
        assert key in payload


# --------------------------------------------------- 生成器的自洽性


def test_tis_estimate_recovers_the_injected_bias(run):
    """从 188 个点的两次测量反推 TIS，误差应远小于量级 —— 证明双测组合真的消掉了偏移。"""
    dx_true, dy_true = run.truth.tis_vector_nm
    assert run.tis_estimate_nm[0] == pytest.approx(dx_true, abs=0.02)
    assert run.tis_estimate_nm[1] == pytest.approx(dy_true, abs=0.02)


def test_generate_is_bit_reproducible(scenario):
    """A-05 的前置条件：同一 seed 两次生成，文件级内容一致。"""
    a = generate(scenario)
    b = generate(scenario)
    assert np.array_equal(a.measurements.dx_nm, b.measurements.dx_nm)
    assert np.array_equal(a.measurements.dy_nm, b.measurements.dy_nm)
    assert np.array_equal(a.truth_values, b.truth_values)


def test_measurement_noise_level_matches_sigma_over_sqrt2(scenario, run):
    """组合后的残差应当是 σ/√2 ≈ 0.0707nm（§3.4），不是 σ。"""
    residual = run.measurements.dx_nm - run.truth_values[:, 0]
    assert residual.std() == pytest.approx(0.1 / np.sqrt(2.0), rel=0.15)
    assert abs(residual.mean()) < 0.02
