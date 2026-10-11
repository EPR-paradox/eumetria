"""正向生成：采样点表与真值指纹（§6.2.4 的第 2–7 步）。

顺序不可交换：先抽 field、再取 mark、算全场坐标、转极坐标、**过滤掉边缘排除区内的
mark**，最后才求指纹。把过滤放到最后会让 Zernike 在 ρ>1 上被调用（§9.2 直接报错）。

与详设的一处有意偏离：`make_fingerprint` 的入参从 `(layout, marks_per_field, …)`
改为 `(points: SamplePoints, …)`。理由：写 `measurements.csv` 需要每个点的
field_id / 场内坐标 / ρ / θ，若 make_fingerprint 内部自己重算一遍采样点，采样逻辑
就会存在两份实现，迟早漂移。现在 `sample_points` 是点表的唯一来源，
`make_fingerprint` 只负责在给定点上求值。
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

from ..basis import Basis
from ..geometry import Field, field_marks, normalize_field, to_polar

Arr: TypeAlias = NDArray[np.float64]
IntArr: TypeAlias = NDArray[np.int64]
BoolArr: TypeAlias = NDArray[np.bool_]

# 空间分块的网格数（5×5 = 25）。regular / stratified 两种确定性采样共用。
SAMPLE_BLOCKS = 5
_RHO_TOL = 1e-12


@dataclass(frozen=True)
class SamplePoints:
    """一次采样的全部点，连同它们的来历。这是 measurements.csv 的唯一来源。"""

    wafer_id: str
    field_id: IntArr  # 沿用 build_layout 的 field_id（保留物理来历，不重编号）
    point_id: IntArr
    x_mm: Arr
    y_mm: Arr
    xf_mm: Arr
    yf_mm: Arr
    rho: Arr
    theta_rad: Arr
    is_edge: BoolArr  # 该点所属 field 是边缘场（见 sample_points 文档）
    field_width_mm: float
    field_height_mm: float

    def __len__(self) -> int:
        return int(self.x_mm.size)

    def xi_eta(self) -> tuple[Arr, Arr]:
        """场内归一化坐标 (ξ, η) ∈ [-1,1]^2，field-level 基函数的自变量。"""
        return normalize_field(self.xf_mm, self.yf_mm, self.field_width_mm, self.field_height_mm)


def sample_fields(
    layout: list[Field],
    n_fields: int,
    mode: str,
    rng: np.random.Generator,
    blocks: int = SAMPLE_BLOCKS,
) -> list[Field]:
    """从 layout 里抽 n_fields 个 field。

    mode='regular'    空间分块，每块取距块心最近的 field（确定性）
    mode='stratified' 空间分块，每块在块内随机取一个（由 rng 决定，同 seed 同结果）
    mode='random'     全局无放回随机

    为什么 `regular` 不是「按列表索引等距抽取」：layout 的 field 列表是网格遍历序
    （先 y 后 x），等距索引会把抽到的点全部赶到某一条带上（实测 25 点等距抽出的
    ρ 全部 > 0.75），完全失去空间代表性。分块法保证每一块都有代表。
    """
    if not layout:
        raise ValueError("layout 为空，无法采样")
    if not 1 <= n_fields <= len(layout):
        raise ValueError(f"n_fields 必须在 [1, {len(layout)}] 内，得到 {n_fields}")
    if mode not in ("regular", "stratified", "random"):
        raise ValueError(f"未知的采样模式 {mode!r}")
    if blocks < 1:
        raise ValueError(f"blocks 必须 ≥ 1，得到 {blocks}")

    if mode == "random":
        idx = rng.choice(len(layout), size=n_fields, replace=False)
        return [layout[int(i)] for i in sorted(idx)]

    cx = np.array([f.cx_mm for f in layout], dtype=np.float64)
    cy = np.array([f.cy_mm for f in layout], dtype=np.float64)

    # 在 field 中心的包围盒上做 blocks×blocks 等宽分块
    x_edges = np.linspace(cx.min(), cx.max(), blocks + 1)
    y_edges = np.linspace(cy.min(), cy.max(), blocks + 1)
    block_of = np.full(len(layout), -1, dtype=np.int64)
    block_centers: list[tuple[float, float]] = []
    for jy in range(blocks):
        for ix in range(blocks):
            bid = jy * blocks + ix
            lo_x, hi_x = x_edges[ix], x_edges[ix + 1]
            lo_y, hi_y = y_edges[jy], y_edges[jy + 1]
            block_centers.append(((lo_x + hi_x) / 2.0, (lo_y + hi_y) / 2.0))
            # 边界的归属：左闭右闭，最外一层收进最后一格（避免浮点把边界点丢掉）
            in_block = (cx >= lo_x) & (cx <= hi_x) & (cy >= lo_y) & (cy <= hi_y)
            block_of[in_block & (block_of < 0)] = bid

    chosen: list[int] = []
    taken = np.zeros(len(layout), dtype=bool)
    while len(chosen) < n_fields:
        progressed = False
        for bid, (bx, by) in enumerate(block_centers):
            if len(chosen) >= n_fields:
                break
            members = np.flatnonzero((block_of == bid) & ~taken)
            if members.size == 0:
                continue
            dist = np.hypot(cx[members] - bx, cy[members] - by)
            if mode == "stratified":
                pick = int(rng.choice(members))
            else:
                # 距块心最近者；平局按 field_id 小者优先（确定性）
                order = np.lexsort((members, dist))
                pick = int(members[order[0]])
            taken[pick] = True
            chosen.append(pick)
            progressed = True
        if not progressed:  # 理论上不可达：n_fields ≤ len(layout) 且每轮必有空位可补
            break
    return [layout[i] for i in chosen[:n_fields]]


def sample_points(
    fields: list[Field],
    marks_per_field: int,
    wafer_radius_mm: float,
    edge_exclusion_mm: float,
) -> SamplePoints:
    """把抽中的 field 展开成 mark 点表，并按 edge exclusion 过滤。

    `is_edge` 的语义（与详设 §6.2 的写法有别，已按实现修订）：**该点所属 field 是
    边缘场**（`Field.is_edge`，即该场有 mark 落在排除区内）。落在排除区内的 mark 本身
    已被剔除，所以「该点自己在排除区内」恒为 False，列会失去信息量。
    """
    if not fields:
        raise ValueError("fields 为空，无法展开 mark")
    usable_radius_mm = wafer_radius_mm - edge_exclusion_mm
    if usable_radius_mm <= 0:
        raise ValueError(
            f"可用半径非正：R={wafer_radius_mm}，edge_exclusion={edge_exclusion_mm}"
        )
    usable_rho = usable_radius_mm / wafer_radius_mm
    return _collect(fields, marks_per_field, wafer_radius_mm, usable_rho)


def _collect(
    fields: list[Field],
    marks_per_field: int,
    wafer_radius_mm: float,
    usable_rho: float,
) -> SamplePoints:
    xs: list[float] = []
    ys: list[float] = []
    xfs: list[float] = []
    yfs: list[float] = []
    fids: list[int] = []
    edges: list[bool] = []
    for field in fields:
        for mark in field_marks(field, marks_per_field):
            xs.append(mark.x_mm)
            ys.append(mark.y_mm)
            xfs.append(mark.xf_mm)
            yfs.append(mark.yf_mm)
            fids.append(field.field_id)
            edges.append(field.is_edge)

    x_arr = np.array(xs, dtype=np.float64)
    y_arr = np.array(ys, dtype=np.float64)
    rho_all, theta_all = to_polar(x_arr, y_arr, wafer_radius_mm)
    keep = rho_all <= usable_rho + _RHO_TOL
    if not np.any(keep):
        raise ValueError("过滤后没有任何可用 mark；检查 edge exclusion 与 mark 内缩量")

    xf_arr = np.array(xfs, dtype=np.float64)[keep]
    yf_arr = np.array(yfs, dtype=np.float64)[keep]
    n = int(keep.sum())
    return SamplePoints(
        wafer_id="",
        field_id=np.array(fids, dtype=np.int64)[keep],
        point_id=np.arange(n, dtype=np.int64),
        x_mm=x_arr[keep],
        y_mm=y_arr[keep],
        xf_mm=xf_arr,
        yf_mm=yf_arr,
        rho=rho_all[keep],
        theta_rad=theta_all[keep],
        is_edge=np.array(edges, dtype=np.bool_)[keep],
        field_width_mm=fields[0].width_mm,
        field_height_mm=fields[0].height_mm,
    )


def _sorted_terms(coeffs: dict[str, float]) -> list[str]:
    """按数值含义排序项名（§9.5：禁止依赖 dict 迭代顺序影响数值结果）。

    字符串序会把 "Z11" 排到 "Z4" 前面 —— 那不只是顺序问题，浮点求和次序变了，
    结果就不再是 bit 级可复现的。
    """

    def key(term: str) -> tuple[int, int, int]:
        if term.startswith("Z"):
            return (0, int(term[1:]), 0)
        p_str, q_str = term.split(",")
        return (1, int(p_str), int(q_str))

    return sorted(coeffs, key=key)


def _eval(basis: Basis, coords: tuple[Arr, Arr], coeffs: dict[str, float], n_points: int) -> Arr:
    if not coeffs:
        return np.zeros(n_points, dtype=np.float64)
    terms = _sorted_terms(coeffs)
    design = basis.values(*coords, terms=terms)
    return design @ np.array([coeffs[t] for t in terms], dtype=np.float64)


def make_fingerprint(
    points: SamplePoints,
    basis: Basis,
    wafer_terms_x: dict[str, float],
    wafer_terms_y: dict[str, float],
    field_basis: Basis,
    field_terms_x: dict[str, float],
    field_terms_y: dict[str, float],
) -> Arr:
    """返回未加噪的真实 (dx, dy)，形状 (n_points, 2)，单位 nm。

        dx_true = Σ_j a_j^x · Z_j(ρ,θ) + Σ_{p,q} b_{pq}^x · ξ^p η^q

    wafer 项与 field 项**相加**：同一张片上每个场的场内畸变相同、全片畸变随位置变化，
    两者在物理上可分离，正好对应两级拟合的两层。
    """
    n = len(points)
    if n == 0:
        raise ValueError("points 为空")
    xi, eta = points.xi_eta()
    wafer_coords = (points.rho, points.theta_rad)
    field_coords = (xi, eta)

    dx = _eval(basis, wafer_coords, wafer_terms_x, n) + _eval(
        field_basis, field_coords, field_terms_x, n
    )
    dy = _eval(basis, wafer_coords, wafer_terms_y, n) + _eval(
        field_basis, field_coords, field_terms_y, n
    )
    return np.column_stack([dx, dy])
