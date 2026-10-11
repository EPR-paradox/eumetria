# eumetria

光刻**量测**侧的 wafer 场分析工具：把一片硅片上的量测数据分解为可解释、可判定、可诊断的成分。

主线是 **overlay**（场分解 → 可校正性 → 采样设计 → 诊断）；focus 复用同一套框架以证明框架可复用；CDU 只做 CD→(dose, focus) 反演一条。

`eumetria` = εὐμετρία：匀称、比例恰当的度量。

## 当前状态

**Phase 1 进行中。** 已完成实现顺序检查表的前 3 步：

| 模块 | 状态 |
|---|---|
| `geometry.coordinates`（`to_polar` / `normalize_field`） | 已实现 |
| `geometry.layout`（`build_layout` / `field_marks`） | 已实现，4 点与 8 点密度（19 点待定，见详设 U-10） |
| `basis.ZernikeBasis`（Noll 索引，单位方差归一化） | 已实现；第一个里程碑：设计矩阵最小二乘可恢复已知系数 |
| `synth` / `fit` / `metrics` / `cli` / `report` | 未开始 |

79 个测试全绿，其中包含与 `zernpy` 的数值对照（j=1..36 最大偏差 4.4e-16）与单位圆正交性
（j≤28 的 406 个内积最大偏离 3.3e-14）。

- 详细设计：[`docs/software-design-v0.1.md`](docs/software-design-v0.1.md)（Draft v0.2）
  - 验收条件 A-01..A-11、技术问题 Q-01..Q-10、交付切片 Phase 1..4、实现顺序检查表

## 开发

```bash
# 环境（本机 python3.10-venv 缺 ensurepip，用 uv 建 venv）
uv venv --python /usr/bin/python3.10 .venv
uv pip install --python .venv/bin/python -i https://pypi.tuna.tsinghua.edu.cn/simple \
    numpy scipy pyyaml pydantic matplotlib pytest zernpy

# 测试
.venv/bin/python -m pytest
```

## 非目标

不做缺陷检测/分类（属 inspection）、不接真实 fab 数据、不做 Web 服务、不引入 C++/native 层。
