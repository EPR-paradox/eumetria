# eumetria

光刻**量测**侧的 wafer 场分析工具：把一片硅片上的量测数据分解为可解释、可判定、可诊断的成分。

主线是 **overlay**（场分解 → 可校正性 → 采样设计 → 诊断）；focus 复用同一套框架以证明框架可复用；CDU 只做 CD→(dose, focus) 反演一条。

`eumetria` = εὐμετρία：匀称、比例恰当的度量。

## 当前状态

**Phase 1 闭环（`run` 一条命令跑通 synth → fit → score）。**

| 模块 | 状态 |
|---|---|
| `geometry.coordinates` / `geometry.layout` | 已实现（76 口径B / 64 口径A；4 点与 8 点 mark） |
| `basis.ZernikeBasis` / `CartesianBasis` | 已实现（Noll 索引；场内基不含常数项，D-09） |
| `synth`（采样 / 指纹 / TIS / 噪声 / 流水线） | 已实现（默认场景 188 点，同 seed bit 级可复现） |
| `io`（scenario / measurements / truth 读写） | 已实现（§6.2.5 四条强制校验） |
| `fit`（联合最小二乘 + 可辨识性诊断） | 已实现（A-03/A-04 通过；CPE 阶数上限有告警） |
| `metrics`（系数恢复 / 残差一致性） | 已实现（唯一读 truth.json 的模块） |
| `cli`（`synth` / `fit` / `score` / `run`） | 已实现（退出码 0/1/2/3，`--json` 走 stdout） |
| `report`（图 + HTML） | 未开始（Phase 4） |

148 个测试全绿。默认场景实测：188 点 / 17 参数 / 条件数 7.60 / 残差 RMS 0.0646nm；
A-03（前 6 阶主项 < 10%）通过，A-04（残差一致性 < 0.01nm）通过。

- 进度：[`docs/progress.md`](docs/progress.md) —— **现在到哪 / 明天从哪接 / 什么还悬着**
- 详细设计：[`docs/software-design-v0.1.md`](docs/software-design-v0.1.md)（Draft v0.2）
  - 验收条件 A-01..A-11、技术问题 Q-01..Q-10、交付切片 Phase 1..4、实现顺序检查表
- 实验记录：[`docs/experiments.md`](docs/experiments.md) —— **数字的唯一出处**，
  每节的表格都由 `experiments/` 下的脚本产出
- 设计问答：[`docs/design-faq.md`](docs/design-faq.md) —— 为什么这么选、证据在哪、边界在哪

## 开发

```bash
# 环境（本机 python3.10-venv 缺 ensurepip，用 uv 建 venv）
uv venv --python /usr/bin/python3.10 .venv
uv pip install --python .venv/bin/python -i https://pypi.tuna.tsinghua.edu.cn/simple \
    numpy scipy pyyaml pydantic matplotlib pytest zernpy
uv pip install --python .venv/bin/python -e .     # 装出 eumetria 命令入口

# 跑整条流水线（synth → fit → score → report），产物在 out/base/
.venv/bin/python main.py                       # 零参数，PyCharm 里直接点绿三角
.venv/bin/eumetria run -c scenarios/base-300mm.yaml -o out/base   # 等价的命令行形式
.venv/bin/python -c "from eumetria.pipeline import run_all; run_all('scenarios/base-300mm.yaml','out/base')"

# 测试
.venv/bin/python -m pytest

# 复现实验（docs/experiments.md 里的数字都出自这些脚本）
.venv/bin/python experiments/q10_estimator_comparison.py
```

## 非目标

不做缺陷检测/分类（属 inspection）、不接真实 fab 数据、不做 Web 服务、不引入 C++/native 层。
