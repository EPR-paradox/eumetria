# 进度

**这个文件回答三个问题：现在到哪了 / 明天从哪接 / 什么还悬着。**
细节不在这里 —— 规格看 `software-design-v0.1.md`，数字看 `experiments.md`，
"为什么这么选"看 `design-faq.md`。这里只放状态与入口。

---

## 现状（2026-10-11）

**Phase 1 完成，Phase 2 未开始。** 156 个测试全绿，端到端一条命令跑通：

```bash
.venv/bin/eumetria run -c scenarios/base-300mm.yaml -o out/base
# → measurements.csv / truth.json / fit_result.json / metrics.json / report/{4 图 + report.html}
```

## 最近一次改动（2026-10-11 收尾）

新增零参数入口 `main.py`，并把流水线编排从 `cli.py` 提到 `eumetria/pipeline.py`：

| 入口 | 命令 | 说明 |
|---|---|---|
| 零参数 | `python main.py` | PyCharm 绿三角直接跑；默认场景 + `out/base` |
| 命令行 | `eumetria run -c ... -o ...` | 完整 CLI，五个子命令 |
| 库调用 | `pipeline.run_all(config, out_dir)` | notebook / 测试用；不经过 argparse |

分层的理由：编排原本埋在 `cli._cmd_run` 里，于是任何新前端都得先构造一个
`argparse.Namespace` 才能跑流水线。现在 `cli.py` 只剩参数解析 + 退出码映射，
`main.py` 只剩默认值与人类可读输出。**改流水线顺序只需改一处。**

## 下一步（明天从这里接）

**Phase 2 结论层** —— 这是项目第一次产出"我发现了什么"，而不是"我实现了什么"：

| 顺序 | 做什么 | 对应问题 | 依赖 |
|---|---|---|---|
| 1 | `fit/model_select.py`：AIC 前向逐步 + `selection='aic'` | Q-02 全项 vs AIC 精选的泛化 | `fit` 已就绪 |
| 2 | `fit/crossval.py`：K 折，**按 wafer 分组** | Q-02/Q-05 | model_select |
| 3 | `cli sweep`：点数 × 阶数网格 + `metrics/sampling_budget.py` | Q-03 采样预算曲线 | crossval |
| 4 | `fit/cascade.py`：正交性等价性实验（Zernike 基下 cascade 与 joint 系数差 < 1e-8） | Q-01 为什么用 Zernike | fit 已就绪 |
| 5 | 把 1–4 的结论写进 `experiments.md` 并在 README 汇总 | — | 全部 |

每个实验的脚本位置与状态见 `experiments.md` 的索引表（Q-01/02/03/05 目前是"待做"）。

## Phase 进度

| Phase | 内容 | 状态 |
|---|---|---|
| 1 | 闭环：geometry + basis + synth + io + fit + metrics + cli + report | **完成** |
| 2 | 结论层：model_select + crossval + sweep + 正交性等价实验 | 未开始 |
| 2.5 | 统计推断层：inference(FIM/CRLB/CI) + design(D-optimal/序贯) + sample-plan + inverse | 未开始 |
| 3 | 深度层：可校正性判定 decompose + significance + budget + tis | 未开始 |
| 4 | 门面：README 汇总 + 负面结论专节（report 已在 Phase 1 提前做完基础图） | 未开始 |

## 验收状态

| 编号 | 要求 | 状态 | 实测 |
|---|---|---|---|
| A-01 | 一条命令生成 measurements.csv + truth.json | 通过 | 188 行 / 19 列 |
| A-02 | 无噪声 + 点数 ≥ 8× 参数 → 全部系数 < 1% | **不适用** | 数据带噪；判定已按 `noise_sigma_nm==0` 加前提 |
| A-03 | 前 6 阶主要项 < 10% | 通过 | wafer 最大相对误差 0.867% |
| A-04 | 拟合残差 RMS 与真值差 < 0.01nm | 通过 | 0.06457 vs 0.06925（差 0.00468） |
| A-05 | 同 seed bit 级一致 | 通过 | 文件字节相同 |
| A-06 | 采样不足时告警并继续，不静默给错结果 | 部分 | `InsufficientPointsError` + `cpe_order_ceiling` 告警已实现；A-06 的完整断言属 Phase 2 |
| A-07 | 采样预算曲线 | 未做 | Phase 2（对应 Q-03） |
| A-08 | 全程不读真实数据、无网络 | 通过 | 无网络调用；真实数据无入口 |
| A-09 | 系数置信区间覆盖率 ∈ [93%, 97%] | 未做 | Phase 2.5 |
| A-10 | 10 万点走 QR/SVD + 报条件数 | 部分 | 条件数已报；10 万点规模属 Phase 2.5 |
| A-11 | 可校正性判定必须输出显著性 | 未做 | Phase 3（`decompose`） |

## 默认场景基线（改动后拿它对照，别凭记忆）

    点位          188（25 field × 8 mark，4 个边缘场各剩 5 mark）
    参数          17（wafer Z1..Z12 + 场内 5 项，场内基不含常数项 D-09）
    条件数        7.60
    残差 RMS      0.0646 nm（≈ σ/√2 = 0.0707，吃掉 17 个自由度后略低）
    CPE 阶数上限  1（§3.1.1 的规则在 fit 里落地并告警）
    TIS 反推      (0.2603, 0.1431) nm，真值 (0.25981, 0.15000)

## 悬着的事

**未决项（详设 §0.3）**

- **U-10** 19 点 mark 密度的定义自相矛盾（"5×5 去掉中心 6 个"既非 90° 也非 180° 对称）。
  当前实现显式拒绝该密度。4/8 点够跑完全部 Phase。
- **U-11** 边缘场 mark 是否自适应布置（往可用区内收），使每场凑够 CPE 所需点数。
  不被 Q-10 判定为瓶颈前不做。

**待量化的第一批数据点**

- 场内二阶项 `(2,0)` 在 TMU 0.1nm 下相对误差 **12.5%**，是全场最差（主项都在 1.5% 以内）。
  原因：真值小（0.5nm）+ ξ² 列在采样点上的杠杆最差（场内基条件数 6.61 由它拉起来）。
  属 **Q-04**（噪声水平 vs 各阶恢复精度），Phase 3。

**依赖与环境的坑（踩过的，别重复踩）**

- 本机 `python3.10-venv` 未装（缺 ensurepip）→ 用 `uv venv`；装包要 unset 代理走清华镜像。
- 包必须 `uv pip install -e .` 才有 `eumetria` 命令入口（pytest 靠 `pythonpath=src` 蒙过去了）。
- 环境变量 `HTTP_PROXY` 指向 Mihomo 7890，它没运行时 pip 必失败。
- 绘图中文必须显式挑 CJK 字体，否则标题画成豆腐块且只报 warning。

**自我更正记录（留着，避免别人照抄错的东西）**

- 初版 §7.4「两级顺序解、顺序不可颠倒」被 Q-10 实验推翻 → 改联合求解（D-09）。
- Q-10a 一度报"池化顺序解偏 0.03nm"，那是真值里 Z1=0 掩掉了病理；含低阶真值后它崩到 2.3nm。
  → `experiments.md` 里已标注更正。
- `fit` 一度给 `(2,0)` 的 12.5% 误差留了"留给 Q-04"，但那是我先看到数字才补的理由；
  记在这里以免下次反过来当成既有结论用。

---

## 更新记录

### 2026-10-11

- 设计：主线定 overlay（D-07）；不碰 C++、难度轴改统计推断（D-08）；Phase 3 升级为可校正性判定深度层。
- 规格纠错 8 处：`Z9=spherical`（Noll 下是 Z11）→ 已修；`79 field` → 实为 76/64；
  `regular` 采样说明与 §6.2.4 冲突；`make_fingerprint` 签名不可用；`measurements` 列数 18→19；
  `is_edge` 语义按实现修订；`score`/`report` 契约缺参数（A-04 算不出来）；19 点 mark 定义矛盾（U-10）。
- 实验：Q-10a/b/c 三个实验，七种估计器的效果汇总（`experiments/q10_estimator_comparison.py`）；
  定下"数字必须由仓库内脚本产出"的规矩。
- 实现：Phase 1 第 1–10 步全部完成（156 测试）。
- 新增文档：`docs/experiments.md`（实验记录）、`docs/design-faq.md`（设计问答）、本文件。
