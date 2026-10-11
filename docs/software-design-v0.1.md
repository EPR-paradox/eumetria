# eumetria 详细设计规格

> 项目：eumetria — 光刻量测场分析工具
> 文档状态：Draft v0.2
> 创建日期：2026-10-03
> 修订日期：2026-10-11
> 目标读者：实现者本人（+ 协作助手）
> v0.1 → v0.2 变更：主线定为 overlay（D-07）；不引入 C++，难度轴改为统计推断（D-08）；
> 新增 Q-08、A-09/A-10、Phase 2.5（FIM/CRB + D-optimal 采样设计）；U-04 定做、U-08 定不做；
> Phase 3 收窄、Phase 4 砍到最小；§3.5 量级表改为可执行断言
> 同日续：Phase 3 升级为「可校正性判定」深度层（不砍），新增 Q-09/Q-10、A-11
> 同日续：开工 Phase 1 第 1–3 步；据实现修正三处规格错误 —— §9.2 的「Z9=spherical」（Noll 下应为 Z11）、
> §7.1 的「79 field」（实为 76/64）、19 点 mark 密度定义矛盾（新增 U-10）
> 同日续：Phase 1 第 4–6 步（synth）；§7.3 regular 采样说明与 §6.2.4 冲突已修正、make_fingerprint 改收 SamplePoints、
> §3.4 补 TIS 组合的 σ/√2 结论、§6.2 列数 18→19 与 is_edge 语义按实现修订
> 同日续：据实测补 §3.1.1「CPE 阶数上限由最差场决定」（默认方案只支持一阶 CPE）+ 三条例律，
> 新增 U-11（边缘场 mark 是否自适应布置），Q-10 加边缘场三臂对照
> 同日续：**框架内决策反转 D-09** —— 实测推翻「两级顺序解」（顺序解偏 0.032–0.358nm，联合解 0.0000nm），
> 同时查出 Zernike Z1 与场内 Cartesian (0,0) 混叠（条件数 3.3e15 → 去掉场内常数项后 7.6）；
> 新增 docs/experiments.md（实验记录）与 docs/design-faq.md（设计问答），experiments/ 放可复现脚本
> 同日收尾：**Phase 1 第 1–10 步全部完成**（geometry → basis → synth → io → fit → metrics → cli → report），
> 验收 A-01/A-03/A-04/A-05 通过、A-02 判为「不适用（带噪）」；156 个测试全绿；下一步进 Phase 2 结论层
> 前置调研：SPIE 9778 / 9424（Zernike-Legendre overlay 建模）、SPIE 848200（YieldStar 量测原理）、
> MDPI 9(1)1（采样策略综述）、JM3 24(2025)024201（EUVlitho）、arXiv 2505.07576（SEM VAD 基准）

---

## 0. 约定与决策状态

### 0.1 已确认决策（D = Decision；含用户明确表达与用户批准的助手提案，反转项标注证据）

| 编号 | 决策 | 说明 |
|---|---|---|
| D-01 | 项目名 `eumetria`，仓库 `/home/j/myprojects/eumetria` | PyPI 可用、GitHub 同名 0 个 |
| D-02 | 主题：光刻**量测**侧的 wafer 场分析，三个指标 overlay / focus / CDU | 三者同源，可共用一套场分解框架 |
| D-03 | 不做 defect / inspection 方向 | WM-811K 类缺陷网格分类明确排除 |
| D-04 | 数据来源为自建正向合成器，不依赖任何真实 fab 数据 | 真实 overlay/focus/CDU 数据不公开 |
| D-05 | 实现语言 Python；本项目是**计算库 + CLI 工具**，不是 Web 服务 | 不引入 FastAPI / DB 服务 |
| D-06 | 定位：KLA/ASML 量测算法岗面试作品 + Python 数值计算练手 | 判断标准是"能否讲深"，不是功能多 |
| D-07 | **主线是 overlay**：overlay 做全链路，focus 用同一框架跑通即止，CDU 只做反演一条 | 深度集中在量测侧最核心、公式最厚的指标上；三线平权 = 三线都浅 |
| D-08 | **不引入 C++/native 层**；难度轴走统计推断（FIM/CRB + 最优采样设计） | 纯 Python 可达；工程能力由 Voxora/kairos 承担，本项目只负责算法深度 |
| D-09 | 两级拟合用**联合求解**（`X=[Z\|C]` 一次 lstsq，再按基块切开系数）；场内基不含常数项 | **决策反转**：初版 §7.4 的「先场内后全片、顺序不可颠倒」被实验推翻 —— 顺序解有偏 0.032–0.358nm，联合解 0.0000nm。证据：`experiments/q10_fit_order.py`、`docs/experiments.md` → Q-10a |

### 0.2 本版默认（P = Proposal，可改，不阻塞开工）

| 编号 | 默认选择 | 影响面 |
|---|---|---|
| P-01 | Python 3.10+，`src/` 布局，包名 `eumetria`，CLI 入口 `eumetria` | 与 Voxora 环境一致 |
| P-02 | 核心依赖仅 `numpy` + `scipy`；`matplotlib` 画图；`pydantic v2` 校验 schema；`pytest` 测试 | 依赖越少越好 |
| P-03 | Zernike **自己实现**（Noll 索引），不硬依赖 `zernpy`；`zernpy` 仅作对照后端 | 避免第三方 API 锁定，且自实现本身是加分项 |
| P-04 | 数据交换用 CSV + JSON 文件，**不引入数据库** | 分析工具，无并发写需求 |
| P-05 | 报告输出 PNG + 单文件 HTML，不做 Web UI | 面试演示够用 |
| P-06 | 单位：几何用 mm，误差量用 nm；Zernike 自变量为归一化坐标 ρ = r / R_wafer | 见 §9.1 |
| P-07 | 300mm wafer（R=150mm）、26×33mm field、3mm 边缘排除为默认场景 | 见 §3.6 |
| P-08 | 一切随机过程必须固定 seed，可 bit 级复现 | 见 §9.5 |

### 0.3 未决项（U = Unresolved）

| 编号 | 问题 | 安全默认 | 是否阻塞 Phase 1 |
|---|---|---|---|
| U-01 | 是否接入开源光刻仿真器（openuv / EUVlitho）生成 CD↔(focus,dose) 数据 | 先用解析 Bossung 模型；预留 adapter 接口（§11.2） | 否（阻塞 Phase 3） |
| U-02 | 是否把 WM-811K 作为第二输入源 | 不做。若做，只能作为独立的数据加载模块，不进主线 | 否 |
| U-03 | 是否发布到 PyPI | 不发。仓库公开即可 | 否 |
| U-04 | 是否实现自适应采样（adaptive sample plan） | **已定：做**。静态最优布局用 D-optimal 贪心；序贯加点按 FIM 信息增益逐步加，Phase 2.5 实现 | 否 |
| U-05 | 是否支持 200mm / 450mm 硅片与不同 field 尺寸 | 配置化即可，默认场景 300mm（§3.6） | 否 |
| U-06 | 是否做严格的**矢量** Zernike 分解（而非 x/y 分量独立分解） | 不做。x/y 独立分解是行业常规，矢量基作为 Phase 4 可选（§3.1、§15.3） | 否 |
| U-07 | field-level 用 Legendre 替代 Cartesian 幂基 | 默认 Cartesian（便于与扫描机系数对话）；Legendre 已在 §7.2 预留，切换成本低。**实测动因（Q-10c）**：幂基的偶次列（ξ²、η²）在采样点集上均值不为零，会把全片常量吸走 —— 含低阶真值时「池化顺序解」的 Z1 误差达 2.0nm；Legendre 的列按构造零均值，与常量的混叠显著更小 | 否 |
| U-08 | 是否实现带约束的校正优化（QP） | **已决：不做**。拟合不用 QP（§9.3）；§7.10 的校正优化明确划入"未做"，写进 §15.3 | 否 |
| U-09 | 鲁棒估计默认级别 | 默认 MAD 剔除（L1）；Huber/LAD 作为可配置升级（§7.9） | 否 |
| U-10 | mark 密度 19 点的定义（§7.1 原表述自相矛盾，无法同时满足「25-6」与「对称」） | 先用 4 / 8 点（已实现）；19 点要么改成 5×5 去掉中心 1 点与 4 个正交邻居（20 点），要么另定 | 否（4/8 点够跑完全部 Phase） |
| U-11 | 边缘场的 mark 是否自适应布置（往可用区内收），使每场都凑够 CPE 所需点数 | 先按固定 4/8 点布局；若 Q-10 显示边缘场是瓶颈，再把 mark 位置做成场景可配 | 否（阻塞 CPE 的最优解，不阻塞一阶 CPE） |

---

## 1. 目标与验收

### 1.1 项目定位

eumetria 是一个把**量测场数据分解为可解释成分**的工具。给一张硅片上的 overlay / focus / CD 测量点，它回答四个问题：

1. 这片硅片上的误差由哪些指纹成分构成？（场分解）
2. 哪些成分是扫描机可以校正的，哪些是残差？（可校正性）
3. 从测量指纹能反推什么物理原因？（诊断）
4. 换一套采样方案，结论还成立吗？（可辨识性）

**主线是 overlay。** 三个指标不平权：overlay 走完整链路（场分解 → 可校正性 → 采样设计 → 诊断），做到最深；focus 用同一套框架跑通即止，用来证明"框架可复用"这个架构论点；CDU 只做 CD→(dose, focus) 反演那一条 —— 逆问题 + 可辨识性是这一岗最能区分候选人的内容。理由：overlay 是量测侧最核心、公式最厚、与扫描机校正对话最直接的指标，且三线平权的结果是三线都浅。

**Non-goals（明确不做）：**
- 不做缺陷检测/分类（属 inspection，不是 metrology）
- 不接真实 fab 数据，不做设备通讯（SECS/GEM 等）
- 不做实时控制回路（APC 只建模其接口语义，不实现控制）
- 不做 Web 服务、不做多用户、不做权限
- 不追求物理仿真精度到"能替代商品软件"

### 1.2 核心闭环

```
scenario.yaml
    │
    ▼
[1] 正向生成 ──► measurements.csv  +  truth.json     ← 真值在此，唯一可信来源
    │            （模拟真实采样点，带 TMU/TIS 噪声）
    ▼
[2] 拟合     ──► fit_result.json
    │            （恢复的系数、残差、条件数、可辨识性警告）
    ▼
[3] 评分     ──► metrics.json
    │            （系数恢复误差 vs truth、残差 RMS nm）
    ▼
[4] 报告     ──► report/*.png + report.html
```

**关键设计约束：只有 [1] 知道真值。** [2][3][4] 只能看到 `measurements.csv`，不允许窥探 `truth.json`。这条隔离是项目可信度的基础 —— 面试时要能指着代码说明"评分器不知道真值，是无偏的"。

### 1.3 验收条件（可观测、可测试）

| 编号 | 条件 | 验证方式 |
|---|---|---|
| A-01 | 给定 `scenario.yaml` + seed，一条命令生成 `measurements.csv` 与 `truth.json` | `eumetria synth -c s.yaml -o out/` 退出码 0，两文件存在 |
| A-02 | 无噪声、采样点 ≥ 8× 系数个数时，全部系数恢复相对误差 < 1% | 单元测试断言 |
| A-03 | TMU = 0.1 nm 噪声下，前 6 阶主要项恢复相对误差 < 10% | 单元测试断言 |
| A-04 | 拟合输出的残差 RMS 与真值残差的差 < 0.01 nm | 断言 |
| A-05 | 同一 seed 两次运行结果 bit 级一致（或差 < 1e-12） | 回归测试 |
| A-06 | 采样点不足时**输出可辨识性警告并继续**，不崩溃、不静默给错结果 | 断言 warning 字段非空 |
| A-07 | 能产出采样预算曲线：点数 → 可稳定辨识的最高阶数 | 产物存在且单调 |
| A-08 | 全程不读取任何真实数据文件 | 代码审查 + 无网络调用 |
| A-09 | 系数置信区间覆盖率：蒙特卡洛 500 次下，95% CI 覆盖真值的频率 ∈ [93%, 97%] | 统计测试断言（覆盖率是可自动化的硬指标） |
| A-11 | 可校正性判定必须输出每项的显著性（p 值 / 置信区间）与所用 α，**禁止静默把项塞进 correctable**；换 α 后权衡表可复现 | 断言 + 产物检查（A-06 的加强版） |
| A-10 | 规模与数值稳定性：10 万点设计矩阵用 QR/SVD 求解（不用正规方程，否则条件数被平方），输出条件数，解与 `numpy.linalg.lstsq` 一致 | 数值断言 + 条件数随点数变化的记录 |

### 1.4 技术问题清单（项目要回答的，不是要实现的）

项目的价值在结论，不在代码量。以下问题各自需要一组对比实验，结论写进 README：

| 编号 | 问题 | 实验设计 |
|---|---|---|
| Q-01 | 正交基（Zernike）vs 非正交基（Cartesian）对系数稳定性的影响 | 同样拟合能力下，比较两组基的系数的 lot-to-lot 方差与条件数 |
| Q-02 | 全项拟合 vs AIC 精选：哪个泛化更好 | Z1–Z30 全项 vs AIC 逐步回归，比 CV RMSE |
| Q-03 | 采样点数与可辨识阶数的定量关系 | 点数梯度 × 阶数梯度，输出热图 |
| Q-04 | 噪声水平如何影响各阶系数的恢复精度 | TMU ∈ {0.02, 0.05, 0.1, 0.2, 0.5} nm |
| Q-05 | 采样布局（随机 vs 规则 vs field 分层）对可辨识性的影响 | 同点数不同布局对比 |
| Q-06 | CDU → (dose, focus) 反演的可辨识性边界 | 反演误差 vs 噪声/采样密度 |
| Q-07 | TIS 未校准会给 overlay 分解带来什么偏差 | 有/无 TIS 校准对比 |
| Q-08 | 采样布局离最优有多远：实测系数方差 vs Cramér-Rao 下界 | 对给定布局算 FIM = ΦᵀΦ/σ²，把 σ²·(ΦᵀΦ)⁻¹ 的理论下界与蒙特卡洛实测方差对比；再做 D-optimal 贪心布局作对照，报增益 |
| Q-09 | 可校正性判定的两类错误：漏校（真 CPE 被塞进 residual）vs 过校（噪声被当成 CPE） | 对每个候选项做系数显著性检验（t 检验 / F 检验），扫描显著性阈值 α，画漏校-过校权衡曲线；真值已知，两类错误都能直接计数 |
| Q-10 | 可校正集与 CPE 阶数怎么权衡：残差降多少、过拟合付多少 | 扫描扫描机支持项数（8/12/16/30）与 CPE 每场阶数，记录 residual RMS、拟合残差 vs CV 残差的过拟合交叉点、`cpe_ratio`（参数数/点数）；产出代价曲线。**另加三臂对照**：① 含边缘场、全局降到一阶 ② 排除 4 个边缘场、保二阶 ③ 逐场自适应降阶（预期最差，作反面对照）—— 各报 residual RMS + CV 残差 + 跨场系数变异，用来验证 §3.1.1「不许逐场降阶」这条纪律 |

### 1.5 本期不做

- 不做多硅片批量统计（Phase 4 可选：lot 级 APC 模拟）
- 不做 reticle / lens 级别的独立指纹分离
- 不做真实扫描机校正文件的导出格式（只做语义等价的 JSON）

---

## 2. 术语表

| 术语 | 含义 |
|---|---|
| **Overlay** | 当前层图案相对下层的套刻位置偏差。向量场 (dx, dy)，nm |
| **Focus** | 曝光时的焦面偏差。标量场，nm |
| **CD** | Critical Dimension，关键尺寸（线宽）。标量场，nm |
| **CDU** | CD 的全片均匀性指标（通常是 CD 的 3σ 或 range） |
| **inter-field / wafer-level** | 场间误差：把每个 field 视为一个点，描述整片形变 |
| **intra-field / field-level** | 场内误差：单个曝光场内部的畸变 |
| **correctable** | 扫描机能直接施加校正的项（线性 6 参数 + 设备支持的高阶） |
| **CPE** | Correction Per Exposure，逐场单独校正 |
| **residual** | 扣除模型后剩下的不可校正部分 |
| **TMU** | Total Measurement Uncertainty，量测总不确定度 |
| **TIS** | Tool Induced Shift，量测工具自身引入的系统偏移 |
| **Bossung 曲线** | CD 随 focus 变化的抛物线（不同 dose 为一族） |
| **FEM** | Focus-Exposure Matrix，焦面-剂量矩阵，用于标定 |
| **SWA** | Side Wall Angle，侧壁角，随 focus 单调变化，可反推 focus |
| **LER / LWR** | Line Edge / Line Width Roughness，随机粗糙度 |
| **Noll 索引** | Zernike 项的标准编号：Z1=piston, Z4=defocus, Z5/Z6=astigmatism, Z7/Z8=coma, Z9=spherical |
| **sample plan** | 采样方案：选哪些 field、每场哪些 mark 位置 |
| **edge exclusion** | 硅片边缘不可用环带（默认 3mm） |
| **field** | 一次曝光覆盖的区域（默认 26×33 mm） |
| **die** | 单个芯片；一个 field 可含多个 die |
| **mark** | 量测标记（光栅），置于 scribe line（切割道）上 |
| **DBO / IBO** | Diffraction / Image Based Overlay，两类 overlay 量测原理 |
| **APC** | Advanced Process Control，把模型系数反馈给扫描机的控制环路 |
| **k1 / k2 / k3** | 光刻 Rayleigh 缩放系数：R = k1·λ/NA、DOF = k2·λ/NA²，k3 为 DOF 近轴修正项。**属成像物理层，与 CPE 无直接关系**（详见 §3.2.1） |
| **DOF** | Depth of Focus，焦深 —— focus 预算的物理来源 |

---

## 3. 领域模型（物理与数学）

本章是实现的依据。所有公式在代码中必须有对应注释与单元测试。

### 3.1 Overlay 的误差分解

在硅片坐标 (x, y)（mm，原点为硅片圆心）与场内坐标 (xf, yf)（mm，原点为 field 中心）上，overlay 位移为：

```
dx(x,y) = dx_wafer(x,y)  +  dx_field(xf,yf)  +  ε_meas
dy(x,y) = dy_wafer(x,y)  +  dy_field(xf,yf)  +  ε_meas
```

**wafer-level 项（用 Zernike，定义域为单位圆）：**

```
ρ = r / R_wafer ,  θ = atan2(y, x)
dx_wafer = Σ_j  a_j^Zx · Z_j(ρ, θ)
dy_wafer = Σ_j  a_j^Zy · Z_j(ρ, θ)
```

其中 `Z_j` 为 Noll 第 j 项。x/y 分量各有一组独立系数（这是矢量场的常规做法；若要做矢量 Zernike 分解，见 U-06，本期不做）。

**field-level 项（用 Cartesian 多项式，定义域为归一化场内坐标）：**

```
ξ = xf / (W/2) ,  η = yf / (H/2)         ∈ [-1, 1]
dx_field = Σ_{p,q}  b_{pq}^x · ξ^p η^q      (1 ≤ p + q ≤ order_f，**不含常数项 (0,0)**)
dy_field = 同理
```

理由：field 是矩形，Legendre 在该域正交。本期默认用 Cartesian 幂基（行业惯用、便于与扫描机系数对话），Legendre 作为可选基（U-07）。

**场内基必须去掉常数项 (0,0)。** 它与 wafer-level 的 Z1（piston，恒等于 1）是**完全相同的列**，
两者物理上不可区分（"全片整体平移" ≡ "每个场都有同样大小的常量偏移"）；常量归 wafer-level 的
Tx = Z1（§3.1.1 的线性 6 参数表里就是这么写的）。实测：保留常数项时联合设计矩阵条件数
**3.3e15**（降秩），去掉后 **7.6**，而解不变。两个基单独看条件数都很正常（Zw 4.7 / Cf 6.6），
叠起来才炸 —— 所以"每个基单独没问题"推不出"合起来没问题"。证据：`docs/experiments.md` → Q-10a。

**线性 6 参数（wafer-level 的低阶特例，须能单独取出）：**

| 参数 | 数学形式 | 物理含义 |
|---|---|---|
| Tx | 常数（Z1） | 平移 X |
| Ty | 常数（Z1） | 平移 Y |
| Mag | ρcosθ（沿径向缩放） | 放大/倍率 |
| Rot | 旋转项 | 旋转 |
| NonOrth | 剪切项 | 非正交性 |

实现要求：从 wafer-level 系数出发，能换算成 industry-standard 的 `(Tx, Ty, Mag, Rot, NonOrth, MagY)` 形式。换算函数 `zernike_to_linear_terms()` 必须有单元测试（已知系数 → 换算 → 反算）。

#### 3.1.1 校正的两个正交维度：类型 × 粒度

行业里这两套术语经常被混着讲，其实是**两个独立的维度**。

**维度一：类型 —— 补偿什么误差**

即"线性 6 参数"。标准形式与命名（本文档全程只用这一套）：

    dx = Tx + Mag·x − Rot·y
    dy = Ty + Rot·x + Mag·y

| 参数 | 物理含义 |
|---|---|
| Tx / Ty | x / y 向平移 |
| Mag | 放大（倍率） |
| Rot | 旋转 |
| NonOrth | 非正交性（剪切） |
| MagY | 各向异性放大（与 Mag 成对使用） |

**术语规则（重要）：不要用 k1/k2/… 这类数字编号指代这些系数。** 两个原因：

1. **撞车** —— 光刻领域里 k1 已经是分辨率因子（R = k1·λ/NA）、k2 是焦深系数（DOF = k2·λ/NA²），见 §3.2.1
2. **不统一** —— 这套编号在不同工具和文献里指代不一致（例如把 k1/k2 当 x/y 平移、k3/k4 当旋转与放大），读者极易错位

本项目一律使用 `Tx / Ty / Mag / Rot / NonOrth / MagY`；描述粒度时用「线性 6 参数」这个说法。

**维度二：粒度 —— 在哪一层施加**

| 粒度 | 名称 | 参数规模 | 说明 |
|---|---|---|---|
| 全片 | wafer-level linear | 1 套/片（linear 10 terms） | 标准 APC run-to-run |
| 全片高阶 | HOPC / HOWA | 1 套/片，高阶 | 高阶场间校正 / 高阶对准 |
| 场内高阶 | i-HOPC | 场内 | 高阶场内校正 |
| **逐场** | **CPE = FxFc** | **6 参数/场** | 先进节点标配 |

**关键：CPE 的定义就是"每场一套 6 参数校正"。** ASML XT:1900i 的 GridMapper 原文（SPIE 813628）：

> "The CPE determines a 6-parameter correction for every exposed field giving the best interfield correction possible."

于是：

    全片线性校正   = 整片 1 套线性 6 参数
    CPE (FxFc)    = 每个 field 各自 1 套线性 6 参数

（FxFc = field-by-field correction，与 CPE 是同一概念的两个名字）

**所以线性参数与 CPE 既不是并列也不是包含，而是"类型"与"粒度"的正交组合 —— 任何粒度都可以装任何类型。**

**为什么需要 CPE：** 它补偿的是标准 APC 的 linear 10 terms 与高阶校正**都压不掉**的那部分系统性工艺指纹（wafer process signature）—— 各个 field 被不同局部条件影响，只有逐场处理才能吃掉。

**CPE 的代价（关键工程约束，容易忽略）：** SPIE 2258353 指出，从高阶 wafer-level 模型转向 CPE 会**大幅降低「参数数 / 测量点数」比率**（80 场 × 6 参数 = 480 个未知量），后果三条：

1. 必须先有**全片 field 采样**（full field sample plan），量测成本暴涨
2. 单点的随机误差（mark 保真度受随机效应影响）**不再被平均掉**，直接放大到校正系数上
3. 因此 CPE 控制环无法逐 lot 更新，只能按工艺稳定性低频更新

**对本项目的直接影响：** 这正是 §9.3 条件数/可辨识性问题的最坏情形 —— 参数数逼近点数时 LS 解失稳。所以 `fit_result.json` 的 `identifiability` 必须包含 **CPE 参数/点数比**（实现见 §7.5）。

**CPE 的阶数上限由「最差的那个场」决定**（实测数据，不是理论顾虑）。CPE 是逐场一套系数，
所以每个场各自的点数必须 ≥ 该场参数数 —— 平均点数够没有用：

| 场内 CPE 阶数 | 项数 | 未知量（x/y 两套） | 满点场（8 mark） | 边缘场（5 mark） |
|---|---|---|---|---|
| 一阶 `(0,0)(1,0)(0,1)` | 3 | 6 | 可辨识，设计矩阵条件数 1.50 | 可辨识，条件数 2.95 |
| 二阶 `+(2,0)(1,1)(0,2)` | 6 | 12 | 可辨识，条件数 6.61 | **降秩：rank 5 < 6，不可辨识** |

实测（seed 20261003，默认 25 场采样 / 8 mark）：21 个场满点、4 个边缘场只剩 5 mark
（field 6 / 15 / 60 / 69）。后果：**默认采样方案下二阶 CPE 直接不可辨识，CPE 上限是一阶。**

三条纪律：

1. **不许逐场自适应降阶。** CPE 的判据是跨场系数变异（Q-09/Q-10），各场阶数不同则系数
   不可比，整个可校正性判定作废；而且扫描机施加的是**同一套**模型，物理上不存在
   「第 6 场用二阶、第 60 场用一阶」这种校正。降阶还有个隐藏代价：降阶场的残差会被系统性
   放大，而残差正是 residual（不可校正）的估计 —— 等于把自己的降阶误读成「这片不可校正」。
2. **阶数上限 = min over 各场「可辨识的最高阶」**，全局统一档位。
3. **要在采样层解决，不是在拟合层**：给边缘场自适应布置 mark（往可用区内收，见 U-11），
   或者把边缘场排除出 CPE 拟合、并在 `identifiability.warnings` 里显式列出（A-06：不静默给错结果）。

**`np.linalg.lstsq` 在欠定系统上不报错** —— 5 点对 6 项时它静默返回最小范数解，
一个物理上无意义的「解」。所以 §7.4 的 `InsufficientPointsError` 必须由**显式的
`n_points ≥ n_params` 与秩检查**触发，不能依赖求解器抛异常。

**顺带一个数字（不是边缘场特有）**：8 点 scribe-line 布局的 mark 落在 ξ = ±0.769、η = ±0.818
处，只覆盖 field 归一化面积的 0.63×0.63 —— 满点场也在外圈 23%/18% 的区域外推。
所以二阶的条件数 6.61 是这个布局的**地板**，换采样布局才可能更低。

### 3.2 Focus

标量场，同一个 Zernike 基：

```
dz(x,y) = Σ_j  c_j · Z_j(ρ, θ)  +  ε_meas        [nm]
```

物理来源映射（用于 §3.7 诊断）：

| Zernike 模式 | 典型物理原因 |
|---|---|
| Z4 (defocus) | leveling 系统误差、硅片整体高度偏差 |
| 低阶径向（Z4, Z11, Z22…） | 硅片形貌（bow / warp / 夹持形变） |
| Z5 / Z6 (astigmatism) | 镜头像散 |
| Z7 / Z8 (coma) | 镜头彗差、对中误差 |
| Z9 (spherical) | 镜头球差 |
| 高阶（Z10 以上） | 场曲、高阶像差 |
| θ 相关的低阶（Z2/Z3、Z7/Z8） | 硅片旋转/对中 |

#### 3.2.1 与 k1 / k2 / k3 的关系（防术语混淆）

k1 / k2 / k3 是**成像物理层**的 Rayleigh 缩放系数，跟本书的场分解不是一回事：

    R   = k1 · λ / NA          分辨率（最小可分辨尺寸）
    DOF = k2 · λ / NA²         焦深
    k3                          DOF 的近轴修正系数（SPIE 4000.759, 2000）

- k1 物理极限 0.25（ASML 官方口径），老工艺 0.4–0.9；k2 典型 0.5–1.0
- k3 的来历：Rayleigh 的 DOF 式是近轴近似，严格解出光程差后焦深依赖**数值半孔径的平方**。用 k3 修正后，原先在 NA = 0.6 / 0.8 处的 DOF 被高估了 10% / 20%，并消除浸没式（n > 1）的歧义
- 它们回答"光刻机能印多小、焦深多深"，**既不是量测参数，也不是校正参数**

**与 CPE 的关系：没有直接关系，只有需求传导。**

| 层级 | 概念 | 回答的问题 |
|---|---|---|
| 成像物理 | k1 / k2 / k3 | 能印多小、焦深剩多少 |
| 量测分析 | Zernike 系数、残差 RMS | 这片硅片的误差是什么（本项目） |
| 校正控制 | 线性 6 参数 / 高阶 / CPE | 扫描机该施加什么、粒度多细 |

传导链：节点推进 → k1 逼近 0.25、DOF 按 λ/NA² 收紧 → on-product overlay（~4nm）与 CDU（~1nm）预算同时收紧 → 单靠 wafer-level 线性校正不够 → 必须升级到高阶 + 逐场校正（CPE）。

一句话：**k1 代表"要求"，CPE 是"手段"**，不在同一层。

**术语红线：** 本项目的模型系数命名是 `a_j`（Zernike）、`b_pq`（Cartesian）—— **不要用 k1/k2/k3**，那三个符号在这个领域已被分辨率与焦深占用。

### 3.3 CD 与 (dose, focus) 的耦合

本节是 CDU 分析的全部依据。

**Bossung 关系（解析形式，本期默认）：**

```
CD(d, f) = CD_0 + k_d · (d - d_0)  +  k_f · (f - f_0)^2  +  k_df · (d - d_0)(f - f_0)^2  +  ...
```

- 对 dose `d`：一阶（单调）
- 对 focus `f`：二阶（抛物线，极小值在最佳焦面）
- 不同 dose 给出同一形状、水平平移的一族抛物线 → Bossung 曲线族

实现约定：
- `BossungModel` 参数（CD_0, d_0, f_0, k_d, k_f）必须可配置，且提供从"物理直觉值"生成参数的默认工厂
- 单调区间约束：反演只在 CD 对 f 的**单侧支**上进行（左右支二义性必须显式处理，见 §10.3）

**CDU 的场分解（Phase 3 核心）：**

```
CD(x,y) = CD_mean
        + Σ_j p_j · Z_j(ρ,θ)        ← dose 场的一阶贡献
        + Bossung(field_focus(x,y)) ← focus 场的二阶贡献
        + LWR(x,y)                  ← 随机项，不可分解
```

反演目标：从 CD map 估计 dose 场与 focus 场的低阶系数。这是一个**病态逆向问题**（dose 一阶、focus 二阶，二者在低阶有耦合），必须报告可辨识性而非硬给答案。

### 3.4 量测误差模型

```
measured = true + TIS + N(0, σ_tmu)
```

| 项 | 默认值 | 说明 |
|---|---|---|
| σ_tmu (overlay) | 0.1 nm | 对标 YieldStar TMU < 0.1 nm |
| σ_tmu (focus) | 1.4 nm (3σ) → σ ≈ 0.47 nm | 对标 YieldStar focus 动态精度 1.4nm 3σ |
| σ_tmu (CD-SEM) | 0.3 nm | 由 `scenario` 配置 |
| TIS | 常数向量 (dx, dy)，默认幅值 0.3 nm，方向固定 | 校准方式见下 |

**TIS 校准模型：** 真实流程对同一 mark 做 0° 与 180° 两次测量：

```
m_0   = true + TIS + ε_0
m_180 = -true + TIS + ε_180          （180° 旋转使真实值反号）
TIS_hat = (m_0 + m_180) / 2
true_hat = (m_0 - m_180) / 2
```

合成器要**同时输出两次测量**（列 `dx_nm_0/dy_nm_0` 与 `dx_nm_180/dy_nm_180`），拟合器默认做该组合。未校准模式由配置开关控制（用于 Q-07 实验）。

**组合的副产品：噪声降到 σ/√2。** `(m_0 − m_180)/2 = true + (ε_0 − ε_180)/2`，方差
`(σ² + σ²)/4 = σ²/2`。所以 TIS 校准后的**等效 TMU 是 σ_tmu/√2**（默认 0.1nm → 0.0707nm）——
「每个 mark 测两次」的收益不只是消掉 TIS，还把量测精度提了 √2。写验收阈值（A-03）时要按
0.0707nm 算，不是 0.1nm。

### 3.5 真实量级参数表（合成器的标尺）

以下数值来自公开文献与厂商资料，作为 `synth` 默认值，**不得随意改动**（改则失去真实性锚点）：

| 量 | 真实量级 | 来源 |
|---|---|---|
| overlay TMU | < 0.1 nm | SPIE 848200 |
| focus 动态精度 | 1.4 nm (3σ) | SPIE 848200 |
| 高密度采样 | 1932 点 / 硅片，20 分钟内 | SPIE 848200 |
| 产线常规采样 | 20–25 field × 4–8 点 = 80–200 点/片 | 行业常规 |
| overlay 残差 RMS（全项 Z1–Z30） | 0.58 nm 拟合 / 1.21 nm CV | 中科院 Jos 25120038 表 1 |
| overlay 残差 RMS（AIC 精选 8 项） | 0.88 nm 拟合 / 0.91 nm CV | 同上 |
| overlay 残差 RMS（经验 4 项） | 3.32 nm 拟合 / 3.57 nm CV | 同上 |
| 先进节点 on-product overlay 要求 | ~4 nm | SPIE 881428 |
| 先进节点 on-product CDU 要求 | ~1 nm | SPIE 881428 |
| Zernike 基带来的改善 | 系数变异 −22%，overlay Mean+3σ −5%，良率 +0.1% | SPIE 9778 (HVM 实测) |
| Zernike 常用阶数 | 全项 Z1–Z30 / 精选 4–8 项 | Jos 25120038 |

**这些数字是 A-02/A-03 里误差阈值的依据。** 合成器的目标不是"看起来像"，而是"量级与统计特征与文献一致"。

**该表必须变成可执行断言**：`tests/property/test_magnitudes.py` 对生成器输出做区间校验（TMU 噪声 σ、focus 3σ、残差 RMS 随阶数的变化、Zernike 基带来的系数变异改善落在文献区间内），越界即失败。这样"锚定文献"是 CI 里跑的东西，而不是文档里的一句话 —— 被质疑合成器真实性时给断言，不给辩解。

### 3.6 几何模型

默认场景（300mm 硅片）：

| 参数 | 默认 | 说明 |
|---|---|---|
| 硅片直径 | 300 mm | R = 150 mm |
| edge exclusion | 3 mm | 可用区 R = 147 mm |
| field 尺寸 | 26 × 33 mm | 12×10 = 120 个候选网格位置 |
| field 布局 | 规则网格，中心对齐 | 数量口径见下，`build_layout` 的 `include_partial` 控制 |
| 每 field 的 mark 位置 | scribe line 上的固定位置集 | 默认 4 / 8 / 19 点三种密度 |
| 每 field 的 die 数 | 由 field 尺寸与 die 尺寸决定（默认不建模 die） | 本期不细分到 die |

参考量（已实测核算，300mm / 26×33 / 3mm 排除，12×10 = 120 个候选网格位置）：

| 口径 | 判据 | field 数 | rho_max |
|---|---|---|---|
| A 完全可用 | 四角都在可用圆内 | **64** | 0.847 |
| B 含 partial（**默认**） | field 中心在可用圆内 | **76** | 0.960 |
| C 有交叠 | 与可用区有任何几何交叠 | 104 | — |

面积等效值 82.4（不含排除）/ 79.1（含排除）**不是 field 数**，只作面积参考 —— 早期草稿曾误当作 field 数使用，已更正。200mm 硅片（R=100）面积等效约 37 个。

`wafer_layout` 必须输出每个 field 的中心坐标、四个角坐标、是否 partial、场内 mark 相对坐标。

### 3.7 诊断规则

把 §3.2 的模式-原因映射实现为可查询表，并支持"给出系数 → 输出候选原因列表（带置信度）"。这是把数学变成工艺语言的接口，也是面试最容易讲的部分。

规则形式：

```
{"rule_id": "FOCUS_FIELD_CURVATURE",
 "condition": {"basis": "zernike", "component": "focus", "dominant": ["Z11", "Z22"]},
 "cause": "field curvature / 高阶像差",
 "confidence": "medium",
 "notes": "与扫描机 field 校正能力比对后可判定是否可校正"}
```

---

## 4. 总体架构

### 4.1 数据流

```
             scenario.yaml
                  │
      ┌───────────┴───────────┐
      ▼                       ▼
  synth/  ────────────►  dataset/
  (正向生成)              ├── measurements.csv
      │                   ├── truth.json
      │                   └── layout.json
      │                        │
      │                        ▼
      │                   models/  ────────►  fit_result.json
      │                   (拟合)                  │
      │                                          ▼
      │                                     metrics.py  ──►  metrics.json
      │                                          │
      └──────────────►  report/  ◄───────────────┘
                        (图 + HTML)
```

**单向依赖，无环。** `models` 不得 import `synth`（否则评分器能拿到真值，破坏 §1.2 的隔离）。这条用测试强制：`test_no_synth_import_in_models.py`。

### 4.2 模块与职责

| 模块 | 职责 | 依赖 |
|---|---|---|
| `eumetria.geometry` | wafer / field / mark 几何布局、坐标变换、edge exclusion | 无（纯 numpy） |
| `eumetria.basis` | 正交基：Zernike(Noll) / Legendre / Cartesian；归一化坐标；导数 | geometry |
| `eumetria.synth` | 指纹生成、采样、噪声、TIS 模拟、Bossung 生成 | geometry, basis |
| `eumetria.fit` | 最小二乘拟合、AIC 选阶、交叉验证、条件数、可辨识性 | geometry, basis |
| `eumetria.robust` | 离群点检测（MAD）、IRLS/Huber、LAD；级联拟合 `fit_cascade` | fit |
| `eumetria.decompose` | correctable / CPE / residual 分离；线性 6 参数换算 | fit |
| `eumetria.correction` | 能力约束下的最优校正（QP，可选，Phase 4） | fit, decompose |
| `eumetria.inverse` | CDU → (dose, focus) 反演 | fit, basis |
| `eumetria.metrics` | 系数恢复误差、残差 RMS、预算曲线 | fit |
| `eumetria.report` | 图与 HTML | 全部（只读） |
| `eumetria.io` | CSV / JSON 读写与 schema 校验 | pydantic |
| `eumetria.cli` | 命令行 | 全部 |

### 4.3 明确不引入的东西

- 无数据库、无服务、无消息队列、无缓存层
- 无深度学习框架（本期不用 torch；理由见 §15.2）
- 无容器化要求（本地 `venv` 跑）
- 无网络调用（离线可运行）

---

## 5. 仓库结构

```
eumetria/
├── README.md                      # 项目说明 + Q-01..Q-07 的结论摘要（这是门面）
├── main.py                        # 零参数入口：python main.py 跑完整条流水线
├── pyproject.toml                 # 构建 + 依赖 + CLI entry point
├── .gitignore
├── docs/
│   ├── software-design-v0.1.md    # 本文档
│   ├── physics-notes.md           # 公式推导与文献摘录（实现依据）
│   ├── experiments.md             # Q-01..Q-10 的实验记录与结论（**数字的唯一出处**）
│   ├── design-faq.md              # 设计问答：为什么这么选、证据在哪、边界在哪
│   └── progress.md                # 进度：现在到哪 / 明天从哪接 / 什么还悬着
├── experiments/                   # 可复现的实验脚本（experiments.md 里的数字由它们产出）
├── scenarios/
│   ├── base-300mm.yaml            # 默认场景
│   ├── sparse-sampling.yaml       # 稀疏采样（Q-03）
│   ├── no-tis-calibration.yaml    # Q-07
│   └── cdu-inversion.yaml         # Q-06
├── src/eumetria/
│   ├── __init__.py
│   ├── geometry/{__init__,layout,coordinates}.py
│   ├── basis/{__init__,zernike,legendre,cartesian,registry}.py
│   ├── synth/{__init__,fingerprint,sampling,noise,tis,bossung}.py
│   ├── fit/{__init__,linear_fit,result,writers,cascade,model_select,crossval,diagnostics}.py
│   ├── robust/{__init__,outliers,irls,lad}.py
│   ├── decompose/{__init__,correctable,linear_terms}.py
│   ├── correction/{__init__,constrained_qp}.py
│   ├── inverse/{__init__,cdu.py}
│   ├── metrics/{__init__,recovery,budget}.py
│   ├── report/{__init__,plots,html}.py
│   ├── io/{__init__,schemas,readers,writers}.py
│   ├── pipeline.py                # 编排：synth→fit→score→report（三个前端共用这一份）
│   └── cli.py
├── tests/
│   ├── unit/                      # 基函数、几何、数值
│   ├── property/                  # 正交性、恢复性
│   ├── regression/                # 固定 seed 金标准
│   └── integration/               # 端到端 CLI
└── out/                           # 生成物（.gitignore）
```

---

## 6. 数据契约

所有文件格式在此固定。实现时用 pydantic 模型校验，读写失败必须报错而非静默补默认值。

内存表示是 `io.schemas.MeasurementSet`（numpy 数组 + 元数据的 dataclass），**不引入 pandas**：
188 行的表用不上 DataFrame，且 P-02 要求依赖最小。

### 6.1 scenario.yaml

```yaml
version: 1
seed: 20261003
wafer_id: W0001            # 补充字段：measurements.csv 需要它
lot_id: null               # 补充字段：允许为空（单片）
wafer:
  diameter_mm: 300
  edge_exclusion_mm: 3.0
field:
  width_mm: 26.0
  height_mm: 33.0
  marks_per_field: 8          # 4 | 8 | 19
sampling:
  mode: static                # static | adaptive
  fields: 25                  # 抽取多少个 field
  selection: regular          # regular | random | stratified
  tmu_nm: 0.1
tis:
  enabled: true
  amplitude_nm: 0.3
  angle_deg: 30.0
  calibrated: true            # true → 输出 0°/180° 两次测量并组合
fingerprint:
  basis: zernike              # zernike | cartesian
  max_order: 12               # Noll 索引上限（12 → Z1..Z12）
  terms:                      # 显式指定非零项与其幅值（nm）
    Z4:  2.5
    Z5: -1.8
    Z6:  1.1
    Z7:  0.9
    Z11: -2.2
  field_order: 2              # 场内 Cartesian 阶数
  field_terms:
    "1,0": 1.2
    "0,1": -0.8
    "2,0": 0.5
focus:
  enabled: false
  terms: {}
cdu:
  enabled: false
  bossung: {cd0_nm: 32.0, dose0_mj: 20.0, focus0_nm: 0.0, k_dose: 0.8, k_focus: -0.05}
  dose_terms: {}
  lwr_nm: 1.2
```

约束：`terms` 的键必须是合法 Noll 索引（`Z<n>`，n ≥ 1）；`max_order` 必须 ≥ 最大指定项；
`sampling.fields` ≤ 可用 field 数；`marks_per_field` ∈ {4, 8}（19 点见 U-10，尚未实现）。

`wafer_id` / `lot_id` 是相对初版的补充字段：CSV 的这两列按 §6.2 必须有值（lot 允许空），
而 scenario 是唯一的输入，所以它们必须来自 scenario。schema 用 `extra="forbid"`，
拼错的键会报错而不是被忽略。

### 6.2 measurements.csv

**这是 fit 的唯一输入。** 列定义：

| 列 | 类型 | 单位 | 说明 | 空值 |
|---|---|---|---|---|
| `wafer_id` | str | — | 硅片标识 | 不允许 |
| `lot_id` | str | — | lot 标识 | 允许（单片刻为 null） |
| `point_id` | int | — | 点序号 | 不允许 |
| `field_id` | int | — | field 序号 | 不允许 |
| `x_mm` | float | mm | 全场坐标（圆心为原点） | 不允许 |
| `y_mm` | float | mm | 同上 | 不允许 |
| `xf_mm` | float | mm | 场内坐标（field 中心为原点） | 不允许 |
| `yf_mm` | float | mm | 同上 | 不允许 |
| `rho` | float | — | 归一化径向坐标 = r / R_wafer | 不允许 |
| `theta_rad` | float | rad | 极角，[-π, π) | 不允许 |
| `dx_nm_0` | float | nm | 0° 测量 X 分量 | 无 TIS 时为 null |
| `dy_nm_0` | float | nm | 0° 测量 Y 分量 | 同上 |
| `dx_nm_180` | float | nm | 180° 测量 X 分量 | 同上 |
| `dy_nm_180` | float | nm | 180° 测量 Y 分量 | 同上 |
| `dx_nm` | float | nm | 最终使用的 X（已做 TIS 组合，或原始值） | 不允许 |
| `dy_nm` | float | nm | 最终使用的 Y | 不允许 |
| `dz_nm` | float | nm | focus 测量值（标量场场景） | 允许 |
| `cd_nm` | float | nm | CD 测量值（CDU 场景） | 允许 |
| `is_edge` | bool | — | **该点所属 field 是边缘场**（场内有 mark 落进排除区）| 不允许 |

**强约束：** `dx_nm` 与 `dy_nm` 之外的所有量都是"测量者可知"的；`xf_mm/yf_mm/rho/theta_rad` 由 `io.readers` 在载入时**重新计算并校验**（不允许信任文件里的值，防止真值泄漏）。

#### 6.2.1 rho 与 theta 是什么

Zernike 多项式**只定义在单位圆上**，所以必须先归一化：

    r     = sqrt(x_mm^2 + y_mm^2)     到硅片圆心的距离 [mm]
    rho   = r / R_wafer               归一化径向坐标, ∈ [0, 1]
    theta = atan2(y_mm, x_mm)         极角 [rad], ∈ [-π, π)

分母是**硅片半径 R_wafer（150mm）**，不是排除后的可用区半径（147mm）—— 归一化基准是物理硅片。

为什么必须归一化：

- Zernike 的正交性**只在 ρ ≤ 1 上成立**。把 mm 坐标直接代入径向多项式（含 ρⁿ 项）会得到 150ⁿ 量级的数，双精度下高阶项溢出，并且各项量纲严重失衡使设计矩阵病态
- ρ 表达的是"相对位置"：ρ=0 圆心、ρ=0.5 半径 75mm、ρ=1 硅片边缘。不同硅片尺寸（200/300/450mm）下同一 ρ 对应同一物理相位，模型可移植

验算（实测 point_id = 0）：

    (x, y) = (-107.000, -96.000) mm
    r      = sqrt(107² + 96²)  = 143.757 mm
    rho    = 143.757 / 150     = 0.9584                    ← 靠边缘的 partial field
    theta  = atan2(-96, -107)  = -2.4103 rad = -138.1°     （第三象限）

| rho | r (mm) | 位置 |
|---|---|---|
| 0.00 | 0.0 | 硅片圆心 |
| 0.50 | 75.0 | 半径过半 |
| 0.98 | 147.0 | 可用区边缘（R − edge_exclusion） |
| 1.00 | 150.0 | 硅片物理边缘 |

**设计域上限是 0.98，不是 1.0。** edge exclusion 内的 mark 没有工艺价值、不可测量，必须在采样阶段剔除。§9.2 要求拟合函数对 rho > 1 抛错是最外层防线；正常路径下不应出现 rho > 0.98。

#### 6.2.2 维度

三层维度，不要混：

| 层 | 形状 | 说明 |
|---|---|---|
| 文件（表） | (n_points, 19) | 一行 = 一个 field 内一个 mark 的一次测量 |
| 数值 | (n_points, 2) | 每点 2 通道 (dx, dy)；focus/CD 场景为 1 通道 |
| 拟合 | y: (2·n_points,)；X: (n_points, n_terms) | 见 §9.3 |

    n_points = Σ_被选 field (该场可用 mark 数)

**不是 field 数 × marks_per_field 的整齐乘积。** 边缘 partial field 的部分 mark 会落进 edge exclusion 被剔除（§6.2.4 第 6 步）。实测：25 field × 8 mark 理论上限 200，实际 **188**（4 个边缘场各只剩 5 个可用 mark）。

| 采样配置 | 理论上限 | 实测（含过滤） |
|---|---|---|
| 25 field × 8 mark（默认） | 200 | **188** |
| 76 field × 8 mark（全采样） | 608 | 边缘场更多，剔除比例更高 |

参数余量（默认）：188 点 vs 12 阶 wafer（2×12 = 24 项）+ field 2 阶（2×6 = 12 项）= 36 未知量 → 余量 5.2 倍。

#### 6.2.3 实测样例（seed 20261003，默认场景，由当前实现生成）

**文件里的浮点写全精度**（`repr(float)`，保证精确往返），下面表格里的 4 位小数只是显示 ——
文件要是只写 4 位，`readers` 的 §6.2.5 第 1 条重算校验（容差 1e-9）必然失败。

`field_id` 沿用 `build_layout` 的编号、不重编号，所以从 6 开始而不是 0。前 8 行：

    pt fid     x_mm     y_mm  xf_mm  yf_mm     rho    theta
     0   6 -127.000   -69.000 -10.00  13.50  0.9636  -2.6439
     1   6 -107.000   -96.000  10.00 -13.50  0.9584  -2.4103
     2   6 -107.000   -69.000  10.00  13.50  0.8488  -2.5688
     3   6 -107.000   -82.500  10.00   0.00  0.9007  -2.4848
     4   6 -117.000   -69.000   0.00  13.50  0.9055  -2.6087
     5   8   -75.000   -96.000 -10.00 -13.50  0.8122  -2.2340
     6   8   -75.000   -69.000 -10.00  13.50  0.6794  -2.3978
     7   8   -55.000   -96.000  10.00 -13.50  0.7376  -2.0911

field 6 的中心是 (-117, -82.5)：8 个候选 mark 有 3 个落进边缘排除区（ρ > 0.98），
只剩 5 个可用，所以前 5 行都属于它 —— **一行就是一个 mark 的一次测量，不是一场一行**。
（本版样例的行序由 mark 枚举顺序决定，与旧版本文档的样例顺序不同，ρ/θ 逐点一致。）

完整一行的 18 列（point_id = 0）：

```csv
wafer_id,lot_id,point_id,field_id,x_mm,y_mm,xf_mm,yf_mm,rho,theta_rad,dx_nm_0,dy_nm_0,dx_nm_180,dy_nm_180,dx_nm,dy_nm,dz_nm,cd_nm,is_edge
W0001,,0,0,-107.000,-96.000,10.000,-13.500,0.9584,-2.4103,3.3124,4.4941,-2.9888,-4.1188,3.1506,4.3064,,,false
```

注意第 11–14 列的过程：0° 与 180° 两次测量的**真值反号、TIS 同号**，所以组合后

    dx_nm = (3.3124 − (−2.9888)) / 2 = 3.1506
    TIS   = (3.3124 + (−2.9888)) / 2 = 0.1618     ← 被丢弃的那一半，用于校准验证

#### 6.2.4 正向生成步骤（顺序不可交换）

1. `build_layout` → 76 个 field（口径 B，带 `is_partial` 标记）
2. `sample_fields(n=25, mode='regular')` → **按空间均匀**抽取（5×5 分块各取一块；禁止按列表索引等距抽）
3. `field_marks(8)` → 每场 8 个候选 mark 的场内坐标（距边 3mm，落在 scribe line）
4. 全场坐标 = field 中心 + mark 场内偏移
5. `to_polar` → (rho, theta)
6. **mark 级过滤**：剔除 rho > (R − edge_exclusion)/R = 0.98 的 mark
7. `make_fingerprint` → 真实值

       dx_true(rho,theta) = Σ_j a_j^x·Z_j(rho,theta) + Σ_{p,q} b_{pq}^x·ξ^p·η^q
       dy_true(rho,theta) = 同构，x/y 各有独立系数集

8. `apply_tis(calibrated=True)` → 0° 与 180° 两组读数

       0°   : m_0   = +true + TIS + ε
       180° : m_180 = −true + TIS + ε

9. 组合 `dx_nm = (dx_nm_0 − dx_nm_180)/2`，TIS 被消掉
10. 写 `truth.json`：全部 a、b、TIS 向量、σ、seed

第 8 步产出的 `_0` / `_180` 四列保留在文件里 —— 它是量测流程的痕迹，也是 Q-07（TIS 未校准的偏差）的实验素材。拟合只使用组合后的 `dx_nm`。

#### 6.2.5 读取时的一致性校验（强制）

`io.readers` 载入 CSV 后必须：

1. 重算 rho/theta 与文件值比对（容差 1e-9），不一致即报错 —— 防真值泄漏
2. 断言 rho ≤ 1.0，超出抛 `InvalidMeasurementError`
3. 若四列 `_0`/`_180` 都在，断言 `dx_nm == (dx_nm_0 − dx_nm_180)/2`（容差 1e-9）
4. 断言同一 `field_id` 内 `(xf_mm, yf_mm)` 组合不重复

### 6.3 truth.json

```json
{
  "version": 1,
  "seed": 20261003,
  "wafer_id": "W0001",
  "basis": "zernike",
  "max_order": 12,
  "wafer_terms": {"Z4": 2.5, "Z5": -1.8, "Z6": 1.1, "Z7": 0.9, "Z11": -2.2},
  "wafer_terms_y": {"Z4": 1.9, "Z7": -1.3},
  "field_terms_x": {"1,0": 1.2, "0,1": -0.8, "2,0": 0.5},
  "field_terms_y": {"1,0": -0.6, "0,2": 0.4},
  "tis_vector_nm": [0.26, 0.15],
  "noise_sigma_nm": 0.1,
  "generated_at": "2026-10-03T16:40:00+08:00",
  "eumetria_version": "0.1.0"
}
```

**只有 `metrics` 模块可以读这个文件。** 代码层面通过参数隔离：`metrics.compare(fit_result, truth)` 里 truth 是显式入参，`fit` 的函数签名中不存在 truth。

### 6.4 fit_result.json

`fit` 只写它自己产出的键（version … identifiability、field_diagnostics）；
`linear_terms` / `correctable` / `cpe_terms` / `uncorrectable_residual_rms_nm` 由 Phase 3 的
`decompose` 补齐 —— **缺字段而不是写 0**（用 0 冒充"已知值"比缺字段危险）。
`field_diagnostics` 逐场给出秩与条件数，是 §3.1.1「CPE 阶数上限由最差场决定」的证据。

```json
{
  "version": 1,
  "basis": "zernike",
  "max_order": 12,
  "selected_terms": ["Z4", "Z5", "Z6", "Z7", "Z11"],
  "selection_method": "aic",              // aic | fixed | none
  "wafer_coeffs_x": {"Z4": 2.47, "Z5": -1.76, "Z6": 1.08, "Z7": 0.86, "Z11": -2.15},
  "wafer_coeffs_y": {"Z4": 1.85, "Z7": -1.27},
  "field_coeffs_x": {"1,0": 1.18, "0,1": -0.79, "2,0": 0.48},
  "field_coeffs_y": {"1,0": -0.58, "0,2": 0.39},
  "residual_rms_nm": 0.14,
  "residual_rms_nm_x": 0.11,
  "residual_rms_nm_y": 0.17,
  "n_points": 200,
  "n_params": 13,
  "condition_number": 42.7,
  "identifiability": {
    "max_identifiable_order": 12,
    "cpe_order_ceiling": 1,
    "warnings": ["以下场在 2 阶场内基下不可辨识（mark 数或秩不足）：6, 15, 60, 69 —— 逐场校正（CPE）不得使用该阶数（§3.1.1），也不许逐场降阶",
                 "CPE 阶数上限被最差的场限制为 1 阶（配置要的是 2 阶）：全局统一降档，或在采样层解决，不逐场降阶"],
    "rank_deficient": false
  },
  "field_diagnostics": {
    "6": {"n_marks": 5, "configured_order": 2, "rank": 5, "n_params": 6, "condition_number": null, "identifiable": false, "max_identifiable_order": 1},
    "8": {"n_marks": 8, "configured_order": 2, "rank": 6, "n_params": 6, "condition_number": 6.61, "identifiable": true,  "max_identifiable_order": 2}
  },
  "crossval": {"k_folds": 5, "group_by": "wafer", "rmse_nm": 0.13},
  "linear_terms": {"Tx_nm": 0.12, "Ty_nm": -0.08, "Mag_ppm": 1.4, "Rot_urad": 0.6, "NonOrth_urad": 0.2},
  "correctable": ["Tx_nm", "Ty_nm", "Mag_ppm", "Rot_urad", "NonOrth_urad", "Z4"],
  "cpe_terms": ["Z5", "Z6", "Z7"],
  "uncorrectable_residual_rms_nm": 0.12
}
```

### 6.5 metrics.json

```json
{
  "version": 1,
  "coefficient_recovery": {
    "Z4": {"truth": 2.5, "fit": 2.47, "abs_err": 0.03, "rel_err": 0.012},
    "Z5": {"truth": -1.8, "fit": -1.76, "abs_err": 0.04, "rel_err": 0.022}
  },
  "max_rel_err": 0.022,
  "residual_rms_nm": 0.14,
  "residual_agreement_nm": 0.004,
  "pass_A02": false,
  "pass_A03": true,
  "notes": []
}
```

### 6.6 版本与兼容

- 所有 JSON 带 `version`，跨版本不兼容时读取必须报错（`SchemaVersionError`）
- `measurements.csv` 增列允许（向后兼容），删列或改单位必须升 `version`
- 生成物（`out/`）不入 git

---

## 7. 模块接口契约

签名以类型标注形式给出。`NDArray[np.float64]` 简写为 `Arr`。所有角度弧度，所有坐标 mm，所有误差 nm（除非名称另有后缀）。

### 7.1 geometry

```python
@dataclass(frozen=True)
class Field:
    field_id: int
    cx_mm: float; cy_mm: float          # field 中心（全场坐标）
    width_mm: float; height_mm: float
    is_partial: bool
    is_edge: bool

@dataclass(frozen=True)
class Mark:
    mark_id: int
    xf_mm: float; yf_mm: float          # 场内坐标
    x_mm: float; y_mm: float            # 全场坐标

def build_layout(wafer_diameter_mm: float,
                 edge_exclusion_mm: float,
                 field_width_mm: float,
                 field_height_mm: float,
                 include_partial: bool = True) -> list[Field]:
    """规则网格铺满候选区（cs 网格 = ceil(2R/W) x ceil(2R/H)），中心对齐。
       include_partial=True（默认）→ 口径 B：field 中心在可用圆内即纳入，标 is_partial
       include_partial=False      → 口径 A：仅四角都在可用圆内的 field
       默认场景（300mm/26x33/3mm）应分别得到 76 / 64 个。
       任一 mark 落进 edge exclusion 则该 mark 标 is_edge=True。
       确定性，无随机。"""

def field_marks(field: Field, marks_per_field: int) -> list[Mark]:
    """返回场内 mark 位置。marks_per_field ∈ {4, 8, 19}：
       4  → 四角内侧固定偏移（默认距边 3mm）
       8  → 4 角 + 4 边中点
       19 → **待定（U-10）**：原表述「5×5 去掉中心 6 个」自相矛盾（25-6=19 既非 90° 也非
            180° 对称），当前实现显式拒绝该密度，不静默降级
    未知密度抛 ValueError，不得静默降级。"""

def to_polar(x_mm: Arr, y_mm: Arr, wafer_radius_mm: float) -> tuple[Arr, Arr]:
    """返回 (rho, theta)，rho = r / R_wafer ∈ [0,1]，theta ∈ [-π, π)。"""

def normalize_field(xf_mm: Arr, yf_mm: Arr,
                    field_width_mm: float, field_height_mm: float) -> tuple[Arr, Arr]:
    """返回 (xi, eta) ∈ [-1, 1]^2。"""
```

### 7.2 basis

```python
class Basis(Protocol):
    name: str
    def terms(self) -> list[str]: ...
    def values(self, coords: Arr, terms: list[str]) -> Arr:
        """返回形状 (n_points, n_terms) 的设计矩阵。"""

class ZernikeBasis:
    """Noll 索引。Z1=piston, Z2/Z3=tilt, Z4=defocus, Z5/Z6=astigmatism,
       Z7/Z8=coma, Z9/Z10=trefoil, Z11=spherical, ...（**Noll；j≥9 起与 ANSI/Fringe 编号错位**）
       归一化：单位方差 N_n^m = sqrt(2(n+1)/(1+δ_m0))，即 (1/π)∬Z_jZ_k ρdρdθ = δ_jk。
       max_order 是最大 Noll 索引 j（§9.4），不是最大径向阶数 n。"""
    def __init__(self, max_order: int): ...
    def values(self, rho: Arr, theta: Arr, terms: list[str]) -> Arr: ...
    def radial(self, n: int, m: int, rho: Arr) -> Arr: ...   # R_n^m
    def noll_to_nm(self, j: int) -> tuple[int, int]: ...      # 索引换算

class CartesianBasis:
    """幂基 ξ^p η^q，(p+q) ≤ order。terms 形如 "p,q"。"""
    def __init__(self, order: int): ...
    def values(self, xi: Arr, eta: Arr, terms: list[str]) -> Arr: ...

class LegendreBasis:
    """可选后端，定义域 [-1,1]^2。"""
```

**必须通过的性质测试（§13.2）：** 正交性（非对角内积 ≈ 0，相对容差 1e-10）、在单位圆外的定义域处理、与 `zernpy` 的数值一致性（对照后端，容差 1e-8）。

### 7.3 synth

```python
def sample_points(fields: list[Field], marks_per_field: int,
                  wafer_radius_mm: float, edge_exclusion_mm: float) -> SamplePoints:
    """把抽中的 field 展开成 mark 点表，并按 edge exclusion 过滤（ρ ≤ 0.98）。
       返回的 SamplePoints 是 measurements.csv 的唯一来源：field_id、点号、全场/场内
       坐标、ρ、θ 一次算清 —— 否则采样逻辑会有两份实现，迟早漂移。"""

def make_fingerprint(points: SamplePoints,
                     basis: Basis,
                     wafer_terms_x: dict[str, float],
                     wafer_terms_y: dict[str, float],
                     field_basis: Basis,
                     field_terms_x: dict[str, float],
                     field_terms_y: dict[str, float]) -> Arr:
    """返回未加噪的真实 (dx, dy)，形状 (n_points, 2)，单位 nm。"""

def sample_fields(layout: list[Field], n_fields: int, mode: str,
                  rng: np.random.Generator, blocks: int = 5) -> list[Field]:
    """mode ∈ {regular, random, stratified}，抽 n_fields 个 field。
       regular    = 在 field 中心的包围盒上做 5×5 空间分块，每块取距块心最近者（确定性）
       stratified = 同样分块，每块在块内随机取一个（由 rng 决定）
       random     = 全局无放回随机
       **不是「按列表索引等距抽取」** —— layout 是网格遍历序（先 y 后 x），等距索引
       会把抽到的点全赶到某一条带上（实测 25 点等距抽出的 ρ 全部 > 0.75），
       完全失去空间代表性。见 §6.2.4 第 2 步。"""

def add_noise(true: Arr, sigma_nm: float, rng: np.random.Generator) -> Arr: ...

def apply_tis(true: Arr, tis_vector_nm: tuple[float, float],
              calibrated: bool, sigma_nm: float,
              rng: np.random.Generator) -> dict[str, Arr]:
    """calibrated=True 时返回 0°/180° 两组测量；False 时只返回带 TIS 的单组。"""

def bossung(dose_mj: Arr, focus_nm: Arr, params: BossungParams) -> Arr: ...
```

### 7.4 fit

```python
@dataclass(frozen=True)
class FitConfig:
    basis: str = "zernike"
    max_order: int = 12
    selection: str = "none"         # none | fixed；aic 属 Phase 2（Q-02）
    fixed_terms: list[str] | None = None
    field_order: int = 2

def field_terms_for(order: int) -> list[str]:
    """场内项集：order 阶幂基，**去掉常数项 (0,0)**（D-09）。"""

def build_design(measurements: MeasurementSet,
                 config: FitConfig) -> tuple[Arr, list[str], list[str]]:
    """返回 (X, wafer_terms, field_terms)；X = [Z(ρ,θ) | C(ξ,η)] 按列拼接。"""

def fit_field(measurements: MeasurementSet,
              config: FitConfig) -> FitResult:
    """两级模型，**一次联合求解**（D-09）：
       X = [Z(ρ,θ) | C(ξ,η)]，np.linalg.lstsq(X, y, rcond=None)，再按基块切开系数
       分别报告 wafer_coeffs_* / field_coeffs_*。

       **签名里没有 layout**：`MeasurementSet` 自带场中心、场尺寸与场内坐标，且 (ρ,θ)
       在读取时已重算校验过；再传一份 layout 只会多一个能和它矛盾的信息源。

       **不要用「先场内、再用残差拟合全片」的一次扫描顺序解**（初版 §7.4 的规定，已被实验
       推翻）：两块基在同一批采样点上不正交（列中心化后最大 |cos| = 0.173），顺序解有偏，
       且偏差方向取决于先后。同一份无噪声数据（seed 20261003，188 点）实测：

           联合求解 / backfitting 迭代到收敛    wafer 非零项误差 0.0000 nm，残差 RMS 0.0000
           先场内后全片（一次扫描）             0.0316 nm / 0.1270 nm
           先全片后场内（一次扫描）             0.3579 nm / 0.2470 nm
           逐场各自拟合二阶后再拟合全片          2.6985 nm / 3.8153 nm   ← 灾难，见 Q-10 反面对照

       分层是可解释性的需求（wafer 系数给 APC、场内系数给 i-HOPC/CPE），不是估计器的限制：
       切开系数向量就能分层。若坚持交替结构，必须迭代到收敛（线性模型下与联合解等价）。

       记录条件数；点数 < 参数数或降秩时抛 InsufficientPointsError（显式检查，见 §3.1.1）。
       证据与复现：experiments/q10_fit_order.py、docs/experiments.md → Q-10a。"""

def select_terms(X: Arr, y: Arr, max_order: int, method: str) -> list[str]:
    """method='aic' 时做前向逐步回归，AIC 准则；'none' 返回全部项。"""

def cross_validate(measurements, layout, config) -> CrossValResult:
    """按 group_by 分组做 K 折。group_by='wafer' 时同片点不跨折（防泄漏）。"""
```

### 7.5 decompose

```python
def zernike_to_linear_terms(coeffs: dict[str, float],
                            wafer_radius_mm: float) -> LinearTerms:
    """把 wafer-level 系数换算为 (Tx_nm, Ty_nm, Mag_ppm, Rot_urad, NonOrth_urad)。
    必须有往返测试：linear → zernike → linear 误差 < 1e-9。"""

def split_correctable(coeffs: dict[str, float],
                      scanner_spec: ScannerSpec,
                      cpe_terms: list[str]) -> CorrectableSplit:
    """按两个正交维度切分（§3.1.1）：
         类型维度：线性 6 参数（Tx/Ty/Mag/Rot/NonOrth/MagY）→ 高阶 Zernike
         粒度维度：wafer-level → field-level(i-HOPC) → CPE(FxFc，逐场)
       返回 (correctable, cpe, residual)：
         correctable = 扫描机能全片施加的（线性 + spec.max_order 内的高阶）
         cpe         = 需逐场施加的（每场一套 6 参数）
         residual    = 两者都压不掉的
       ScannerSpec 还必须提供 cpe_enabled 与每场的参数上限。
       必须同时输出 cpe_ratio = (n_fields x 6) / n_points，超过阈值（默认 0.3）
       时写入 identifiability.warnings —— 这是 CPE 最现实的工程约束（§3.1.1）。"""
```

### 7.6 inverse

```python
def invert_cdu(cd_nm: Arr, coords: dict[str, Arr],
               bossung: BossungParams,
               dose_basis: Basis, focus_basis: Basis,
               dose_terms: list[str], focus_terms: list[str],
               branch: str = "negative") -> InverseResult:
    """从 CD 场反演 dose/focus 场系数。
       branch: 'negative' | 'positive' | 'auto'（自动选支需显式报告选择了哪一支）。
       输出必须含 covariance 与 identifiability 字段；不可辨识时给警告而非结果。"""
```

### 7.7 metrics

```python
def coefficient_recovery(fit: dict[str, float],
                         truth: dict[str, float]) -> dict[str, dict[str, float]]: ...

def residual_rms(predicted: Arr, measured: Arr) -> float: ...

def sampling_budget_curve(scenario_path: str,
                          point_counts: list[int],
                          max_order_scan: list[int]) -> pd.DataFrame:
    """输出列：n_points, order, coeff_rel_err_median, coeff_rel_err_p95,
       residual_rms_nm, condition_number。用于 Q-03。"""
```

### 7.8 report

```python
def plot_wafer_map(values: Arr, x_mm: Arr, y_mm: Arr, kind: str, **kw) -> Figure:
    """kind ∈ {scalar, quiver, categorical}。复用 matplotlib；
       如需交互式可 import 外部包（wafermap/wfmap），但不作为硬依赖。"""

def render_html(figures: list[Path], metrics: dict, out_path: Path) -> None: ...
```

### 7.9 分级拟合与鲁棒估计

#### （1）低阶/高阶分开讨论：取决于基是否正交

**这是本节最重要的一条结论。**

- **正交基（Zernike / Legendre）下：一次全项拟合 ≡ 级联拟合。**
  先低阶拟合、再对残差拟合高阶，与一次性全项拟合结果几乎一致 —— 这正是正交性的直接推论，也是 SPIE 9424 那篇论文的结论（它明确指出该性质让 dual-pass / cascading 这类技巧变得**不必要**）。
- **非正交基（Cartesian 幂基）下：级联（cascading / dual-pass）是行业常规做法。**
  一次性拟合全阶会因项间共线导致系数量纲不稳、lot 间抖动；级联让低阶系数先稳定下来，高阶只吃残差。

所以"要不要分开"不是一个风格问题，而是基选择的结果。项目要做的是**把这件事量化出来**：

```python
@dataclass(frozen=True)
class CascadeConfig:
    stages: list[list[str]]      # 每一级拟合的项，如 [["Z2","Z3","Z4"], ["Z5","Z6","Z7","Z8"]]
    propagate_residual: bool = True

def fit_cascade(measurements, layout, basis, config: CascadeConfig) -> CascadeResult:
    """逐级拟合：每级只对上一级的残差拟合。
       输出：各级系数、最终残差、以及与一次性全项拟合的系数差向量。
       CascadeResult 必须含 single_pass_delta（用于 Q-01 的等价性证明）。"""
```

必须实现的等价性实验（Q-01 的核心证据）：Zernike 基下 `fit_cascade` 与 `fit_field` 的系数差 < 1e-8；再把同一个对比换成 Cartesian 基，差异显著。**这一组对比就是"为什么用 Zernike"的数值证据**，比任何口头解释都硬。

**别把"分级拟合"和"语义分层"混为一谈：**

| | 目的 | 正交基下是否必要 |
|---|---|---|
| 分级拟合（cascading） | 数值手段，保证系数稳定 | **不必要**（等价于一次拟合） |
| 语义分层（§7.5 可校正/CPE/残差） | 物理与控制含义，给扫描机下参数 | **始终必要** |

#### （2）鲁棒估计：真实数据里一定有离群点

真实 overlay / focus 测量必然有坏点：mark 被前道工艺破坏、量测失效、局部缺陷。LS 对离群点极敏感（一个 10σ 点足以带偏整片系数），工业流程里剔除离群点是标准步骤。

三级策略（`eumetria.robust` 模块）：

| 级别 | 方法 | 默认 |
|---|---|---|
| L1 | MAD 判据剔除：`|r_i − med(r)| > k·1.4826·MAD`，k 默认 4 | **默认开启** |
| L2 | IRLS + Huber 权重（δ = 1.345σ） | 配置开启 |
| L3 | LAD（L1 损失）→ `scipy.optimize.linprog` | 离群点比例 > 10% 时 |

```python
def detect_outliers(residuals: Arr, k: float = 4.0) -> np.ndarray:
    """返回 bool mask（True = 判为离群点）。基于中位数与 MAD，不用均值方差。"""

def fit_robust(X: Arr, y: Arr, method: str = "mad", max_iter: int = 50) -> tuple[Arr, dict]:
    """返回 (coeffs, diagnostics)。diagnostics 必须含：
       n_outliers, outlier_ratio, iterations, converged, method。"""
```

规则：**默认先跑一次 LS，再跑鲁棒版，两者的系数差写进 `fit_result.json`。差异大本身就是诊断信号**（说明数据存在问题，而不是说明鲁棒解"更好"）。

### 7.10 校正优化（可选，Phase 4）

**先明确：这不是拟合。** 拟合回答"这片硅片的误差是什么"；校正回答"扫描机该施加什么参数"。区别就在**约束**：

```python
def optimize_correction(coeffs: dict[str, float],
                        scanner_spec: ScannerSpec,
                        weight: Arr | None = None) -> CorrectionResult:
    """求解能力约束下的最优校正参数：

           min   || y − X·c ||²
           s.t.  lo_j ≤ c_j ≤ hi_j        ← 盒式约束：每条校正通道的物理上下限

       这是一个 QP。实现用 scipy.optimize.minimize(method='SLSQP')（盒式约束用 bounds 表达，
       不必引入 cvxpy）。返回 chosen_c、unconstrained_c、active_constraints、
       predicted_residual_rms_nm。"""
```

为什么值得做：**未约束的最优解经常超出设备能力**（比如算出需要 5nm 的高阶校正，设备只能给 2nm）。此时真正的问题是"在能力边界内如何最优"，这才是 QP 的正确位置。能讲清"拟合 ≠ 校正"这件事，比会用 QP 有价值得多。

`ScannerSpec` 由 `scanner_spec.yaml` 配置（虚构但量级合理：线性 6 参数范围宽、高阶项收紧），**必须在 README 声明这些是假设值**，不是从真实设备文档来的。

---

## 8. CLI 契约

```
eumetria synth   -c scenarios/base-300mm.yaml -o out/base
eumetria fit     -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
eumetria score   --fit out/base/fit_result.json --truth out/base/truth.json \
                 -i out/base/measurements.csv -c scenarios/base-300mm.yaml -o out/base
eumetria report  -i out/base -o out/base/report
eumetria run     -c scenarios/base-300mm.yaml -o out/base      # 上述四步串联
eumetria sweep   -c scenarios/base-300mm.yaml --sweep points=20,40,80,160,320
                 --sweep order=4,6,8,10,12 -o out/sweep    # Q-03 采样预算
eumetria sample-plan -c scenarios/base-300mm.yaml --terms 16 --points 25
                 --optimize d-optimal -o out/plan            # Q-08 采样设计（附 CRLB 报告）
eumetria sample-plan -c scenarios/base-300mm.yaml --terms 16 --sequential --budget 40
                 -o out/plan_seq                             # U-04 序贯加点（按 FIM 信息增益）
```

前端分三层，**编排只有一份**（都在 `pipeline.py`），互不复制串联逻辑：

| 入口 | 用途 | 做什么 |
|---|---|---|
| `python main.py` | 零参数跑全流程；PyCharm 绿三角 / 新克隆后先试这个 | 默认值 + 人类可读输出 |
| `eumetria <子命令>` | 完整 CLI，可单独跑四步中任一步 | 参数解析 + 退出码映射 |
| `pipeline.run_all(...)` | notebook / 脚本 / 集成测试 | 直接调 |

约定：
- 退出码：0 成功；1 输入/schema 错误；2 计算错误（如点数不足）；3 内部错误
- **`score` 相对初版增加了 `-i/--input` 与 `-c/--config`**：A-04（拟合残差与真值残差之差
  < 0.01nm）必须重建真值场才能判定，而它需要点坐标（measurements）与读 CSV 时的硅片半径
  （scenario，§6.2.5 强制）。原契约只给 fit/truth 两个路径，A-04 无从计算。
- `run` 的四步共用同一个 `-o` 目录，`-i/--fit/--truth` 由它按固定文件名推出，不必重复输入
- `--seed` 可覆盖 scenario 中的 seed
- 所有命令支持 `--json` 输出机器可读摘要到 stdout
- 日志走 stderr，不污染 stdout

---

## 9. 数值与实现约定

### 9.1 单位

| 类别 | 单位 | 强制方式 |
|---|---|---|
| 坐标 | mm | 变量名后缀 `_mm` |
| 误差量 | nm | 后缀 `_nm` |
| 角度 | 弧度 | 后缀 `_rad`；角度值用 `_deg` |
| 线性项 | Tx/Ty 用 nm；Mag 用 ppm；Rot/NonOrth 用 μrad | 见 `LinearTerms` |

禁止在函数内部隐式换单位；所有换单位必须走 `geometry.coordinates` 里的显式函数。

### 9.2 Zernike 约定

- 索引：**Noll**（Z1 = piston, Z4 = defocus, **Z11 = spherical**；Z9 = trefoil）。
  **不许把 ANSI / Fringe 的「球差 = 9」混进来** —— 两套编号从 j=9 起错位，混用会让高阶像差全部对错项。
  执行版本：`basis/zernike.py` 模块 docstring 的映射表 + `tests/test_basis.py::test_noll_convention_traps_the_ansi_confusion`
- 归一化：文档与代码注释必须写明所用归一化（本版采用单位方差约定，`∫ Z_j Z_k dA / π = δ_jk`）
- 定义域：`rho ∈ [0,1]`。`rho > 1` 的点（edge exclusion 外）必须在采样阶段就被排除，拟合函数对 `rho > 1` 抛 `ValueError`（不允许静默截断）。
  这不是纯理论约束：口径 B 的 partial field 会产生**片外 mark**（实测 `rho_max = 1.061`，即 r ≈ 159mm > 150mm），
  所以「按 edge exclusion 过滤 mark」是 `synth.sampling` 的必需步骤，不做就会在这里炸
- 实现方式：径向多项式递推 + `noll_to_nm` 映射表；不使用硬编码系数表（避免抄错）

### 9.3 线性代数与求解器选择

**默认：无约束最小二乘。** 这是 overlay/field 建模的行业标准，也是本项目的基线。

- 统一用 `np.linalg.lstsq(X, y, rcond=None)`
- 每次都计算并记录 `condition_number = np.linalg.cond(X)`
- 条件数 > 1e6 时写入 `identifiability.warnings`
- 项间相关系数矩阵也要算，> 0.9 的项对写进警告

**求解器决策表（按需求升级，不要预先过度设计）：**

| 需求 | 方法 | 实现 | 位置 |
|---|---|---|---|
| 无约束拟合（基线，绝大多数场景） | 普通最小二乘 | `np.linalg.lstsq` | §7.4 |
| 离群点鲁棒（mark 被工艺破坏、量测失效） | MAD 剔除 + LS；或 IRLS + Huber | 手写 IRLS（约 15 行） | §7.9 |
| 重尾噪声下的鲁棒拟合 | LAD（L1 损失）→ 线性规划 | `scipy.optimize.linprog` | §7.9（可选） |
| 项选择（抗过拟合） | AIC 前向逐步（默认）/ L1(LASSO) | 手写逐步；L1 用坐标下降 | §7.4、Q-02 |
| **扫描机能力约束下的最优校正** | **QP（盒式约束）** | `scipy.optimize.minimize(method='SLSQP')` | §7.10（可选） |

**关键区分：拟合（fitting）≠ 校正（correction）。**

- **拟合**：从测量数据估计系数。模型对系数线性、无物理约束 → 无约束 LS 就是正确答案。
- **校正**：在扫描机**能施加的范围**内（每条校正系数有上下限、高阶校正能力有限）选一组值最小化残差 → 这是带盒式约束的优化问题，**QP 的正确位置在这里**，不在拟合。

**不引入 QP 做拟合的理由（重要）：**
1. overlay 模型对系数线性且无约束，用 QP 解无约束 LS 是纯粹的复杂度浪费（引入 `cvxpy`/`SLSQP` 依赖，换来 0 精度提升）
2. 面试中会立刻被问"你的约束条件是什么？"—— 答不上来就暴露了"为了显得高级而堆工具"
3. 真正的数值难点在别处：病态矩阵（条件数）、离群点、采样不足 —— 这些用 SVD/IRLS/信息准则解决，不是 QP

### 9.4 阶数上限

默认 `max_order = 12`（Z1..Z12），覆盖 defocus / astigmatism / coma / spherical 及一阶更高阶项。理由：field-level 默认 2 阶、wafer-level 更高阶在采样点少时不可辨识（Q-03 会量化）。上限由 scenario 控制，不硬编码。

### 9.5 可复现

- 所有随机数走 `np.random.default_rng(seed)`，seed 来自 scenario
- 生成器函数**不接受**全局随机状态
- 回归测试固定 seed 后比对金标准数值（容差 1e-12）
- 禁止依赖 dict/set 迭代顺序影响数值结果（需要遍历时显式 `sorted()`）

### 9.6 浮点容差

| 场景 | 容差 |
|---|---|
| 正交性检验 | 相对 1e-10 |
| 系数往返换算 | 1e-9 |
| 与 zernpy 对照 | 1e-8 |
| 回归金标准 | 1e-12 |

---

## 10. 计算流水线与失败路径

### 10.1 阶段与失败模式

| 阶段 | 失败模式 | 处理 |
|---|---|---|
| layout 构建 | field 尺寸 > 硅片直径 | `ValueError`（配置错误，退出码 1） |
| sampling | 请求 field 数 > 可用数 | `ValueError`，报出可用上限 |
| 指纹生成 | term 名非法 / max_order 不足 | `SchemaValidationError`（退出码 1） |
| 噪声/TIS | sigma < 0 | `ValueError` |
| fit | 点数 < 参数数，或设计矩阵降秩 | `InsufficientPointsError`（退出码 2），信息含"需要 ≥ N 点"。**必须显式检查，不能依赖 lstsq —— 欠定系统上 lstsq 不报错，会静默返回最小范数解** |
| fit | 设计矩阵秩亏 | 不抛异常；置 `rank_deficient=true` + 警告，返回最小范数解 |
| fit | 条件数过大 | 同上，给警告 |
| crossval | 折内点数不足 | 降折数并警告（最低 2 折），2 折仍不足则跳过 CV 置 null |
| inverse | 反演病态 / 选支不确定 | 给警告与候选解，不硬给唯一答案 |
| report | 无数据 | 跳过该图，记录到 notes，不失败 |
| io | version 不匹配 | `SchemaVersionError`（退出码 1） |

### 10.2 事务性

本项目无外部副作用，无事务需求。所有输出写到**临时目录后原子 rename**，避免半成品文件被后续步骤读取：

```
out/base/.tmp-measurements.csv  →  out/base/measurements.csv
```

### 10.3 Bossung 反演的二义性（重点）

CD 对 focus 是抛物线 → **同一个 CD 值对应两个 focus**（左右支）。设计上必须：

1. `branch` 参数显式指定，或 `auto` 时依据"同一 field 内相邻 mark 的 CD 单调性"判支
2. 无法判支时输出**两个解**并标记 `ambiguous=true`
3. 绝不静默选一支 —— 这是本项目最容易出错、也最容易在面试中被问到的点

### 10.4 幂等

同 scenario + 同 seed → 输出字节级一致。`run` 命令可重复执行，覆盖旧输出（不追加）。

---

## 11. 外部依赖边界

### 11.1 硬依赖

| 包 | 用途 | 封装位置 | 锁定策略 |
|---|---|---|---|
| `numpy` | 数值基础 | 全局 | `>=1.26,<3` |
| `scipy` | 统计/插值（如需） | `fit`, `inverse` | `>=1.11,<2` |
| `pandas` | CSV 读写与表操作 | `io` | `>=2.0,<3` |
| `pydantic` | schema 校验 | `io.schemas` | `>=2.0,<3` |
| `matplotlib` | 出图 | `report` | `>=3.8,<4` |
| `pyyaml` | 配置 | `io` | `>=6.0,<7` |
| `pytest` | 测试（dev） | `tests` | 任意 |

**`zernpy` 不是硬依赖** —— 只用于对照测试（`tests/property/test_zernike_vs_zernpy.py`，标记 `@pytest.mark.optional`，未安装则 skip）。

### 11.2 可选外部（预留接口，本期不实现）

```python
class LithographySimBackend(Protocol):
    """未来接入开源光刻仿真器，用于生成物理真实的 CD↔(focus,dose) 关系。"""
    def bossung_from_simulation(self, dose_range, focus_range, config) -> Arr: ...
```

候选实现（U-01，接入前必须先自己验证其输出合理性）：
- `EUVlitho`（takahashi-edalab，MIT，有 JM3 论文支撑）— 最可信
- `LithographySimulator`（quarterwave0）— 支持任意像差
- `high-na-euv-sim` — 含 Zernike 波前与 focus stack
- `openuv` / `euvsimulator` — 描述可疑（README 含 LLM 生成痕迹），**用前必须人工验证**

本期用解析 Bossung（§3.3）。接仿真属于 Phase 3 的可选升级，且必须在 `docs/experiments.md` 里记录验证过程。

### 11.3 无网络

运行期零网络调用。任何下载数据的脚本必须是显式、独立、可选的，且不进 `eumetria` 包。

---

## 12. 输出与可视化清单

| 图 | 内容 | 用于 |
|---|---|---|
| `field_layout.png` | 硅片 + field 网格 + 采样点 + edge exclusion 环 | 说明采样方案 |
| `overlay_quiver.png` | overlay 向量场（箭头）+ 残差向量 | 主结果 |
| `overlay_scalar_{x,y}.png` | dx / dy 的热图 | 分解前 |
| `zernike_terms.png` | 恢复系数 vs 真值的棒图（成对） | A-02/A-03 证据 |
| `residual_map.png` | 拟合后残差的空间分布 | 判断是否漏了系统性成分 |
| `sampling_budget.png` | 点数 vs 可辨识阶数（热图或折线） | Q-03 结论 |
| `crossval.png` | 各阶数下的 CV RMSE | Q-02 结论 |
| `bossung_family.png` | CD-focus 曲线族 | Phase 3 |
| `cdu_decomposition.png` | CD map = dose 场 + focus 场 + 残差（四联图） | Phase 3 主结果 |
| `report.html` | 上述图的汇总 + metrics 表 + 结论文字 | 交付/演示 |

命名规范：`{phase}_{name}.png`，全部落在 `out/<scenario>/report/`。

---

## 13. 测试策略

### 13.1 单元

- `noll_to_nm` 映射对 j=1..36 的已知对照
- 径向多项式 `R_n^m` 在 ρ=0、ρ=1 的解析值
- `to_polar` / `normalize_field` 在四象限与原点、边界值
- `build_layout`：默认场景断言确切值 —— `include_partial=True` → 76；`False` → 64。不接受范围断言
- `sample_fields(mode='regular')`：必须**按空间均匀**抽取（网格分块后每块取样），不得按列表索引等距抽取 —— 后者会因网格遍历顺序把点集中到边缘（实测踩过：25 点等距抽取得到 rho 全部 > 0.75）
- `field_marks` 三种密度的点数与对称性
- Bossung：已知参数下的极值位置、单调区间

### 13.2 性质测试（重点）

| 测试 | 断言 |
|---|---|
| 正交性 | `Z^T Z / n → I`，相对误差 < 1e-10（在充足均匀采样点上） |
| 恢复性（核心） | 生成 → 拟合 → 系数相对误差 < 1e-3（无噪声，点数 = 10× 项数） |
| 有噪声恢复 | σ = 0.1 nm 时主导项相对误差 < 10% |
| 往返换算 | `zernike → linear → zernike` 误差 < 1e-9 |
| 隔离性 | `models/*` 不 import `synth`（AST 扫描断言，§4.1） |
| 确定性 | 同 seed 两次运行输出差 < 1e-12 |
| 正交基等价性 | Zernike 基下 `fit_cascade` vs `fit_field` 系数差 < 1e-8（证明正交性 ⇒ 级联不必要，§7.9） |
| 非正交基对照 | Cartesian 基下同一对比，系数差显著（量级差异要记录进 experiments.md） |
| mark 级过滤 | 构造一个 mark 落在 edge exclusion 内的 field，断言该点被剔除且 rho_max ≤ 0.98 |
| 鲁棒性 | 手工注入 5% 离群点（10σ），LS 系数被带偏 ≥ 10×，MAD+Huber 解恢复误差 < 5% |
| CPE 可辨识性 | 构造 n_fields×6 / n_points > 0.3 的场景，断言 `cpe_ratio` 警告被触发且不抛出异常（§3.1.1） |
| 类型×粒度独立性 | 同一组测量数据分别用 wafer-level 与 CPE 两种粒度拟合，断言前者恢复的线性 6 参数与后者各场线性参数的均值一致（验证 §3.1.1 的正交性论述） |

### 13.3 边界

- 点数恰好等于参数数（应可解但条件数差，需给警告）
- 点数 = 参数数 − 1（抛 `InsufficientPointsError`）
- σ = 0；σ 极大（系数恢复必然失败，断言不崩溃且 rel_err 被正确报告）
- 单项指纹（只有 Z4）→ 恢复应完美
- 空 terms（全零指纹）→ 拟合结果应近似全零，残差 RMS ≈ σ
- `rho` 略大于 1（edge exclusion 处理）

### 13.4 回归

固定 seed（20261003）的完整流水线金标准：`truth.json` 内容哈希、`fit_result.json` 的关键数值、`metrics.json` 全字段。

### 13.5 集成

`eumetria run` 端到端：退出码 0、四个产物文件存在、`report.html` 含 8 张图引用。

---

## 14. 交付切片

每个切片是端到端可验证的，不是"写了一堆模块"。

### Phase 1 — 闭环（3–5 天）
**交付：** `geometry` + `basis`(Zernike) + `synth`(指纹/采样/噪声) + `fit`(固定项最小二乘) + `metrics` + `cli synth|fit|score` + 基础图
**验证：** A-01、A-02、A-04、A-05 通过；`eumetria run` 一条命令跑出残差 RMS。
**此时项目已可讲。**

### Phase 2 — 结论层（+2 天）
**交付：** `model_select`(AIC/逐步) + `crossval` + `fit_cascade`(正交性等价性实验) + `robust`(MAD/Huber) + `sampling_budget_curve` + `sweep` 命令
**验证：** Q-01、Q-02、Q-03、Q-05 有结论写入 `docs/experiments.md`；`sampling_budget.png` 产出；正交基等价性测试通过。
**此时项目有"我发现了什么"而不仅是"我实现了什么"。**

### Phase 2.5 — 统计推断层（+2 天，纯 Python）
**交付：** `inference`(FIM / CRLB / 系数协方差与置信区间) + `design`(D-optimal 贪心布局 + 序贯加点) + `sample-plan` 命令 + `inverse`(CDU→dose/focus, Bossung)
**验证：** A-09、A-10 通过；Q-06、Q-08 有结论（实测方差 vs CRLB 的差距、最优布局相对常规采样的增益）；`cdu_decomposition.png` 产出；Bossung 二义性按 §10.3 显式输出，禁止静默选支。
**此时项目从"我试出来了"变成"我算出了它不可能更好"。** 这是不碰 C++ 的情况下能达到的最高难度轴（D-08）。

### Phase 3 — 深度层：可校正性判定（+3–4 天）
**第二个卖点，与 Phase 2.5 平级，不砍。** 量测算法真正难的不是"分解出指纹"，而是"哪些成分该归扫描机校正、哪些该认命" —— 判错任何一个方向都是良率代价（漏校 = 白丢良率，过校 = 引入振荡）。这是本项目最能体现领域判断的一层。

**交付：** `decompose`(correctable / CPE / residual 三分 + 线性 6 参换算) + `significance`(每项的 t/F 检验、p 值与置信区间) + `budget`(扫描机支持项数 × CPE 每场阶数的敏感度扫描) + `tis`(校准残余误差的传播分析) + §3.5 量级表的自动断言
**验证：** Q-04、Q-07、Q-09、Q-10 有结论；`overfitting_crossover.png`（CPE 阶数 vs 拟合残差/CV 残差的交叉点）、`correctability_budget.png` 产出；`tests/property/test_magnitudes.py` 通过。
**此时项目能回答"这项该不该校"与"校到几阶为止"，并且给出代价曲线而不是一个数字。**

### Phase 4 — 门面（+1 天，砍到最小）
**交付：** `report`(HTML) + README（含 Q-01..Q-10 结论摘要 + 与文献量级的对比表 + **负面结论专节**）+ 诊断规则表
**验证：** 陌生人（或未来的自己）只读 README + report.html 就能理解项目做了什么、结论是什么、边界在哪。
**明确不做：** §7.10 的 `correction`(QP 校正优化) 不做（U-08）。Phase 4 是把已有的收干净，不是加功能。
**范围红线：** 预算紧张时压缩 Phase 4（报告层）与 Phase 2 的 sweep 网格密度，**不砍 Phase 2.5 与 Phase 3** —— 这两层是项目的两个卖点（统计推断 ≈ 算法深度，可校正性判定 ≈ 领域判断）。

---

## 15. 风险与已知局限

### 15.1 风险

| 风险 | 影响 | 缓解 |
|---|---|---|
| 合成数据的真实性被质疑 | 面试时说服力下降 | §3.5 参数表全部对标文献；README 明确列出假设与来源 |
| 自实现 Zernike 出错 | 全盘结果不可信 | 与 `zernpy` 对照测试（§13.2） |
| 反演问题病态被当成 bug | 误判为"实现有问题" | §10.3 显式设计二义性与警告输出 |
| 范围膨胀到图像处理/CNN | 偏离 metrology 主线 | §1.1 non-goals 明确禁止 |
| 采样预算结论不显著 | 亮点不成立 | 先跑 Q-03 的粗网格（5×5），若无可辨识区间则如实报告（"该布局下不可辨识"本身也是结论） |

### 15.2 明确不用的方法（及理由）

| 不做 | 理由 |
|---|---|
| CNN / 深度学习做指纹分解 | 该领域的方法是线性分解 + 统计推断；端到端 DL 会把可解释性（面试最看重的部分）全部丢掉 |
| 通用降维（PCA/UMAP）当主线 | 正交基分解的物理含义是 PCA 给不了的；PCA 可作为"对照基线"出现，不当主方法 |
| 用真实 fab 数据 | 不存在公开数据；且合规风险 |
| Web 后端包装 | 见 D-05；项目价值在算法结论，不在 API |

### 15.3 已知局限（要在 README 里明说）

1. 矢量特征（vector Zernike / 严格的正交矢量基）未实现，x/y 分量独立分解
2. 不做 reticle heating、lens heating 等时间依赖效应（静态场）
3. APC 只做系数语义与稳定性分析，不实现实际反馈控制
4. Bossung 用解析模型，未包含光刻胶化学放大与随机效应（Phase 3 可选升级）
5. 单硅片为主；lot 级统计只在 Phase 4 可选扩展
6. 未实现带约束的校正优化（§7.10 / U-08）：本项目只做拟合与可校正性分解，不产校正动作本身
7. 不做 native 加速层（D-08）：性能约束只到"10 万点 QR 求解可接受"（A-10），不追吞吐；工程性能不在本项目的主张范围内

---

## 附：实现顺序检查表

**进度、验收状态与下一步都在 [`docs/progress.md`](progress.md)** —— 本文档只放规格，
不放进度；数字的出处是 [`docs/experiments.md`](experiments.md)。

开工时按此顺序，每步都可独立验证：

- [x] `geometry.to_polar` / `normalize_field` + 测试
- [x] `basis.ZernikeBasis` + 正交性测试 + zernpy 对照（j=1..36 最大偏差 4.4e-16）
- [x] `basis` 设计矩阵 → 最小二乘能恢复已知系数（**第一个真正的里程碑**：96 点 / 12 项，相对误差 6.0e-13，条件数 2.3）
- [x] `geometry.build_layout` + `field_marks`（默认场景 76 field 口径B / 64 口径A；已实现 4 点与 8 点，19 点见 U-10）
- [x] `synth.make_fingerprint` 端到端：layout → 指纹（入参改为 `SamplePoints`，见 §7.3）
- [x] `synth.sample_fields` / `sample_points` + `add_noise` + `apply_tis`
      （默认场景复现 25 field / 188 点 / 4 个边缘场各剩 5 mark）
- [x] `io.schemas` + 读写（pydantic `extra='forbid'`；§6.2.5 四条强制校验；浮点写 `repr` 精确往返；
      version 在进 pydantic 之前检查，否则 `SchemaVersionError` 会被包成 `ValidationError`）
- [x] `fit.fit_field` **联合求解**（`X=[Z|C]` 一次 lstsq 后按块切开，D-09）+ `InsufficientPointsError`（显式 `n_points ≥ n_params` + 秩检查）
      + 条件数/相关系数警告 + 逐场秩诊断与 `cpe_order_ceiling`（§3.1.1 的落地）
- [x] `metrics.coefficient_recovery` + `residual_rms` + A-02/A-03 判定（A-02 要判自己的前提：
      `truth.noise_sigma_nm == 0`，否则带噪数据可能被 1% 阈值蒙过）
- [x] `cli synth|fit|score|run` + 退出码映射 + `--json`（stdout 干净 / 日志走 stderr）
      + 集成测试（A-01..A-05 逐条断言）
- [x] `report` 基础图（4 张：实测 dx 场 / 拟合 dx 场 / overlay 矢量 / 双通道残差图
      + 自包含单文件 HTML，图以 base64 内嵌；中文字体按可用性自动挑，避免画成豆腐块）
- [x] → **Phase 1 完成，README 已写**
- [ ] `model_select` + `crossval` + `sweep` + `sampling_budget_curve`（Q-01/02/03/05）
- [ ] `robust`(MAD 剔除) + `fit_cascade` 正交性等价性实验（Q-01 的数值证据）
- [ ] `inference.crlb` / `coefficient_ci` + A-09 覆盖率测试（500 次蒙特卡洛）
- [ ] `design.d_optimal` / `sequential` + `sample-plan` 命令 + Q-08 结论（实测方差 vs CRLB）
- [ ] `inverse`（Bossung，二义性显式输出）+ Q-06
- [ ] `decompose`（correctable / CPE / residual 三分 + 线性 6 参换算 + 每场秩/条件数报告 + CPE 阶数上限判定，见 §3.1.1）
- [ ] `significance`（每项 p 值/CI + 漏校-过校权衡曲线）+ Q-09
- [ ] `budget`（扫描机项数 × CPE 阶数扫描）+ Q-10 + `overfitting_crossover.png`
- [ ] `tis`（校准残余误差传播）+ Q-07
- [ ] `test_magnitudes.py`（§3.5 区间断言）+ Q-04
- [ ] → Phase 2/2.5/3 完成，README 汇总 Q-01..Q-10 + 负面结论专节
