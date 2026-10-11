# eumetria

光刻**量测**侧的 wafer 场分析工具：把一片硅片上的量测数据分解为可解释、可判定、可诊断的成分。

主线是 **overlay**（场分解 → 可校正性 → 采样设计 → 诊断）；focus 复用同一套框架以证明框架可复用；CDU 只做 CD→(dose, focus) 反演一条。

`eumetria` = εὐμετρία：匀称、比例恰当的度量。

## 当前状态

**设计阶段。代码尚未开始**，仓库里目前只有设计文档，没有可运行的实现。

- 详细设计：[`docs/software-design-v0.1.md`](docs/software-design-v0.1.md)（Draft v0.2）
  - 验收条件 A-01..A-11、技术问题 Q-01..Q-10、交付切片 Phase 1..4、实现顺序检查表

## 非目标

不做缺陷检测/分类（属 inspection）、不接真实 fab 数据、不做 Web 服务、不引入 C++/native 层。
