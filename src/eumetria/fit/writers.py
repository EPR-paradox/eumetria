"""写 `fit_result.json`（§6.4）。

放在 fit 包内而不是 io 包内，是为了守住依赖方向（§4.2）：io 不依赖 fit，
fit 依赖 io —— 反过来 import 会成环。

Phase 1 只写「拟合能给出的部分」；`linear_terms` / `correctable` / `cpe_terms` /
`uncorrectable_residual_rms_nm` 由 Phase 3 的 `decompose` 补齐。它们**整键不出现**，
而不是写成 0 —— 用 0 冒充"已知值"比缺字段危险得多。
"""

from __future__ import annotations

import json
from pathlib import Path

from .result import FitResult


def write_fit_result(path: str | Path, result: FitResult) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(result.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return path
