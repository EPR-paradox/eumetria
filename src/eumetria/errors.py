"""eumetria 的异常层次。

统一在这里定义，便于 CLI 把异常映射成退出码（§8）：
输入/schema 错误 → 1；计算错误（点数不足、降秩等）→ 2；内部错误 → 3。
"""

from __future__ import annotations


class EumetriaError(Exception):
    """本项目的所有受控异常的基类。"""


class SchemaVersionError(EumetriaError, ValueError):
    """文件里的 `version` 与当前代码不兼容（§6.6）。"""


class SchemaValidationError(EumetriaError, ValueError):
    """结构/schema 校验失败：缺列、非法项名、越界取值等。"""


class InvalidMeasurementError(EumetriaError, ValueError):
    """`measurements.csv` 的内容自相矛盾（§6.2.5 的四条强制校验）。"""


class InsufficientPointsError(EumetriaError, ValueError):
    """点数不足或设计矩阵降秩：不给解，也不静默返回最小范数解（§3.1.1、§7.4）。"""

    def __init__(self, message: str, *, required: int, available: int) -> None:
        super().__init__(f"{message}（需要 ≥ {required} 点，实际 {available} 点）")
        self.required = required
        self.available = available
