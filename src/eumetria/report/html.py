"""单文件 HTML 报告（§7.8、P-05）。

图以 base64 内嵌，所以产出是**一个** HTML 文件 —— 发给别人、或者面试时打开都不用带 PNG。
摘要区只放判定与关键数字，不放解释：解释在 docs/ 里，报告是指向证据的指针。
"""

from __future__ import annotations

import base64
from collections.abc import Mapping
from html import escape
from pathlib import Path
from typing import Any

_STYLE = """
body { font-family: -apple-system, "Segoe UI", "Noto Sans CJK SC", sans-serif;
       max-width: 1100px; margin: 2rem auto; padding: 0 1rem; color: #1a1a1a; }
h1 { font-size: 1.5rem; } h2 { font-size: 1.15rem; margin-top: 2rem; }
table { border-collapse: collapse; margin: 0.5rem 0 1rem; }
th, td { border: 1px solid #d0d0d0; padding: 0.3rem 0.7rem; font-size: 0.9rem; }
th { background: #f4f4f4; text-align: left; }
.pass { color: #146c2e; font-weight: 600; } .fail { color: #a4262c; font-weight: 600; }
.na { color: #777; }
img { max-width: 100%; border: 1px solid #e0e0e0; margin: 0.5rem 0 1.5rem; }
.note { color: #555; font-size: 0.9rem; }
"""


def _fmt_flag(value: Any) -> str:
    if value is True:
        return '<span class="pass">通过</span>'
    if value is False:
        return '<span class="fail">未过</span>'
    return '<span class="na">—</span>'


def render_html(
    figures: list[Path],
    metrics: Mapping[str, Any],
    out_path: str | Path,
    *,
    title: str = "eumetria 报告",
    fit: Mapping[str, Any] | None = None,
) -> Path:
    """把图与判定写成一个自包含的 HTML。"""
    blocks: list[str] = [
        "<!doctype html>",
        '<html lang="zh"><head><meta charset="utf-8">',
        f"<title>{escape(title)}</title>",
        f"<style>{_STYLE}</style></head><body>",
        f"<h1>{escape(title)}</h1>",
    ]

    if fit is not None:
        ident = fit.get("identifiability", {})
        blocks += [
            "<h2>拟合摘要</h2><table>",
            f"<tr><th>点数 / 参数数</th><td>{fit.get('n_points')} / {fit.get('n_params')}</td></tr>",
            f"<tr><th>设计矩阵条件数</th><td>{fit.get('condition_number')}</td></tr>",
            f"<tr><th>残差 RMS (nm)</th><td>{fit.get('residual_rms_nm')}</td></tr>",
            f"<tr><th>可辨识最高阶</th><td>{ident.get('max_identifiable_order')}</td></tr>",
            f"<tr><th>CPE 阶数上限</th><td>{ident.get('cpe_order_ceiling')}</td></tr>",
            "</table>",
        ]

    blocks += [
        "<h2>验收判定</h2><table>",
        "<tr><th>条件</th><th>结果</th><th>数值</th></tr>",
        f"<tr><td>A-02（无噪声 + 点数 ≥ 8× 参数 → 全部系数 &lt; 1%）</td>"
        f"<td>{_fmt_flag(metrics.get('pass_A02'))}</td>"
        f"<td>最大相对误差 {metrics.get('max_rel_err')}"
        + ("" if metrics.get("a02_applicable") else "（不适用：数据带噪）")
        + "</td></tr>",
        f"<tr><td>A-03（前 6 阶主要项 &lt; 10%）</td>"
        f"<td>{_fmt_flag(metrics.get('pass_A03'))}</td>"
        f"<td>最大相对误差 {metrics.get('max_rel_err')}</td></tr>",
        f"<tr><td>A-04（残差 RMS 与真值差 &lt; 0.01nm）</td>"
        f"<td>{_fmt_flag(metrics.get('residual_agreement_nm') is not None and metrics['residual_agreement_nm'] < 0.01)}</td>"
        f"<td>{metrics.get('residual_agreement_nm')} nm</td></tr>",
        "</table>",
    ]

    notes = list(metrics.get("notes", []))
    if notes:
        blocks.append("<h2>说明</h2><ul>")
        blocks += [f'<li class="note">{escape(str(n))}</li>' for n in notes]
        blocks.append("</ul>")

    if fit is not None:
        warnings = list(fit.get("identifiability", {}).get("warnings", []))
        if warnings:
            blocks.append("<h2>可辨识性告警</h2><ul>")
            blocks += [f'<li class="note">{escape(str(w))}</li>' for w in warnings]
            blocks.append("</ul>")

    blocks.append("<h2>图</h2>")
    for figure in figures:
        path = Path(figure)
        if not path.exists():  # pragma: no cover - 防御
            continue
        encoded = base64.b64encode(path.read_bytes()).decode("ascii")
        blocks.append(f'<figure><img alt="{escape(path.stem)}" '
                      f'src="data:image/png;base64,{encoded}"><figcaption '
                      f'class="note">{escape(path.name)}</figcaption></figure>')

    blocks.append("</body></html>")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(blocks), encoding="utf-8")
    return out_path
