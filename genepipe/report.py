"""
genepipe.report -- automatic report generation.

An automated pipeline is only useful if its results are readable without a
programmer. This module renders the run manifest and every stage's outputs into
both Markdown (for the repository) and a self-contained HTML page (for reading),
and writes the machine-readable JSON manifest that the reproducibility record
depends on.
"""
from __future__ import annotations

import html
import json
import os
from datetime import datetime
from typing import Any

from .utils import __version__, get_logger

log = get_logger("report")

CSS = """
body{font-family:Georgia,serif;max-width:960px;margin:2rem auto;color:#1b1b1b;line-height:1.5}
h1{color:#2F4F6F;border-bottom:2px solid #2F4F6F;padding-bottom:.3rem}
h2{color:#2F4F6F;margin-top:1.6rem}
table{border-collapse:collapse;width:100%;margin:.8rem 0;font-size:.9rem}
th{background:#2F4F6F;color:#fff;text-align:left;padding:.4rem .5rem}
td{border:1px solid #cfd8e3;padding:.35rem .5rem}
tr:nth-child(even) td{background:#f4f7fa}
code,pre{font-family:"DejaVu Sans Mono",monospace;background:#f2f4f7;padding:.1rem .3rem;border-radius:3px}
pre{padding:.8rem;overflow-x:auto;font-size:.85rem}
.ok{color:#2E7D5B;font-weight:bold}.fail{color:#B03030;font-weight:bold}.skip{color:#C9762B;font-weight:bold}
.meta{color:#666;font-size:.85rem}
"""


def _table(headers: list[str], rows: list[list[Any]]) -> str:
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(c))}</td>" for c in r) + "</tr>"
                   for r in rows)
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def _md_table(headers: list[str], rows: list[list[Any]]) -> str:
    out = ["| " + " | ".join(str(h) for h in headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def _status_span(status: str) -> str:
    return f'<span class="{status}">{status}</span>'


def write_manifest(manifest: dict, path: str) -> str:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2, default=str)
    log.info("wrote run manifest -> %s", path)
    return path


def build_markdown(manifest: dict) -> str:
    """Render the run manifest as Markdown."""
    md = [f"# genePipe run report",
          "",
          f"*Generated {manifest['finished_at']} by genePipe v{manifest['environment']['pipeline_version']}*",
          "",
          "## Run configuration",
          "",
          f"- **Input:** `{manifest['input']}`",
          f"- **Output directory:** `{manifest['output_dir']}`",
          f"- **Reference sequence:** `{manifest.get('reference_id') or 'auto (longest sequence)'}`",
          f"- **Alignment:** {manifest.get('alignment_method', 'n/a')}",
          f"- **Total runtime:** {manifest['total_seconds']} s",
          "",
          "## Stage summary",
          "",
          _md_table(["Stage", "Status", "Seconds", "Notes"],
                    [[s["name"], s["status"], s["seconds"], s["message"]]
                     for s in manifest["stages"]]),
          ""]

    qc = manifest.get("stages_metrics", {}).get("qc", {})
    if qc:
        md += ["## Quality control", "",
               _md_table(["Metric", "Value"], [[k, v] for k, v in qc.items() if not isinstance(v, dict)]),
               ""]
    vs = manifest.get("stages_metrics", {}).get("variants", {})
    if vs:
        md += ["## Variant summary", "",
               _md_table(["Metric", "Value"], [[k, v] for k, v in vs.items()]), ""]
    assoc = manifest.get("stages_metrics", {}).get("association", {})
    if assoc.get("n_significant") is not None:
        md += ["## Association testing", "",
               _md_table(["Metric", "Value"],
                         [["Sites tested", assoc.get("n_tested")],
                          ["Significant at FDR 0.05", assoc.get("n_significant")],
                          ["Most significant site", assoc.get("top_site") or "none"]]), ""]
    sup = manifest.get("stages_metrics", {}).get("bootstrap", {})
    if sup:
        md += ["## Phylogenetic support", "",
               _md_table(["Clade", "Support (%)"],
                         [[k, v] for k, v in sup.items()]), ""]
    md += ["## Artefacts", ""]
    for path_ in manifest.get("artefacts", []):
        md.append(f"- `{os.path.relpath(path_, manifest['output_dir'])}`")
    return "\n".join(md) + "\n"


def build_html(manifest: dict) -> str:
    """Render the run manifest as a standalone HTML page."""
    parts = [f"<!doctype html><html><head><meta charset='utf-8'>",
             f"<title>genePipe run report</title><style>{CSS}</style></head><body>",
             "<h1>genePipe run report</h1>",
             f"<p class='meta'>Generated {html.escape(str(manifest['finished_at']))} by "
             f"genePipe v{manifest['environment']['pipeline_version']} on "
             f"{html.escape(manifest['environment']['platform'])}</p>",
             "<h2>Run configuration</h2>",
             _table(["Setting", "Value"],
                    [["Input", manifest["input"]],
                     ["Output directory", manifest["output_dir"]],
                     ["Reference sequence", manifest.get("reference_id") or "auto (longest sequence)"],
                     ["Alignment method", manifest.get("alignment_method", "n/a")],
                     ["Bootstrap replicates", manifest.get("bootstrap_replicates", "n/a")],
                     ["Total runtime (s)", manifest["total_seconds"]]]),
             "<h2>Stage summary</h2>",
             _table(["Stage", "Status", "Seconds", "Notes"],
                    [[s["name"], _status_span(s["status"]), s["seconds"], s["message"] or "&mdash;"]
                     for s in manifest["stages"]])]

    metrics = manifest.get("stages_metrics", {})
    if metrics.get("qc"):
        parts += ["<h2>Quality control</h2>",
                  _table(["Metric", "Value"],
                         [[k, v] for k, v in metrics["qc"].items() if not isinstance(v, dict)])]
    if metrics.get("variants"):
        parts += ["<h2>Variant summary</h2>",
                  _table(["Metric", "Value"], [[k, v] for k, v in metrics["variants"].items()])]
    assoc = metrics.get("association", {})
    if assoc.get("n_tested") is not None:
        parts += ["<h2>Association testing</h2>",
                  _table(["Metric", "Value"],
                         [["SNV sites tested", assoc.get("n_tested")],
                          ["Significant at FDR 0.05", assoc.get("n_significant")],
                          ["Most significant site", assoc.get("top_site") or "none"]])]
    sup = metrics.get("bootstrap", {})
    if sup:
        parts += ["<h2>Phylogenetic support</h2>",
                  _table(["Clade", "Bootstrap support (%)"], [[k, v] for k, v in sup.items()])]
    parts += ["<h2>Artefacts</h2><ul>"]
    for p in manifest.get("artefacts", []):
        parts.append(f"<li><code>{html.escape(os.path.relpath(p, manifest['output_dir']))}</code></li>")
    parts.append("</ul></body></html>")
    return "\n".join(parts)


def write_reports(manifest: dict, output_dir: str) -> list[str]:
    os.makedirs(output_dir, exist_ok=True)
    paths = []
    md_path = os.path.join(output_dir, "run_report.md")
    with open(md_path, "w", encoding="utf-8") as fh:
        fh.write(build_markdown(manifest))
    paths.append(md_path)
    html_path = os.path.join(output_dir, "run_report.html")
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(build_html(manifest))
    paths.append(html_path)
    log.info("wrote run reports -> %s", ", ".join(os.path.basename(p) for p in paths))
    return paths
