from __future__ import annotations

import html
import json
import math
import re
import statistics
import textwrap
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


@dataclass
class CommandRecord:
    command: str
    timestamp: str
    step: str | None = None


@dataclass
class LengthStats:
    count: int = 0
    min_length: int | None = None
    max_length: int | None = None
    mean_length: float | None = None
    median_length: float | None = None
    histogram: dict[int, int] = field(default_factory=dict)
    bin_width: int = 5

    def to_dict(self) -> dict:
        return {
            "count": self.count,
            "min": self.min_length,
            "max": self.max_length,
            "mean": self.mean_length,
            "median": self.median_length,
            "histogram": self.histogram,
            "bin_width": self.bin_width,
        }


@dataclass
class StepRecord:
    name: str
    file_key: str
    filename: str
    path: str
    count: int | None
    exists: bool
    length_stats: LengthStats | None = None
    domain_counts: dict[str, int] = field(default_factory=dict)
    phylum_counts: dict[str, int] = field(default_factory=dict)


@dataclass
class BuildReport:
    dataset_name: str
    primer_name: str
    working_dir: str
    started_at: str
    suffix: str | None = None
    build_stem: str = ""
    finished_at: str | None = None
    elapsed_seconds: float | None = None
    dry_run: bool = False
    environment: str | None = None
    primer_fwd: str = ""
    primer_rev: str = ""
    primer_min_length: int | None = None
    primer_max_length: int | None = None
    primer_max_n: int | None = None
    primer_mismatch: int | None = None
    primer_pga_percid: float | None = None
    input_files: list[str] = field(default_factory=list)
    commands: list[CommandRecord] = field(default_factory=list)
    steps: list[StepRecord] = field(default_factory=list)
    report_path: str | None = None


PIPELINE_STEPS: list[tuple[str, str, bool]] = [
    ("Import", "dataset_crabs_file", True),
    ("In-silico PCR (prefilter)", "pcr_prefilter_file", True),
    ("In-silico PCR (length filter)", "pcr_file", True),
    ("PGA", "pga_file", False),
    ("Dereplicate", "dereplicate_file", True),
    ("Filter", "filtered_file", True),
    ("SINTAX export", "sintax_file", True),
]


def count_file_records(path: Path) -> int:
    if not path.exists():
        return 0
    if path.suffix == ".fasta":
        count = 0
        with path.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    count += 1
        return count
    result = 0
    with path.open() as handle:
        for line in handle:
            if line.strip():
                result += 1
    return result


def analyze_crabs_txt(path: Path, bin_width: int = 5) -> tuple[LengthStats, dict[str, int], dict[str, int]]:
    stats = LengthStats(bin_width=bin_width)
    domain_counts: dict[str, int] = {}
    phylum_counts: dict[str, int] = {}
    lengths: list[int] = []

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            sequence = parts[10]
            length = len(sequence)
            lengths.append(length)
            stats.histogram[(length // bin_width) * bin_width] = (
                stats.histogram.get((length // bin_width) * bin_width, 0) + 1
            )
            domain = parts[3] if parts[3] not in ("", "NA") else "unassigned"
            phylum = parts[4] if parts[4] not in ("", "NA") else "unassigned"
            domain_counts[domain] = domain_counts.get(domain, 0) + 1
            phylum_counts[phylum] = phylum_counts.get(phylum, 0) + 1

    if not lengths:
        return stats, domain_counts, phylum_counts

    stats.count = len(lengths)
    stats.min_length = min(lengths)
    stats.max_length = max(lengths)
    stats.mean_length = statistics.mean(lengths)
    stats.median_length = statistics.median(lengths)
    return stats, domain_counts, phylum_counts


def analyze_fasta(path: Path, bin_width: int = 5) -> LengthStats:
    stats = LengthStats(bin_width=bin_width)
    lengths: list[int] = []
    sequence_parts: list[str] = []

    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if sequence_parts:
                    length = len("".join(sequence_parts))
                    lengths.append(length)
                    stats.histogram[(length // bin_width) * bin_width] = (
                        stats.histogram.get((length // bin_width) * bin_width, 0) + 1
                    )
                    sequence_parts = []
            else:
                sequence_parts.append(line.strip())
        if sequence_parts:
            length = len("".join(sequence_parts))
            lengths.append(length)
            stats.histogram[(length // bin_width) * bin_width] = (
                stats.histogram.get((length // bin_width) * bin_width, 0) + 1
            )

    if not lengths:
        return stats

    stats.count = len(lengths)
    stats.min_length = min(lengths)
    stats.max_length = max(lengths)
    stats.mean_length = statistics.mean(lengths)
    stats.median_length = statistics.median(lengths)
    return stats


def top_counts(counts: dict[str, int], limit: int = 15) -> list[tuple[str, int]]:
    return sorted(counts.items(), key=lambda item: (-item[1], item[0]))[:limit]


def step_count(steps: list[StepRecord], file_key: str) -> int | None:
    step = next((item for item in steps if item.file_key == file_key), None)
    return step.count if step else None


def import_count(steps: list[StepRecord]) -> int | None:
    return step_count(steps, "dataset_crabs_file")


def final_sequence_count(steps: list[StepRecord]) -> int | None:
    for file_key in ("sintax_file", "filtered_file", "dereplicate_file"):
        step = next((item for item in steps if item.file_key == file_key and item.count is not None), None)
        if step is not None:
            return step.count
    return steps[-1].count if steps else None


def retention_pct(current: int | None, previous: int | None) -> float | None:
    if current is None or previous is None or previous == 0:
        return None
    return 100.0 * current / previous


def step_note(step: StepRecord) -> str:
    if step.file_key == "pga_file":
        return "PGA expands amplicons via global alignment, so counts usually increase."
    return ""


def human_int(value: int | None) -> str:
    if value is None:
        return "—"
    return f"{value:,}"


def human_float(value: float | None, digits: int = 1) -> str:
    if value is None:
        return "—"
    return f"{value:.{digits}f}"


def format_axis_count(value: float) -> str:
    value = int(round(value))
    if value >= 1_000_000:
        return f"{value / 1_000_000:.1f}M"
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    return str(value)



def select_length_ticks(min_bp: int, max_bp: int, max_ticks: int = 6) -> list[int]:
    if max_bp <= min_bp:
        return [min_bp]
    span = max_bp - min_bp
    raw_step = span / (max_ticks - 1)
    if raw_step <= 5:
        step = 5
    elif raw_step <= 10:
        step = 10
    elif raw_step <= 25:
        step = 25
    elif raw_step <= 50:
        step = 50
    elif raw_step <= 100:
        step = 100
    elif raw_step <= 250:
        step = 250
    elif raw_step <= 500:
        step = 500
    else:
        step = int(10 ** round(math.log10(raw_step)))

    start = ((min_bp + step - 1) // step) * step
    ticks = []
    value = start
    while value <= max_bp and len(ticks) < max_ticks:
        ticks.append(value)
        value += step
    if ticks and ticks[0] != min_bp:
        ticks.insert(0, min_bp)
    if ticks and ticks[-1] != max_bp:
        ticks.append(max_bp)
    return ticks[:max_ticks] if len(ticks) > max_ticks else ticks


def render_histogram_svg(stats: LengthStats, title: str, width: int = 1048, height: int = 320) -> str:
    if not stats.histogram:
        return f'<p class="muted">No length data for {html.escape(title)}.</p>'

    bins = sorted(stats.histogram.items())
    max_count = max(count for _, count in bins)
    margin_left, margin_right, margin_top, margin_bottom = 56, 20, 24, 52
    chart_width = width - margin_left - margin_right
    chart_height = height - margin_top - margin_bottom
    bar_gap = 2
    bar_width = max(1, (chart_width - bar_gap * (len(bins) - 1)) / len(bins))
    chart_bottom = margin_top + chart_height

    bars = []
    for index, (bin_start, count) in enumerate(bins):
        bar_height = 0 if max_count == 0 else (count / max_count) * chart_height
        x = margin_left + index * (bar_width + bar_gap)
        y = chart_bottom - bar_height
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" '
            f'class="bar" data-count="{count}" data-bin="{bin_start}">'
            f"<title>{bin_start}-{bin_start + stats.bin_width - 1} bp: {count:,}</title></rect>"
        )

    y_tick_values = [0]
    if max_count > 0:
        for fraction in (0.25, 0.5, 0.75, 1.0):
            y_tick_values.append(max_count * fraction)
    y_ticks = []
    for tick_value in y_tick_values:
        tick_y = chart_bottom - (0 if max_count == 0 else (tick_value / max_count) * chart_height)
        y_ticks.append(
            f'<line x1="{margin_left}" y1="{tick_y:.1f}" x2="{width - margin_right}" y2="{tick_y:.1f}" class="grid" />'
            f'<line x1="{margin_left - 4}" y1="{tick_y:.1f}" x2="{margin_left}" y2="{tick_y:.1f}" class="axis" />'
            f'<text x="{margin_left - 8}" y="{tick_y + 3:.1f}" text-anchor="end" class="tick-label">'
            f"{html.escape(format_axis_count(tick_value))}</text>"
        )

    x_ticks = []
    if stats.min_length is not None and stats.max_length is not None:
        length_span = max(stats.max_length - stats.min_length, 1)
        for tick_bp in select_length_ticks(stats.min_length, stats.max_length):
            x = margin_left + ((tick_bp - stats.min_length) / length_span) * chart_width
            x_ticks.append(
                f'<line x1="{x:.1f}" y1="{chart_bottom}" x2="{x:.1f}" y2="{chart_bottom + 4}" class="axis" />'
                f'<text x="{x:.1f}" y="{chart_bottom + 16}" text-anchor="middle" class="tick-label">'
                f"{html.escape(str(tick_bp))}</text>"
            )

    return f"""
    <div class="chart-block">
      <h3>{html.escape(title)}</h3>
      <p class="muted">
        n={human_int(stats.count)} |
        min={human_int(stats.min_length)} |
        median={human_float(stats.median_length)} |
        mean={human_float(stats.mean_length)} |
        max={human_int(stats.max_length)} bp
      </p>
      <svg viewBox="0 0 {width} {height}" class="histogram" role="img" aria-label="{html.escape(title)}">
        {''.join(y_ticks)}
        <line x1="{margin_left}" y1="{margin_top}" x2="{margin_left}" y2="{chart_bottom}" class="axis" />
        <line x1="{margin_left}" y1="{chart_bottom}" x2="{width - margin_right}" y2="{chart_bottom}" class="axis" />
        {''.join(bars)}
        {''.join(x_ticks)}
        <text x="{(margin_left + width - margin_right) / 2:.1f}" y="{height - 8}" text-anchor="middle" class="axis-label">Sequence length (bp)</text>
        <text x="14" y="{(margin_top + chart_bottom) / 2:.1f}" text-anchor="middle" class="axis-label" transform="rotate(-90 14 {(margin_top + chart_bottom) / 2:.1f})">Count</text>
      </svg>
    </div>
    """


def render_taxonomy_table(title: str, counts: dict[str, int], total: int) -> str:
    if not counts:
        return f'<p class="muted">No taxonomy data for {html.escape(title)}.</p>'

    rows = []
    for label, count in top_counts(counts):
        pct = 100.0 * count / total if total else 0
        rows.append(
            f"<tr><td>{html.escape(label)}</td><td class='num'>{count:,}</td>"
            f"<td class='num'>{pct:.1f}%</td>"
            f"<td><div class='inline-bar' style='width:{pct:.1f}%'></div></td></tr>"
        )

    return f"""
    <div class="table-block">
      <h3>{html.escape(title)}</h3>
      <table>
        <thead><tr><th>Taxon</th><th>Count</th><th>%</th><th></th></tr></thead>
        <tbody>{''.join(rows)}</tbody>
      </table>
    </div>
    """


def human_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.0f} s"
    if seconds < 3600:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.1f} h"


def render_funnel_table(steps: list[StepRecord]) -> str:
    rows = []
    import_baseline = import_count(steps)
    for step in steps:
        from_import = retention_pct(step.count, import_baseline)
        note = step_note(step)
        rows.append(
            "<tr>"
            f"<td>{html.escape(step.name)}{f'<div class=\"muted\">{html.escape(note)}</div>' if note else ''}</td>"
            f"<td><code>{html.escape(step.filename)}</code></td>"
            f"<td class='num'>{human_int(step.count)}</td>"
            f"<td class='num'>{'' if from_import is None or step.file_key == 'dataset_crabs_file' else f'{from_import:.1f}%'}</td>"
            "</tr>"
        )

    import_note = (
        "<p class='muted'>Counts after in-silico PCR are amplicons, not source records. "
        "PGA can expand one source hit into multiple aligned amplicons, so later steps may exceed the PCR count.</p>"
        if import_baseline is not None
        else "<p class='muted'>Import file not found — retention percentages are unavailable.</p>"
    )
    return f"""
    {import_note}
    <table class="funnel">
      <thead>
        <tr>
          <th>Step</th>
          <th>Output file</th>
          <th>Records</th>
          <th>% of import</th>
        </tr>
      </thead>
      <tbody>{''.join(rows)}</tbody>
    </table>
    """


def format_command(command: str) -> str:
    cmd = " ".join(textwrap.dedent(command).split())
    return re.sub(r" (--)", r" \\\n\1", cmd)


def render_commands(commands: list[CommandRecord]) -> str:
    if not commands:
        return "<p class='muted'>No commands recorded. Run <code>builder.build()</code> to capture CRABS commands.</p>"

    grouped: dict[str, list[CommandRecord]] = {}
    for record in commands:
        step = record.step or "Other"
        grouped.setdefault(step, []).append(record)

    sections = []
    for step_name, step_commands in grouped.items():
        items = []
        for index, record in enumerate(step_commands, start=1):
            items.append(
                f"<div class='command-block'>"
                f"<div class='command-label'>Command {index} <span class='muted'>{html.escape(record.timestamp)}</span></div>"
                f"<pre>{html.escape(record.command)}</pre>"
                f"</div>"
            )
        sections.append(
            f"<div class='command-step'><h3>{html.escape(step_name)}</h3>{''.join(items)}</div>"
        )
    return "".join(sections)


def render_html(report: BuildReport) -> str:
    imported_count = import_count(report.steps)
    final_count = final_sequence_count(report.steps)
    pcr_count = step_count(report.steps, "pcr_file")
    pcr_prefilter_step = next((step for step in report.steps if step.file_key == "pcr_prefilter_file"), None)
    pcr_step = next((step for step in report.steps if step.file_key == "pcr_file"), None)
    filtered_step = next((step for step in report.steps if step.file_key == "filtered_file"), None)

    charts = []
    if pcr_prefilter_step and pcr_prefilter_step.length_stats:
        charts.append(render_histogram_svg(pcr_prefilter_step.length_stats, "Amplicon length after in-silico PCR"))
    if pcr_step and pcr_step.length_stats:
        charts.append(render_histogram_svg(pcr_step.length_stats, "Amplicon length after PCR length filter"))
    if filtered_step and filtered_step.length_stats:
        charts.append(render_histogram_svg(filtered_step.length_stats, "Amplicon length in final filtered set"))

    taxonomy = ""
    if filtered_step and filtered_step.phylum_counts and filtered_step.count:
        taxonomy = render_taxonomy_table(
            "Top phyla in final filtered set",
            filtered_step.phylum_counts,
            filtered_step.count,
        )

    display_name = report.build_stem or f"{report.dataset_name}_{report.primer_name}"
    metadata = {
        "dataset": report.dataset_name,
        "primer_set": report.primer_name,
        "suffix": report.suffix,
        "build_stem": display_name,
        "working_dir": report.working_dir,
        "started_at": report.started_at,
        "finished_at": report.finished_at,
        "elapsed_seconds": report.elapsed_seconds,
        "dry_run": report.dry_run,
        "environment": report.environment,
        "input_files": report.input_files,
        "imported_sequences": imported_count,
        "final_sequences": final_count,
    }

    overall_yield = (
        f"{100.0 * final_count / imported_count:.1f}%"
        if imported_count not in (None, 0) and final_count is not None
        else "—"
    )

    build_time_metric = ""
    if report.elapsed_seconds is not None and not report.dry_run:
        build_time_metric = f"""
        <div class="metric">
          <div class="label">Build time</div>
          <div class="value">{html.escape(human_duration(report.elapsed_seconds))}</div>
        </div>"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Reference DB build report — {html.escape(display_name)}</title>
  <style>
    :root {{
      --bg: #f6f8fb;
      --card: #ffffff;
      --text: #1f2937;
      --muted: #6b7280;
      --border: #e5e7eb;
      --accent: #2563eb;
      --accent-soft: #dbeafe;
      --bar: #3b82f6;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--text);
      line-height: 1.5;
    }}
    main {{
      max-width: 1100px;
      margin: 0 auto;
      padding: 24px;
    }}
    h1, h2, h3 {{ margin: 0 0 12px; }}
    h1 {{ font-size: 1.8rem; }}
    h2 {{ font-size: 1.2rem; margin-top: 28px; }}
    .hero, .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 20px;
      margin-bottom: 16px;
    }}
    .hero-grid, .meta-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin-top: 16px;
    }}
    .metric {{
      background: var(--accent-soft);
      border-radius: 10px;
      padding: 12px;
    }}
    .metric .label {{
      color: var(--muted);
      font-size: 0.85rem;
    }}
    .metric .value {{
      font-size: 1.4rem;
      font-weight: 700;
    }}
    .muted {{ color: var(--muted); }}
    table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.95rem;
    }}
    th, td {{
      border-bottom: 1px solid var(--border);
      padding: 10px 8px;
      text-align: left;
      vertical-align: top;
    }}
    th {{ color: var(--muted); font-weight: 600; }}
    td.num {{ text-align: right; white-space: nowrap; }}
    code {{
      background: #f3f4f6;
      padding: 2px 6px;
      border-radius: 6px;
      font-size: 0.85rem;
      word-break: break-all;
    }}
    pre {{
      white-space: pre-wrap;
      word-break: break-word;
      background: #111827;
      color: #f9fafb;
      padding: 12px;
      border-radius: 8px;
      overflow-x: auto;
      font-size: 0.85rem;
    }}
    details {{
      border: 1px solid var(--border);
      border-radius: 8px;
      padding: 10px 12px;
      margin-bottom: 8px;
      background: #fafafa;
    }}
    summary {{ cursor: pointer; font-weight: 600; }}
    .command-block {{
      margin-bottom: 12px;
    }}
    .command-label {{
      font-size: 0.9rem;
      font-weight: 600;
      margin-bottom: 6px;
    }}
    .charts {{
      display: flex;
      flex-direction: column;
      gap: 16px;
    }}
    .chart-block, .table-block {{
      width: 100%;
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 16px;
    }}
    .command-step {{
      margin-bottom: 16px;
    }}
    .command-step h3 {{
      font-size: 1rem;
      margin-bottom: 8px;
    }}
    .histogram {{
      width: 100%;
      height: auto;
      background: #fbfdff;
      border: 1px solid var(--border);
      border-radius: 8px;
    }}
    .axis {{ stroke: #6b7280; stroke-width: 1; }}
    .grid {{ stroke: #e5e7eb; stroke-width: 1; }}
    .tick-label {{
      font-size: 10px;
      fill: #6b7280;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    .axis-label {{
      font-size: 11px;
      fill: #374151;
      font-weight: 600;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    .bar {{ fill: var(--bar); opacity: 0.85; }}
    .inline-bar {{
      height: 10px;
      background: var(--bar);
      border-radius: 999px;
      min-width: 2px;
    }}
    .primer-table td:first-child {{ width: 180px; color: var(--muted); }}
  </style>
</head>
<body>
  <main>
    <section class="hero">
      <h1>Reference database build report</h1>
      <p class="muted">
        {html.escape(display_name)}
        · generated {html.escape(report.finished_at or report.started_at)}
      </p>
      <div class="hero-grid">
        <div class="metric">
          <div class="label">Imported source sequences</div>
          <div class="value">{human_int(imported_count)}</div>
        </div>
        <div class="metric">
          <div class="label">Final reference amplicons</div>
          <div class="value">{human_int(final_count)}</div>
        </div>
        <div class="metric">
          <div class="label">Yield (final / import)</div>
          <div class="value">{overall_yield}</div>
        </div>
        <div class="metric">
          <div class="label">PCR amplicons</div>
          <div class="value">{human_int(pcr_count)}</div>
        </div>{build_time_metric}
      </div>
    </section>

    <section class="card">
      <h2>Primer set</h2>
      <table class="primer-table">
        <tr><td>Forward</td><td><code>{html.escape(report.primer_fwd)}</code></td></tr>
        <tr><td>Reverse</td><td><code>{html.escape(report.primer_rev)}</code></td></tr>
        <tr><td>Length range</td><td>{human_int(report.primer_min_length)} – {human_int(report.primer_max_length)} bp</td></tr>
        <tr><td>Max N</td><td>{html.escape(str(report.primer_max_n))}</td></tr>
        <tr><td>PCR mismatch</td><td>{html.escape(str(report.primer_mismatch))}</td></tr>
        <tr><td>PGA % identity</td><td>{html.escape(str(report.primer_pga_percid))}</td></tr>
        <tr><td>Input files</td><td>{''.join(f'<div><code>{html.escape(path)}</code></div>' for path in report.input_files) or '<span class="muted">—</span>'}</td></tr>
        <tr><td>Working directory</td><td><code>{html.escape(report.working_dir)}</code></td></tr>
        <tr><td>Environment</td><td>{html.escape(report.environment or 'default')}</td></tr>
        <tr><td>Dry run</td><td>{'yes' if report.dry_run else 'no'}</td></tr>
      </table>
    </section>

    <section class="card">
      <h2>Pipeline retention</h2>
      {render_funnel_table(report.steps)}
    </section>

    <section class="charts">
      {''.join(charts)}
    </section>

    {taxonomy}

    <section class="card">
      <h2>Commands</h2>
      {render_commands(report.commands)}
    </section>

    <section class="card">
      <h2>Raw metadata</h2>
      <pre>{html.escape(json.dumps(metadata, indent=2))}</pre>
    </section>
  </main>
</body>
</html>
"""


class BuildReportCollector:
    def __init__(self, working_dir: str, dry_run: bool = False, environment: str | None = None):
        self.working_dir = working_dir
        self.dry_run = dry_run
        self.environment = environment
        self._current_step: str | None = None
        self.report = BuildReport(
            dataset_name="",
            primer_name="",
            working_dir=working_dir,
            started_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
            dry_run=dry_run,
            environment=environment,
        )
        self._build_started_monotonic: float | None = None

    def begin_build(self) -> None:
        self._build_started_monotonic = time.monotonic()
        self.report.started_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        self.report.finished_at = None
        self.report.elapsed_seconds = None

    def finish_build(self) -> None:
        if self._build_started_monotonic is None:
            return
        self.report.finished_at = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        self.report.elapsed_seconds = time.monotonic() - self._build_started_monotonic
        self._build_started_monotonic = None

    def set_database(self, database) -> None:
        dataset = database.dataset
        primer = database.primer_set
        self.report.dataset_name = dataset.name
        self.report.primer_name = primer.name
        self.report.suffix = database.suffix
        self.report.build_stem = database.build_stem
        self.report.primer_fwd = primer.fwd
        self.report.primer_rev = primer.rev
        self.report.primer_min_length = primer.min_length
        self.report.primer_max_length = primer.max_length
        self.report.primer_max_n = primer.max_n
        self.report.primer_mismatch = primer.mismatch
        self.report.primer_pga_percid = primer.pga_percid
        if dataset.files:
            self.report.input_files = [file.path for file in dataset.files]

    @contextmanager
    def step(self, step_name: str):
        previous = self._current_step
        self._current_step = step_name
        try:
            yield
        finally:
            self._current_step = previous

    def record_command(self, command: str, step: str | None = None) -> None:
        self.report.commands.append(
            CommandRecord(
                command=format_command(command),
                timestamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC"),
                step=step or self._current_step,
            )
        )

    def manifest_path(self) -> Path:
        stem = self.report.build_stem or f"{self.report.dataset_name}_{self.report.primer_name}"
        return Path(self.working_dir) / f"{stem}_build_manifest.json"

    def save_manifest(self) -> None:
        payload = {
            "dataset_name": self.report.dataset_name,
            "primer_name": self.report.primer_name,
            "suffix": self.report.suffix,
            "build_stem": self.report.build_stem,
            "started_at": self.report.started_at,
            "finished_at": self.report.finished_at,
            "elapsed_seconds": self.report.elapsed_seconds,
            "dry_run": self.report.dry_run,
            "environment": self.report.environment,
            "primer_fwd": self.report.primer_fwd,
            "primer_rev": self.report.primer_rev,
            "primer_min_length": self.report.primer_min_length,
            "primer_max_length": self.report.primer_max_length,
            "primer_max_n": self.report.primer_max_n,
            "primer_mismatch": self.report.primer_mismatch,
            "primer_pga_percid": self.report.primer_pga_percid,
            "input_files": self.report.input_files,
            "commands": [
                {
                    "command": record.command,
                    "timestamp": record.timestamp,
                    "step": record.step,
                }
                for record in self.report.commands
            ],
        }
        self.manifest_path().write_text(json.dumps(payload, indent=2), encoding="utf-8")

    def load_manifest(self) -> bool:
        path = self.manifest_path()
        if not path.exists():
            return False
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.report.started_at = payload.get("started_at", self.report.started_at)
        self.report.finished_at = payload.get("finished_at")
        self.report.elapsed_seconds = payload.get("elapsed_seconds")
        self.report.commands = [
            CommandRecord(
                command=item["command"],
                timestamp=item["timestamp"],
                step=item.get("step"),
            )
            for item in payload.get("commands", [])
        ]
        return bool(self.report.commands)

    def record_steps(self, files: dict[str, str]) -> None:
        self.report.steps = []
        for step_name, file_key, analyze in PIPELINE_STEPS:
            filename = files[file_key]
            path = Path(self.working_dir) / filename
            exists = path.exists()
            count = count_file_records(path) if exists else None
            length_stats = None
            domain_counts: dict[str, int] = {}
            phylum_counts: dict[str, int] = {}

            if exists and analyze and count:
                if path.suffix == ".fasta":
                    length_stats = analyze_fasta(path)
                else:
                    length_stats, domain_counts, phylum_counts = analyze_crabs_txt(path)

            self.report.steps.append(
                StepRecord(
                    name=step_name,
                    file_key=file_key,
                    filename=filename,
                    path=str(path),
                    count=count,
                    exists=exists,
                    length_stats=length_stats,
                    domain_counts=domain_counts,
                    phylum_counts=phylum_counts,
                )
            )

    def write_report(self, files: dict[str, str]) -> str:
        if not self.report.commands:
            self.load_manifest()

        self.record_steps(files)

        if self.report.commands:
            self.save_manifest()

        stem = self.report.build_stem or f"{self.report.dataset_name}_{self.report.primer_name}"
        report_name = f"{stem}_build_report.html"
        report_path = Path(self.working_dir) / report_name
        report_path.write_text(render_html(self.report), encoding="utf-8")
        self.report.report_path = str(report_path)
        return str(report_path)
