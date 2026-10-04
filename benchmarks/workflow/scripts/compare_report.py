#!/usr/bin/env python3
"""Write a before-vs-after Markdown report from two pipeline runs.

Joins the `bench_build_summary.tsv`, `bench_query_summary.tsv` and `disk.tsv`
of a "before" and an "after" results directory (see compare.sh) on their cell
keys and reports median metrics, the change, and whether it stands out from
replicate noise.

A change is called "better"/"worse" only when the before and after min-max
ranges do not overlap; otherwise it is "noise". With fewer than 3 replicates
no call is made. The naive `pyranges1` engine runs identical code on both
sides, so its rows are a control: they should sit at about 0%.

Standard library only, so it runs anywhere.
"""

import argparse
import csv
import json
import sys
from pathlib import Path

NAIVE_ENGINE = "pyranges1"
MB = 1024 * 1024

BUILD_KEYS = ("annotation", "build")
QUERY_KEYS = ("annotation", "engine", "query", "param")

# (column stem, label, scale, format) for metrics with median/min/max columns.
BUILD_METRICS = [
    ("wall_s", "wall time (s)", 1, "{:.3f}"),
    ("max_rss_gb", "peak RSS (MB)", 1024, "{:.1f}"),
]
QUERY_TIME_METRICS = [
    ("wall_s", "wall time (s)", 1, "{:.3f}"),
    ("op_wall_s", "read+filter (s)", 1, "{:.3f}"),
]
QUERY_MEMORY_METRICS = [
    ("max_rss_gb", "peak RSS (MB)", 1024, "{:.1f}"),
]


def read_tsv(path: Path) -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def num(row: dict | None, col: str) -> float | None:
    if row is None:
        return None
    value = row.get(col, "")
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def index(rows: list[dict], keys: tuple[str, ...]) -> dict[tuple, dict]:
    return {tuple(r.get(k, "") for k in keys): r for r in rows}


def pct(before: float | None, after: float | None) -> float | None:
    if before is None or after is None or before == 0:
        return None
    return (after - before) / before * 100


def fmt_pct(value: float | None) -> str:
    return "n/a" if value is None else f"{value:+.1f}%"


def verdict(before: dict, after: dict, stem: str, min_reps: int) -> str:
    """better/worse only if the replicate min-max ranges don't overlap."""
    reps = min(num(before, "rep_count") or 0, num(after, "rep_count") or 0)
    if reps < min_reps:
        return f"n/a (<{min_reps} reps)"
    b_lo, b_hi = num(before, f"{stem}_min"), num(before, f"{stem}_max")
    a_lo, a_hi = num(after, f"{stem}_min"), num(after, f"{stem}_max")
    if None in (b_lo, b_hi, a_lo, a_hi):
        return "n/a"
    if a_hi < b_lo:
        return "**better**"
    if a_lo > b_hi:
        return "**worse**"
    return "noise"


def metric_cells(
    before, after, stem, scale, fmt, min_reps
) -> tuple[list[str], float | None]:
    b = num(before, f"{stem}_median")
    a = num(after, f"{stem}_median")
    change = pct(b, a)
    cells = [
        "n/a" if b is None else fmt.format(b * scale),
        "n/a" if a is None else fmt.format(a * scale),
        fmt_pct(change),
        verdict(before, after, stem, min_reps),
    ]
    return cells, change


def table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def comparison_table(before_rows, after_rows, keys, metrics, min_reps, extra=None):
    """One row per cell present on both sides; returns (markdown, changes)."""
    before_idx, after_idx = index(before_rows, keys), index(after_rows, keys)
    header = list(keys)
    for _, label, _, _ in metrics:
        header += [f"{label} before", "after", "Δ", "call"]
    if extra:
        header += [label for label, _ in extra]
    rows, changes = [], []
    for key in sorted(set(before_idx) & set(after_idx)):
        b, a = before_idx[key], after_idx[key]
        row = list(key)
        for stem, _, scale, fmt in metrics:
            cells, change = metric_cells(b, a, stem, scale, fmt, min_reps)
            row += cells
            changes.append((key, stem, change, cells[-1]))
        if extra:
            row += [fn(b, a) for _, fn in extra]
        rows.append(row)
    missing = sorted(set(before_idx) ^ set(after_idx))
    return table(header, rows), changes, missing


def mb_change(col: str):
    def fn(before: dict, after: dict) -> str:
        b, a = num(before, col), num(after, col)
        if b is None or a is None:
            return "n/a"
        return f"{b / MB:.1f} → {a / MB:.1f} ({fmt_pct(pct(b, a))})"

    return fn


def correctness(before_q, after_q, before_disk, after_disk) -> list[str]:
    problems = []
    b_idx, a_idx = index(before_q, QUERY_KEYS), index(after_q, QUERY_KEYS)
    for key in sorted(set(b_idx) & set(a_idx)):
        for col in ("n_rows", "n_cols"):
            if b_idx[key].get(col) != a_idx[key].get(col):
                problems.append(
                    f"`{'/'.join(key)}`: {col} {b_idx[key].get(col)} → {a_idx[key].get(col)}"
                )
    for side, rows in (("before", before_q), ("after", after_q)):
        for r in rows:
            if r.get("any_failed", "").lower() == "true":
                problems.append(
                    f"{side}: `{r['annotation']}/{r['engine']}/{r['query']}/{r['param']}` had a failed run"
                )
            if r.get("n_rows_consistent", "true").lower() == "false":
                problems.append(
                    f"{side}: engines disagree on n_rows for `{r['annotation']}/{r['query']}/{r['param']}`"
                )
    bd = index(before_disk, ("annotation", "artifact"))
    ad = index(after_disk, ("annotation", "artifact"))
    for key in sorted(set(bd) & set(ad)):
        if bd[key].get("bytes") != ad[key].get("bytes"):
            problems.append(
                f"`{'/'.join(key)}` size on disk {bd[key].get('bytes')} → {ad[key].get('bytes')} bytes"
            )
    return problems


def source_check(rows: list[dict], expected_src: str, side: str) -> list[str]:
    """gxf2parquet engines must have imported the intended source tree."""
    if not expected_src:
        return []
    want = str(Path(expected_src).resolve() / "gxf2parquet")
    bad = {
        r.get("gxf2parquet_path")
        for r in rows
        if r.get("engine", "").startswith("gxf2parquet")
        and r.get("gxf2parquet_path") != want
    }
    return [
        f"{side}: gxf2parquet imported from `{p}`, expected `{want}`"
        for p in sorted(bad - {None})
    ]


def summarise(changes, engines_filter) -> list[str]:
    lines = []
    for stem, label in (
        ("wall_s", "wall time"),
        ("max_rss_gb", "peak RSS"),
        ("op_wall_s", "read+filter time"),
    ):
        sel = [c for c in changes if c[1] == stem and engines_filter(c[0])]
        if not sel:
            continue
        better = sum("better" in c[3] for c in sel)
        worse = sum("worse" in c[3] for c in sel)
        vals = sorted(c[2] for c in sel if c[2] is not None)
        median = vals[len(vals) // 2] if vals else None
        lines.append(
            f"- **{label}**: {better} better, {worse} worse, "
            f"{len(sel) - better - worse} noise/n.a. of {len(sel)} cells; "
            f"median change {fmt_pct(median)}"
        )
    return lines


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--before", required=True, type=Path, help="results dir of the 'before' run"
    )
    p.add_argument(
        "--after", required=True, type=Path, help="results dir of the 'after' run"
    )
    p.add_argument("--before-label", default="before")
    p.add_argument("--after-label", default="after")
    p.add_argument("--command", default="", help="command line that produced the runs")
    p.add_argument(
        "--min-reps",
        type=int,
        default=3,
        help="replicates needed before calling a change (default 3)",
    )
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    load = lambda d, f: read_tsv(d / f)  # noqa: E731
    b_build, a_build = (
        load(args.before, "bench_build_summary.tsv"),
        load(args.after, "bench_build_summary.tsv"),
    )
    b_query, a_query = (
        load(args.before, "bench_query_summary.tsv"),
        load(args.after, "bench_query_summary.tsv"),
    )
    b_disk, a_disk = load(args.before, "disk.tsv"), load(args.after, "disk.tsv")
    b_cfg, a_cfg = (
        read_json(args.before / "run_config.json"),
        read_json(args.after / "run_config.json"),
    )
    host = read_json(args.after / "host.json")
    if not (b_build or b_query) or not (a_build or a_query):
        print("ERROR: missing summary TSVs in --before or --after", file=sys.stderr)
        return 1

    first = (a_query or a_build)[0]
    commit = lambda rows: rows[0].get("gxf2parquet_commit", "") if rows else ""  # noqa: E731
    out = [
        "# gxf2parquet: before vs after",
        "",
        "## Setup",
        "",
        table(
            ["", "before", "after"],
            [
                ["source", args.before_label, args.after_label],
                [
                    "commit",
                    commit(b_query or b_build) or "n/a",
                    commit(a_query or a_build) or "n/a",
                ],
                [
                    "replicates",
                    str(b_cfg.get("replicates", "?")),
                    str(a_cfg.get("replicates", "?")),
                ],
                [
                    "cache mode",
                    str(b_cfg.get("cache_mode", "?")),
                    str(a_cfg.get("cache_mode", "?")),
                ],
                [
                    "annotations",
                    ", ".join(b_cfg.get("annotations", [])),
                    ", ".join(a_cfg.get("annotations", [])),
                ],
            ],
        ),
        "",
        f"Host: {host.get('cpu_model', '?')}, {host.get('cpu_count_logical', '?')} logical CPUs, "
        f"{(host.get('total_mem_bytes') or 0) / 1024**3:.1f} GB RAM, {host.get('os', '?')} {host.get('os_release', '')}. "
        f"pandas {first.get('pandas_version', '?')}, pyarrow {first.get('pyarrow_version', '?')}, "
        f"pyranges1 {first.get('pyranges1_version', '?')}, Python {first.get('python_version', '?')}. "
        "Both sides ran in the same pixi environment; only the gxf2parquet source differs.",
        "",
    ]
    if args.command:
        out += ["Command:", "", "```bash", args.command, "```", ""]

    build_md, build_changes, build_missing = comparison_table(
        b_build,
        a_build,
        BUILD_KEYS,
        BUILD_METRICS,
        args.min_reps,
        extra=[("output size (MB)", mb_change("output_bytes"))],
    )
    query_time_md, time_changes, query_missing = comparison_table(
        b_query, a_query, QUERY_KEYS, QUERY_TIME_METRICS, args.min_reps
    )
    query_mem_md, mem_changes, _ = comparison_table(
        b_query,
        a_query,
        QUERY_KEYS,
        QUERY_MEMORY_METRICS,
        args.min_reps,
        extra=[("result in memory (MB)", mb_change("df_mem_bytes_median"))],
    )
    query_changes = time_changes + mem_changes

    problems = correctness(b_query, a_query, b_disk, a_disk)
    problems += source_check(b_query, b_cfg.get("gxf2parquet_src", ""), "before")
    problems += source_check(a_query, a_cfg.get("gxf2parquet_src", ""), "after")

    out += ["## Summary", ""]
    out += (
        ["Build (`gxf2parquet build`):"]
        + summarise(build_changes, lambda k: True)
        + [""]
    )
    out += ["Queries, gxf2parquet engines:"]
    out += summarise(query_changes, lambda k: k[1] != NAIVE_ENGINE) + [""]
    control = [
        abs(c[2])
        for c in query_changes
        if c[0][1] == NAIVE_ENGINE and c[2] is not None and c[1] == "wall_s"
    ]
    if control:
        out += [
            f"Control (naive `{NAIVE_ENGINE}`, same code on both sides): wall-time changes up to "
            f"{max(control):.1f}% (median {sorted(control)[len(control) // 2]:.1f}%). Changes of this size "
            "or smaller elsewhere are within run-to-run variation.",
            "",
        ]
    out += ["Correctness:", ""]
    out += (
        [f"- ⚠️ {msg}" for msg in problems]
        if problems
        else [
            "- Identical `n_rows`/`n_cols` for every query, identical Parquet sizes, no failed runs."
        ]
    )
    out += [""]

    out += ["## Build", "", build_md, ""]
    out += ["## Queries: time", "", query_time_md, ""]
    out += ["## Queries: memory", "", query_mem_md, ""]
    if build_missing or query_missing:
        out += [
            "Cells present on only one side (not compared): "
            + ", ".join("/".join(k) for k in build_missing + query_missing),
            "",
        ]
    out += [
        "## How to read this",
        "",
        "- Medians over replicates; Δ = (after − before) / before. Negative is better for every metric here.",
        f"- **call**: better/worse only when the before and after min–max ranges don't overlap "
        f"(needs ≥ {args.min_reps} replicates); otherwise noise.",
        "- *wall time* and *peak RSS* come from GNU `time -v` for the whole process (start-up, imports, "
        "work). *read+filter* is timed inside the process around the query itself.",
        "- Not measured here: the `gxf2parquet query` CLI output paths (streamed gtf/gff3/bed writing, "
        "no copies for tsv/csv/parquet output). The pipeline calls the Python API.",
        "",
    ]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(out))
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
