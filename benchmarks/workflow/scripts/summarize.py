#!/usr/bin/env python3
"""Collapse replicates to per-cell median + min + max.

For query results, also flags cells where engines disagree on the number of
returned rows (`n_rows_consistent`), a cheap check that every engine answers
the same question.

Adapted from fulcrumgenomics/riker benchmark-pipeline (MIT License,
Copyright (c) 2026 Fulcrum Genomics LLC).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

GROUP_KEYS = {
    "build": ["annotation", "build", "partition_cols", "compression", "cache_mode"],
    "query": [
        "annotation",
        "engine",
        "engine_family",
        "build",
        "query",
        "param",
        "input_artifact",
        "cache_mode",
    ],
}

NUMERIC_COLS = {
    "build": ["wall_s", "user_s", "sys_s", "cpu_percent", "max_rss_kb", "max_rss_gb"],
    "query": [
        "wall_s",
        "user_s",
        "sys_s",
        "cpu_percent",
        "max_rss_kb",
        "max_rss_gb",
        "import_s",
        "op_wall_s",
        "df_mem_bytes",
    ],
}

# Constant within a group; carried through from the first replicate.
CARRY_COLS = {
    "build": ["output_bytes", "output_n_files"],
    "query": ["n_rows", "n_cols", "input_bytes"],
}
CARRY_COMMON = [
    "gxf2parquet_version",
    "pyranges1_version",
    "pyarrow_version",
    "pandas_version",
    "python_version",
    "host_hostname",
    "host_arch",
    "host_cpu_model",
    "host_cpu_count_logical",
    "host_total_mem_bytes",
    "host_os",
    "host_os_release",
]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--kind", required=True, choices=["build", "query"])
    p.add_argument("--in", dest="inp", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    try:
        df = pd.read_csv(args.inp, sep="\t")
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if len(df) == 0:
        Path(args.out).write_text("")
        return 0

    keys = GROUP_KEYS[args.kind]
    df[keys] = df[keys].fillna("")
    grouped = df.groupby(keys, dropna=False)

    agg = {
        c: ["median", "min", "max"] for c in NUMERIC_COLS[args.kind] if c in df.columns
    }
    out = grouped.agg(agg)
    out.columns = [f"{col}_{stat}" for col, stat in out.columns]
    carry = [c for c in CARRY_COLS[args.kind] + CARRY_COMMON if c in df.columns]
    out = out.join(grouped[carry].first())
    out["rep_count"] = grouped.size()
    out["any_failed"] = grouped["exit_status"].apply(lambda s: bool((s != 0).any()))
    out = out.reset_index()

    if args.kind == "query":
        cell = ["annotation", "query", "param"]
        n_distinct = out.groupby(cell)["n_rows"].transform("nunique")
        out["n_rows_consistent"] = n_distinct == 1
        for _, bad in out[~out["n_rows_consistent"]].groupby(cell):
            detail = ", ".join(f"{e}={n}" for e, n in zip(bad["engine"], bad["n_rows"]))
            print(
                f"WARNING: engines disagree on n_rows for "
                f"{tuple(bad[cell].iloc[0])}: {detail}",
                file=sys.stderr,
            )

    out.to_csv(args.out, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
