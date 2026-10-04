#!/usr/bin/env python3
"""Merge per-run GNU time output, query results and disk sizes into TSVs.

Inputs are the exact file lists Snakemake hands over (not a directory glob),
so results left over from runs with a different config are never mixed in.

Outputs:
  bench_build.tsv  one row per (annotation, build, rep)
  bench_query.tsv  one row per (annotation, engine, query, param, rep)
  disk.tsv         one row per (annotation, artifact)

Adapted from fulcrumgenomics/riker benchmark-pipeline (MIT License,
Copyright (c) 2026 Fulcrum Genomics LLC).
"""

import argparse
import csv
import json
import platform
import re
import sys
from importlib import metadata
from pathlib import Path

# parse_gnu_time sits next to this script.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from parse_gnu_time import parse as parse_time  # noqa: E402

# <results_dir>/build/{annotation}/{build}/rep{rep}/time.txt
BUILD_RE = re.compile(
    r"/build/(?P<annotation>[^/]+)/(?P<build>[^/]+)/rep(?P<rep>\d+)/time\.txt$"
)
# <results_dir>/query/{annotation}/{engine}/{query}/{param}/rep{rep}/time.txt
QUERY_RE = re.compile(
    r"/query/(?P<annotation>[^/]+)/(?P<engine>[^/]+)/(?P<query>[^/]+)/"
    r"(?P<param>[^/]+)/rep(?P<rep>\d+)/time\.txt$"
)

GXF2PARQUET_PREFIX = "gxf2parquet-"


def parse_path(regex: re.Pattern, path: str) -> dict:
    m = regex.search(path)
    if not m:
        raise ValueError(f"path doesn't match expected schema: {path}")
    d = m.groupdict()
    d["rep"] = int(d["rep"])
    return d


def timing_columns(time_txt: str) -> dict:
    t = parse_time(time_txt)
    max_rss_kb = int(t.get("max_rss_kb") or 0)
    return {
        "wall_s": float(t.get("wall_s") or 0.0),
        "user_s": float(t.get("user_s") or 0.0),
        "sys_s": float(t.get("sys_s") or 0.0),
        "cpu_percent": float(t.get("cpu_percent") or 0.0),
        "max_rss_kb": max_rss_kb,
        "max_rss_gb": round(max_rss_kb / (1024 * 1024), 4),
        "exit_status": int(t.get("exit_status") or 0),
    }


def package_versions() -> dict:
    out = {"python_version": platform.python_version()}
    for pkg in ("gxf2parquet", "pyranges1", "pyarrow", "pandas"):
        try:
            out[f"{pkg}_version"] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            out[f"{pkg}_version"] = ""
    return out


def host_columns(host_json: str) -> dict:
    with open(host_json) as fh:
        host = json.load(fh)
    keys = (
        "hostname",
        "arch",
        "cpu_model",
        "cpu_count_logical",
        "total_mem_bytes",
        "os",
        "os_release",
    )
    return {f"host_{k}": host.get(k) for k in keys}


def write_tsv(path: str, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    """Write rows as a TSV; an empty file when there are no rows."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        print(f"WARNING: no rows for {path}", file=sys.stderr)
    with open(path, "w", newline="") as fh:
        if not rows and fieldnames is None:
            return
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames or list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--host", required=True, help="host.json")
    p.add_argument(
        "--run-config",
        required=True,
        help="run_config.json (the resolved Snakemake config)",
    )
    p.add_argument("--build-times", nargs="*", default=[])
    p.add_argument("--query-times", nargs="*", default=[])
    p.add_argument("--disk", nargs="*", default=[], help="disk_usage.py JSON files")
    p.add_argument("--out-build", required=True)
    p.add_argument("--out-query", required=True)
    p.add_argument("--out-disk", required=True)
    args = p.parse_args()

    with open(args.run_config) as fh:
        config = json.load(fh)
    builds = config["builds"]
    cache_mode = config["cache_mode"]
    naive_input = config.get("naive_input", "gtf")
    common = {"cache_mode": cache_mode, **package_versions(), **host_columns(args.host)}

    disk_rows = [json.loads(Path(f).read_text()) for f in args.disk]
    disk_bytes = {(r["annotation"], r["artifact"]): r["bytes"] for r in disk_rows}
    disk_files = {(r["annotation"], r["artifact"]): r["n_files"] for r in disk_rows}

    build_rows = []
    for time_txt in args.build_times:
        meta = parse_path(BUILD_RE, time_txt)
        cfg = builds[meta["build"]]
        key = (meta["annotation"], meta["build"])
        build_rows.append(
            {
                **meta,
                "partition_cols": "+".join(cfg.get("partition_cols") or []) or "none",
                "compression": cfg.get("compression", "zstd"),
                **timing_columns(time_txt),
                "output_bytes": disk_bytes.get(key),
                "output_n_files": disk_files.get(key),
                **common,
            }
        )

    query_rows = []
    for time_txt in args.query_times:
        meta = parse_path(QUERY_RE, time_txt)
        engine = meta["engine"]
        if engine.startswith(GXF2PARQUET_PREFIX):
            family, build = "gxf2parquet", engine[len(GXF2PARQUET_PREFIX) :]
            artifact = build
        else:
            family, build = engine, ""
            artifact = naive_input
        result_json = Path(time_txt).with_name("result.json")
        result = json.loads(result_json.read_text()) if result_json.exists() else {}
        query_rows.append(
            {
                **meta,
                "engine_family": family,
                "build": build,
                **timing_columns(time_txt),
                "import_s": result.get("import_s"),
                "op_wall_s": result.get("op_wall_s"),
                "n_rows": result.get("n_rows"),
                "n_cols": result.get("n_cols"),
                "df_mem_bytes": result.get("df_mem_bytes"),
                "input_artifact": artifact,
                "input_bytes": disk_bytes.get((meta["annotation"], artifact)),
                **common,
            }
        )

    write_tsv(args.out_build, build_rows)
    write_tsv(args.out_query, query_rows)
    write_tsv(args.out_disk, disk_rows, ["annotation", "artifact", "bytes", "n_files"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
