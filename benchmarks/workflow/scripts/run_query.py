#!/usr/bin/env python3
"""Timed query driver: one query, one engine, one process.

The pipeline wraps this script in GNU `time -v`, so wall time and peak RSS
cover the whole process (interpreter start-up, imports, read, filter). The
JSON written to --out-json adds finer-grained numbers:

  import_s      time to import the engine's libraries
  op_wall_s     time for the read + filter itself
  n_rows/n_cols shape of the result (cross-engine sanity check)
  df_mem_bytes  deep in-memory size of the result

Both engines end with an in-memory pyranges1 PyRanges:

  pyranges1    pr.read_gtf(..., duplicate_attr=True) on the GTF, then filter
               in memory (the naive workflow)
  gxf2parquet  read_gxf_parquet(...) with column selection / predicate
               pushdown on a Parquet build
"""

import argparse
import json
import sys
import time
from pathlib import Path

# Mirrors gxf2parquet.write.CORE_GXF_COLUMNS. Defined here rather than
# imported so the naive engine doesn't pay for importing gxf2parquet (and
# pyarrow); the gxf2parquet engine asserts the two stay in sync.
CORE_GXF_COLUMNS = (
    "Chromosome",
    "Source",
    "Feature",
    "Start",
    "End",
    "Score",
    "Strand",
    "Frame",
)


def read_gene_names(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text().splitlines() if line.strip()]


def run_pyranges1(args, gene_names):
    t0 = time.perf_counter()
    import pyranges1 as pr

    import_s = time.perf_counter() - t0

    t0 = time.perf_counter()
    gr = pr.read_gtf(str(args.input), duplicate_attr=True)
    if args.query == "chrom":
        gr = gr[gr["Chromosome"] == args.chrom]
    elif args.query == "columns":
        gr = gr[[c for c in (*CORE_GXF_COLUMNS, *args.columns) if c in gr.columns]]
    elif args.query == "gene_set":
        gr = gr[gr["gene_name"].isin(gene_names)]
    op_wall_s = time.perf_counter() - t0
    return gr, import_s, op_wall_s


def run_gxf2parquet(args, gene_names):
    t0 = time.perf_counter()
    from gxf2parquet import read_gxf_parquet
    from gxf2parquet.write import CORE_GXF_COLUMNS as PKG_CORE

    import_s = time.perf_counter() - t0
    assert tuple(PKG_CORE) == CORE_GXF_COLUMNS, "CORE_GXF_COLUMNS out of sync"

    columns = None
    filters = None
    if args.query == "chrom":
        filters = [("Chromosome", "==", args.chrom)]
    elif args.query == "columns":
        columns = [*CORE_GXF_COLUMNS, *args.columns]
    elif args.query == "gene_set":
        filters = [("gene_name", "in", gene_names)]

    t0 = time.perf_counter()
    gr = read_gxf_parquet(args.input, columns=columns, filters=filters)
    op_wall_s = time.perf_counter() - t0
    return gr, import_s, op_wall_s


def main() -> int:
    p = argparse.ArgumentParser(description="Run one benchmark query.")
    p.add_argument("--engine", required=True, choices=["pyranges1", "gxf2parquet"])
    p.add_argument(
        "--input",
        required=True,
        type=Path,
        help="GTF (pyranges1) or Parquet file/dataset dir (gxf2parquet)",
    )
    p.add_argument(
        "--query", required=True, choices=["full", "chrom", "columns", "gene_set"]
    )
    p.add_argument("--chrom", help="chromosome for --query chrom")
    p.add_argument(
        "--columns",
        nargs="+",
        default=[],
        help="attribute columns for --query columns (core columns are always kept)",
    )
    p.add_argument(
        "--gene-names-file",
        type=Path,
        help="one gene_name per line, for --query gene_set",
    )
    p.add_argument("--out-json", required=True, type=Path)
    args = p.parse_args()

    if args.query == "chrom" and not args.chrom:
        p.error("--query chrom requires --chrom")
    if args.query == "columns" and not args.columns:
        p.error("--query columns requires --columns")
    if args.query == "gene_set" and not args.gene_names_file:
        p.error("--query gene_set requires --gene-names-file")

    gene_names = read_gene_names(args.gene_names_file) if args.gene_names_file else []

    runner = run_pyranges1 if args.engine == "pyranges1" else run_gxf2parquet
    gr, import_s, op_wall_s = runner(args, gene_names)

    result = {
        "engine": args.engine,
        "query": args.query,
        "import_s": round(import_s, 6),
        "op_wall_s": round(op_wall_s, 6),
        "n_rows": int(len(gr)),
        "n_cols": int(len(gr.columns)),
        "df_mem_bytes": int(gr.memory_usage(deep=True).sum()),
        "result_type": type(gr).__name__,
        "n_gene_names": len(gene_names),
    }
    if args.engine == "gxf2parquet":
        # Already imported (untimed here): records which source tree was used,
        # so before/after comparisons can confirm the PYTHONPATH switch worked.
        import gxf2parquet

        result["gxf2parquet_path"] = str(Path(gxf2parquet.__file__).resolve().parent)
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
