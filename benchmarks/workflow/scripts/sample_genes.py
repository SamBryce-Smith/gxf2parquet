#!/usr/bin/env python3
"""Draw a deterministic random sample of gene names (untimed staging step).

Reads the gene_name column from a Parquet build, takes the sorted set of
distinct non-null names and samples N of them with a fixed seed, so every
engine and replicate queries the same gene set. If the annotation has fewer
than N distinct names, all of them are written (and a warning printed).
"""

import argparse
import random
import sys
from pathlib import Path

import pyarrow.parquet as pq


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--parquet",
        required=True,
        type=Path,
        help="Parquet file or partitioned dataset directory",
    )
    p.add_argument("--n", required=True, type=int, help="number of gene names")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    column = pq.read_table(str(args.parquet), columns=["gene_name"]).column("gene_name")
    names = sorted({n for n in column.to_pylist() if n})
    if len(names) < args.n:
        print(
            f"WARNING: only {len(names)} distinct gene names (< {args.n}); using all",
            file=sys.stderr,
        )
    sample = random.Random(args.seed).sample(names, min(args.n, len(names)))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("".join(f"{n}\n" for n in sample))
    return 0


if __name__ == "__main__":
    sys.exit(main())
