#!/usr/bin/env python3
"""Make a larger GTF by repeating a small one under renamed chromosomes.

Copy 0 is the input unchanged. Copy i >= 1 renames every chromosome to
``<chrom>_c<i>`` and appends ``_c<i>`` to the ID/name attributes, so genes,
transcripts and exons stay distinct across copies. The result keeps the
input's structure and column mix at N times the size: enough rows for real
memory and runtime differences to show above interpreter start-up, without
downloading a full annotation.

A single-chromosome query on the original name then selects 1/N of the rows,
and a gene-name sample spans all copies.
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

# Attributes made unique per copy (GTF `key "value";` form).
SUFFIXED_ATTRIBUTES = (
    "gene_id",
    "transcript_id",
    "exon_id",
    "gene_name",
    "transcript_name",
    "protein_id",
)
_ATTR_RE = re.compile(r"\b(" + "|".join(SUFFIXED_ATTRIBUTES) + r') "([^"]*)"')


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="", compresslevel=6)
    return open(path, mode, newline="")


def scale_line(line: str, copy: int) -> str:
    if copy == 0 or line.startswith("#"):
        return line
    fields = line.split("\t")
    if len(fields) < 9:
        return line
    suffix = f"_c{copy}"
    fields[0] = fields[0] + suffix
    fields[8] = _ATTR_RE.sub(
        lambda m: f'{m.group(1)} "{m.group(2)}{suffix}"', fields[8]
    )
    return "\t".join(fields)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input", required=True, type=Path, help="GTF (optionally .gz)")
    p.add_argument("--copies", required=True, type=int, help="total number of copies")
    p.add_argument("--out", required=True, type=Path, help="output GTF (.gz to gzip)")
    args = p.parse_args()
    if args.copies < 1:
        p.error("--copies must be >= 1")

    with open_text(args.input, "r") as fh:
        lines = fh.readlines()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.out, "w") as out:
        for copy in range(args.copies):
            # Header comments only once, from copy 0.
            out.writelines(
                scale_line(line, copy)
                for line in lines
                if copy == 0 or not line.startswith("#")
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
