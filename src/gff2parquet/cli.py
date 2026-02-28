"""Command-line interface for GTF to Parquet conversion."""

import argparse
import sys
from pathlib import Path

from .convert import gtf_to_parquet
from .schema import get_preset


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        prog="gff2parquet",
        description="Convert GTF annotation files to Parquet format.",
    )

    parser.add_argument(
        "input",
        type=Path,
        help="Input GTF file (may be gzipped)",
    )
    parser.add_argument(
        "output",
        type=Path,
        help="Output Parquet file or directory (if partitioned)",
    )
    parser.add_argument(
        "--partition-cols",
        nargs="+",
        metavar="COL",
        help="Columns to partition by (e.g., Chromosome Feature)",
    )
    parser.add_argument(
        "--compression",
        choices=["zstd", "snappy", "gzip", "none"],
        default="zstd",
        help="Compression codec (default: zstd)",
    )
    parser.add_argument(
        "--preset",
        choices=["gencode", "ensembl"],
        default="gencode",
        help="Schema preset for categorical columns (default: gencode)",
    )

    args = parser.parse_args(argv)

    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    preset = get_preset(args.preset)

    try:
        gtf_to_parquet(
            args.input,
            args.output,
            preset=preset,
            partition_cols=args.partition_cols,
            compression=args.compression,
        )
        print(f"Converted {args.input} -> {args.output}")
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
