"""Command-line interface for GTF to Parquet conversion."""

import argparse
import sys
from pathlib import Path

from .convert import gtf_to_parquet
from .schema import get_preset


def _cmd_build(args: argparse.Namespace) -> int:
    """Handler for the `build` subcommand."""
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


def _cmd_query(args: argparse.Namespace) -> int:
    """Handler for the `query` subcommand."""
    print("query: not yet implemented", file=sys.stderr)
    return 1


def main(argv: list[str] | None = None) -> int:
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        prog="gff2parquet",
        description="Convert and query GFF/GTF annotation files in Parquet format.",
    )

    subparsers = parser.add_subparsers(dest="subcommand")

    # --- build subcommand ---
    build_parser = subparsers.add_parser(
        "build",
        help="Convert a GTF/GFF3 file to Parquet format.",
        description="Convert a GTF/GFF3 annotation file to Parquet format.",
    )
    build_parser.add_argument(
        "input",
        type=Path,
        help="Input GTF file (may be gzipped)",
    )
    build_parser.add_argument(
        "output",
        type=Path,
        help="Output Parquet file or directory (if partitioned)",
    )
    build_parser.add_argument(
        "--partition-cols",
        nargs="+",
        metavar="COL",
        help="Columns to partition by (e.g., Chromosome Feature)",
    )
    build_parser.add_argument(
        "--compression",
        choices=["zstd", "snappy", "gzip", "none"],
        default="zstd",
        help="Compression codec (default: zstd)",
    )
    build_parser.add_argument(
        "--preset",
        choices=["gencode", "ensembl"],
        default="gencode",
        help="Schema preset for categorical columns (default: gencode)",
    )

    # --- query subcommand (placeholder) ---
    subparsers.add_parser(
        "query",
        help="Query a Parquet annotation file.",
        description="Query a GFF/GTF Parquet file by region, strand, or column filters.",
    )

    args = parser.parse_args(argv)

    if args.subcommand is None:
        parser.print_help()
        return 0

    if args.subcommand == "build":
        return _cmd_build(args)
    if args.subcommand == "query":
        return _cmd_query(args)

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
