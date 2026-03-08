"""Command-line interface for GTF to Parquet conversion."""

import argparse
import sys
from pathlib import Path

from .convert import detect_format, gff_to_parquet, gtf_to_parquet
from .filters import parse_filter
from .query import query_gff_parquet
from .schema import get_preset
from .write import detect_output_format, write_gff3, write_gtf, write_parquet


def _cmd_build(args: argparse.Namespace) -> int:
    """Handler for the `build` subcommand."""
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    try:
        fmt = detect_format(args.input)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    preset = get_preset(args.preset)

    try:
        if fmt == "gtf":
            gtf_to_parquet(
                args.input,
                args.output,
                preset=preset,
                partition_cols=args.partition_cols,
                compression=args.compression,
            )
        else:
            gff_to_parquet(
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
    if not args.input.exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        return 1

    # Parse --filter options (each is [COL, OP, VAL])
    extra_filters = None
    if args.filter:
        try:
            extra_filters = [parse_filter(col, op, val) for col, op, val in args.filter]
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    try:
        df = query_gff_parquet(
            args.input,
            regions=args.region or None,
            strand=args.strand,
            filters=extra_filters,
            columns=args.columns or None,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine output format
    output = args.output
    if args.format:
        fmt = args.format
    else:
        fmt = detect_output_format(output, default="gtf")

    try:
        if fmt == "gtf":
            write_gtf(df, output)
        elif fmt == "gff3":
            write_gff3(df, output)
        elif fmt == "parquet":
            if output is None:
                print(
                    "Error: --output is required for parquet format.", file=sys.stderr
                )
                return 1
            write_parquet(df, output, compression=args.compression)
        else:
            print(f"Error: Unknown output format {fmt!r}.", file=sys.stderr)
            return 1
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    return 0


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
        description=(
            "Convert a GTF/GFF3 annotation file to Parquet format. "
            "The input format is auto-detected from the file extension "
            "(.gtf, .gff, .gff3, with optional .gz)."
        ),
    )
    build_parser.add_argument(
        "input",
        type=Path,
        help="Input GTF or GFF3 file (may be gzipped)",
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
        choices=["base", "gencode", "ensembl"],
        default="gencode",
        help="Schema preset for categorical columns (default: gencode)",
    )

    # --- query subcommand ---
    query_parser = subparsers.add_parser(
        "query",
        help="Query a Parquet annotation file.",
        description=(
            "Query a GFF/GTF Parquet file by region, strand, or column filters "
            "and write the results to a file or stdout."
        ),
    )
    query_parser.add_argument(
        "input",
        type=Path,
        help="Input Parquet file or directory",
    )
    query_parser.add_argument(
        "--region",
        action="append",
        metavar="REGION",
        help=(
            "Genomic region in 'CHROM' or 'CHROM:START-END' format (1-based, closed). "
            "May be repeated; multiple regions are OR-combined."
        ),
    )
    query_parser.add_argument(
        "--strand",
        choices=["plus", "minus", "+", "-"],
        metavar="STRAND",
        help="Strand filter: 'plus' (+) or 'minus' (-)",
    )
    query_parser.add_argument(
        "--filter",
        nargs=3,
        action="append",
        metavar=("COL", "OP", "VAL"),
        help=(
            "Column filter, e.g. --filter gene_type eq protein_coding. "
            "Operators: eq, ne, gt, lt, ge, le, isin, notin. "
            "For isin/notin, VAL is comma-separated. "
            "May be repeated; all filters are AND-combined."
        ),
    )
    query_parser.add_argument(
        "--columns",
        nargs="+",
        metavar="COL",
        help="Columns to include in the output (reads all by default)",
    )
    query_parser.add_argument(
        "--output",
        type=Path,
        metavar="FILE",
        help="Output file path. Omit to write to stdout.",
    )
    query_parser.add_argument(
        "--format",
        choices=["gtf", "gff3", "parquet"],
        metavar="FORMAT",
        help=(
            "Output format: gtf, gff3, or parquet. "
            "Auto-detected from --output extension when omitted (default: gtf)."
        ),
    )
    query_parser.add_argument(
        "--compression",
        choices=["zstd", "snappy", "gzip", "none"],
        default="zstd",
        help="Compression codec for parquet output (default: zstd)",
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
