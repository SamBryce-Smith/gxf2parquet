"""Command-line interface for GTF to Parquet conversion."""

import argparse
import sys
from pathlib import Path

from .convert import detect_format, gff_to_parquet, gtf_to_parquet
from .filters import parse_filter, parse_strand
from .query import query_gxf_parquet
from .read import read_source_format
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

    # Flatten --region lists (action="append" + nargs="+" gives [[r1], [r2, r3], ...])
    regions: list[str] | None = None
    if args.region:
        regions = [r for group in args.region for r in group]

    # Flatten --strand lists and map plus/minus -> +/-
    strands: list[str] = []
    if args.strand:
        raw_strands = [s for group in args.strand for s in group]
        try:
            strands = [parse_strand(s) for s in raw_strands]
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    # Validate strand count: 0 (none), 1 (broadcast), or len(regions) (per-region)
    n_regions = len(regions) if regions else 0
    n_strands = len(strands)
    if n_strands > 1 and n_strands != n_regions:
        print(
            f"Error: {n_strands} strand value(s) given but {n_regions} region(s) "
            f"specified. Provide 0, 1, or exactly {n_regions} strand value(s).",
            file=sys.stderr,
        )
        return 1

    strand_arg: str | list[str] | None
    if n_strands == 0:
        strand_arg = None
    elif n_strands == 1:
        strand_arg = strands[0]
    else:
        strand_arg = strands

    # Parse --filter options (each is a token list [COL, OP, VAL, ...])
    extra_filters = None
    if args.filter:
        try:
            extra_filters = [parse_filter(tokens) for tokens in args.filter]
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    try:
        df = query_gxf_parquet(
            args.input,
            regions=regions,
            strand=strand_arg,
            filters=extra_filters,
            columns=args.columns or None,
            as_pyranges=False,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    # Determine output format: explicit --format > extension > Parquet metadata > "gtf"
    output = args.output
    if args.format:
        fmt = args.format
    elif output is not None and detect_output_format(output) != "gtf":
        fmt = detect_output_format(output)
    else:
        fmt = read_source_format(args.input) or "gtf"
        if output is not None:
            detected = detect_output_format(output)
            if detected == "parquet":
                fmt = "parquet"

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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # GENCODE GTF (default preset)
  gff2parquet build gencode.v47.annotation.gtf.gz gencode.parquet

  # Partition by Chromosome and Feature for faster region/feature queries
  gff2parquet build gencode.v47.annotation.gtf.gz gencode.parquet \\
      --partition-cols Chromosome Feature

  # Ensembl GTF with matching preset
  gff2parquet build Homo_sapiens.GRCh38.113.gtf.gz ensembl.parquet --preset ensembl

  # GFF3 input (auto-detected), snappy compression
  gff2parquet build Homo_sapiens.GRCh38.113.gff3.gz ensembl.parquet \\
      --preset ensembl --compression snappy
""",
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
            "Query a GXF (GTF/GFF) Parquet file by region, strand, or column filters "
            "and write the results to a file or stdout."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
examples:
  # All genes on chr1, written as GTF to stdout
  gff2parquet query gencode.parquet --region chr1 --filter Feature eq gene

  # Region query with strand, select columns, write to GTF file
  gff2parquet query gencode.parquet \\
      --region chr1:11869-14409 --strand plus \\
      --columns Chromosome Start End Strand Feature gene_name \\
      --output chr1_region.gtf

  # Per-region strand: chr1 plus-strand, chr2 minus-strand
  gff2parquet query gencode.parquet \\
      --region chr1 --region chr2 \\
      --strand plus --strand minus

  # Filter by feature set and gene type (space-separated values for isin)
  gff2parquet query gencode.parquet \\
      --filter Feature isin exon CDS \\
      --filter gene_type eq protein_coding \\
      --output coding_exons.parquet

  # Whole chromosome, GFF3 output
  gff2parquet query gencode.parquet --region chr2 --output chr2.gff3
""",
    )
    query_parser.add_argument(
        "input",
        type=Path,
        help="Input Parquet file or directory",
    )
    query_parser.add_argument(
        "--region",
        nargs="+",
        action="append",
        metavar="REGION",
        help=(
            "Genomic region in 'CHROM' or 'CHROM:START-END' format (1-based, closed). "
            "May be repeated; multiple regions are OR-combined."
        ),
    )
    query_parser.add_argument(
        "--strand",
        nargs="+",
        action="append",
        metavar="STRAND",
        help=(
            "Strand filter: 'plus' (+) or 'minus' (-). "
            "A single value is broadcast to all regions. "
            "Multiple values must match the number of --region flags "
            "and are paired positionally."
        ),
    )
    query_parser.add_argument(
        "--filter",
        nargs="+",
        action="append",
        metavar="TOKEN",
        help=(
            "Column filter as 'COL OP VAL [VAL ...]', e.g. --filter Feature eq exon "
            "or --filter gene_type isin protein_coding lncRNA. "
            "Operators: eq, ne, gt, lt, ge, le, isin, notin. "
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
            "Auto-detected from --output extension when omitted; "
            "defaults to the source format stored in Parquet metadata (fallback: gtf)."
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
