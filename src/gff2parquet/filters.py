"""Region and filter parsing utilities for GFF/GTF Parquet queries."""

from __future__ import annotations

_OP_MAP: dict[str, str] = {
    "eq": "==",
    "ne": "!=",
    "gt": ">",
    "lt": "<",
    "ge": ">=",
    "le": "<=",
    "isin": "in",
    "notin": "not in",
}

_STRAND_MAP: dict[str, str] = {
    "plus": "+",
    "minus": "-",
    "+": "+",
    "-": "-",
}


def parse_region(region: str) -> list[tuple]:
    """Parse a genomic region string into pyarrow filter tuples.

    Coordinates are 1-based closed (GTF/GFF convention), matching what is
    stored in Parquet.

    Args:
        region: Region string in one of the formats:
            - ``"chr1"`` — chromosome only
            - ``"chr1:1000-2000"`` — chromosome with 1-based closed interval

    Returns:
        List of filter tuples suitable for pyarrow predicate pushdown.

    Raises:
        ValueError: If the region string is malformed.

    Examples:
        >>> parse_region("chr1")
        [('Chromosome', '==', 'chr1')]
        >>> parse_region("chr1:1000-2000")
        [('Chromosome', '==', 'chr1'), ('Start', '>=', 1000), ('End', '<=', 2000)]
    """
    if ":" not in region:
        return [("Chromosome", "==", region)]

    chrom, interval = region.split(":", 1)

    if "-" not in interval:
        raise ValueError(
            f"Invalid region {region!r}. Expected 'CHROM' or 'CHROM:START-END'."
        )

    start_str, end_str = interval.split("-", 1)
    try:
        start = int(start_str)
        end = int(end_str)
    except ValueError:
        raise ValueError(
            f"Invalid coordinates in region {region!r}: start and end must be integers."
        )

    if start < 1:
        raise ValueError(
            f"Invalid region {region!r}: start must be >= 1 (1-based coordinates)."
        )
    if end < start:
        raise ValueError(f"Invalid region {region!r}: end must be >= start.")

    return [
        ("Chromosome", "==", chrom),
        ("Start", ">=", start),
        ("End", "<=", end),
    ]


def parse_strand(strand: str) -> str:
    """Convert a strand alias to the GTF/GFF strand symbol.

    Args:
        strand: One of ``'plus'``, ``'minus'``, ``'+'``, or ``'-'``.

    Returns:
        GTF strand symbol (``'+'`` or ``'-'``).

    Raises:
        ValueError: If the strand value is not recognised.
    """
    result = _STRAND_MAP.get(strand.lower())
    if result is None:
        raise ValueError(
            f"Unknown strand {strand!r}. Expected 'plus', 'minus', '+', or '-'."
        )
    return result


def parse_filter(col: str, op: str, val: str) -> tuple:
    """Parse a single column filter into a pyarrow filter tuple.

    Args:
        col: Column name (e.g., ``'gene_type'``).
        op: Operator alias (e.g., ``'eq'``, ``'ne'``, ``'isin'``).
            See :data:`_OP_MAP` for the full list.
        val: Value to compare against. For ``'isin'`` / ``'notin'``,
            provide a comma-separated list (e.g., ``'chr1,chr2'``).

    Returns:
        A pyarrow filter tuple ``(col, arrow_op, value)``.

    Raises:
        ValueError: If the operator is not recognised.
    """
    op_lower = op.lower()
    if op_lower not in _OP_MAP:
        valid = ", ".join(_OP_MAP.keys())
        raise ValueError(f"Unknown operator {op!r}. Valid operators: {valid}")

    arrow_op = _OP_MAP[op_lower]

    if op_lower in ("isin", "notin"):
        parsed_val: list[str] | str = [v.strip() for v in val.split(",")]
    else:
        parsed_val = val

    return (col, arrow_op, parsed_val)


def build_filters(
    *,
    regions: list[str] | None = None,
    strand: str | None = None,
    extra_filters: list[tuple] | None = None,
) -> list[list[tuple]] | list[tuple] | None:
    """Combine region, strand, and column filters into pyarrow DNF format.

    Multiple regions are OR-combined; strand and column filters are AND-combined
    into every region group.

    Args:
        regions: List of region strings (e.g., ``["chr1:1000-2000", "chr2"]``).
        strand: Strand alias or None to skip strand filtering.
        extra_filters: Additional ``(col, op, val)`` filter tuples to AND in.

    Returns:
        Filters in pyarrow DNF format, or ``None`` if nothing was specified.
    """
    strand_filter = [("Strand", "==", parse_strand(strand))] if strand else []
    col_filters = list(extra_filters) if extra_filters else []
    common = strand_filter + col_filters

    if regions:
        groups = [parse_region(r) + common for r in regions]
        return groups[0] if len(groups) == 1 else groups

    return common if common else None
