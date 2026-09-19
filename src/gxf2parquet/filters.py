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


def _coerce(value: str) -> int | float | str:
    """Coerce a string to int, then float, falling back to str."""
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def parse_filter(tokens: list[str]) -> tuple:
    """Parse a single column filter token list into a pyarrow filter tuple.

    The token list must be ``[column, op, value, ...]`` where extra tokens
    after the operator are collected as the value set for ``isin``/``notin``.
    Numeric strings are coerced to ``int`` or ``float`` where possible.

    Args:
        tokens: List of tokens, e.g. ``["Feature", "eq", "exon"]`` or
            ``["gene_type", "isin", "protein_coding", "lncRNA"]``.

    Returns:
        A pyarrow filter tuple ``(col, arrow_op, value)``.

    Raises:
        ValueError: If the token list has fewer than 3 elements, the operator
            is not recognised, or ``isin``/``notin`` have no values.
    """
    if len(tokens) < 3:
        raise ValueError(
            f"Filter requires at least 3 tokens (column op value ...), got: {tokens!r}"
        )

    col, op, *rest = tokens
    op_lower = op.lower()
    if op_lower not in _OP_MAP:
        valid = ", ".join(_OP_MAP.keys())
        raise ValueError(f"Unknown operator {op!r}. Valid operators: {valid}")

    arrow_op = _OP_MAP[op_lower]

    if op_lower in ("isin", "notin"):
        if not rest:
            raise ValueError(f"Operator {op!r} requires at least one value.")
        parsed_val: list | int | float | str = [_coerce(v) for v in rest]
    else:
        parsed_val = _coerce(rest[0])

    return (col, arrow_op, parsed_val)


def build_filters(
    *,
    regions: list[str] | None = None,
    strand: str | list[str] | None = None,
    extra_filters: list[tuple] | None = None,
) -> list[list[tuple]] | list[tuple] | None:
    """Combine region, strand, and column filters into pyarrow DNF format.

    Multiple regions are OR-combined; strand and column filters are AND-combined
    into every region group.

    When ``strand`` is a single string it is broadcast to every region.
    When ``strand`` is a list it must have the same length as ``regions`` and
    each value is paired positionally with its corresponding region.

    Args:
        regions: List of region strings (e.g., ``["chr1:1000-2000", "chr2"]``).
        strand: Strand alias(es). A single string broadcasts to all regions;
            a list is paired positionally with ``regions``.  ``None`` skips
            strand filtering entirely.
        extra_filters: Additional ``(col, op, val)`` filter tuples to AND in.

    Returns:
        Filters in pyarrow DNF format, or ``None`` if nothing was specified.

    Raises:
        ValueError: If ``strand`` is a list whose length differs from ``regions``.
    """
    col_filters = list(extra_filters) if extra_filters else []

    if regions:
        # Resolve per-region strand filters
        if isinstance(strand, list):
            if len(strand) != len(regions):
                raise ValueError(
                    f"Number of strand values ({len(strand)}) must match "
                    f"number of regions ({len(regions)}) when providing "
                    f"per-region strands."
                )
            strand_filters_per_region = [
                [("Strand", "==", parse_strand(s))] for s in strand
            ]
        elif strand:
            strand_filters_per_region = [
                [("Strand", "==", parse_strand(strand))] for _ in regions
            ]
        else:
            strand_filters_per_region = [[] for _ in regions]

        groups = [
            parse_region(r) + sf + col_filters
            for r, sf in zip(regions, strand_filters_per_region)
        ]
        return groups[0] if len(groups) == 1 else groups

    # No regions — apply strand and column filters globally
    if strand:
        if isinstance(strand, list):
            strand_filter = [("Strand", "==", parse_strand(s)) for s in strand]
        else:
            strand_filter = [("Strand", "==", parse_strand(strand))]
    else:
        strand_filter = []

    common = strand_filter + col_filters
    return common if common else None
