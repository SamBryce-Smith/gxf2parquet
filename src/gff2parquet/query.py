"""Query utilities for GFF/GTF Parquet files."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from .filters import build_filters


def query_gff_parquet(
    parquet_path: str | Path,
    *,
    regions: list[str] | None = None,
    strand: str | None = None,
    filters: list[tuple] | list[list[tuple]] | None = None,
    columns: list[str] | None = None,
) -> pd.DataFrame:
    """Query a GFF/GTF Parquet file and return a DataFrame.

    Coordinates in the returned DataFrame are 1-based closed (GTF/GFF native),
    exactly as stored in the Parquet file.

    Args:
        parquet_path: Path to a Parquet file or partitioned dataset directory.
        regions: Genomic regions to filter by (e.g., ``["chr1", "chr2:1000-5000"]``).
            Multiple regions are OR-combined. Coordinates are 1-based closed.
        strand: Strand filter. One of ``'plus'``, ``'minus'``, ``'+'``, or ``'-'``.
        filters: Additional pyarrow filter tuples to AND with region/strand filters.
            Can be a flat list like ``[("gene_type", "==", "protein_coding")]``
            or a nested list-of-lists for OR groups.
        columns: Columns to read. ``None`` reads all columns.

    Returns:
        pandas DataFrame with 1-based Start/End coordinates.

    Examples:
        >>> df = query_gff_parquet("annotations.parquet", regions=["chr1"])
        >>> df = query_gff_parquet(
        ...     "annotations.parquet",
        ...     regions=["chr1:1-248956422"],
        ...     strand="plus",
        ...     filters=[("Feature", "==", "gene")],
        ...     columns=["Chromosome", "Start", "End", "Strand", "gene_name"],
        ... )
    """
    combined = build_filters(regions=regions, strand=strand, extra_filters=filters)
    table = pq.read_table(
        str(parquet_path),
        columns=columns,
        filters=combined,
    )
    return table.to_pandas()
