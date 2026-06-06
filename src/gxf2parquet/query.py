"""Query utilities for GFF/GTF Parquet files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import pyranges1 as pr

from .filters import build_filters


def query_gxf_parquet(
    parquet_path: str | Path,
    *,
    regions: list[str] | None = None,
    strand: str | list[str] | None = None,
    filters: list[tuple[str, str, Any]] | None = None,
    columns: list[str] | None = None,
    as_pyranges: bool = True,
) -> pr.PyRanges | pd.DataFrame:
    """Query a GXF (GTF/GFF) Parquet file with optional region, strand, and column filters.

    Args:
        parquet_path: Path to a Parquet file or partitioned dataset directory.
        regions: Genomic regions to filter by (e.g., ``["chr1", "chr2:1000-5000"]``).
            Multiple regions are OR-combined. Coordinates are 1-based closed (GTF/GFF).
        strand: Strand filter. A single string (``'+'``, ``'-'``, ``'plus'``,
            ``'minus'``) is broadcast to all regions. A list must match the length
            of ``regions`` and is paired positionally (one strand per region).
            ``None`` applies no strand filtering.
        filters: List of pyarrow-compatible filter tuples
            ``(column, operator, value)`` where *operator* uses pyarrow syntax
            (``'=='``, ``'!='``, ``'>'``, ``'<'``, ``'>='``, ``'<='``, ``'in'``,
            ``'not in'``). Multiple filters are AND-combined.
        columns: Restrict output to these columns. ``None`` reads all columns.
        as_pyranges: If ``True`` (default), return a PyRanges object with 0-based
            coordinates. If ``False``, return a DataFrame with 1-based coordinates.

    Returns:
        PyRanges (0-based) or DataFrame (1-based).
    """
    parquet_path = Path(parquet_path)

    combined = build_filters(regions=regions, strand=strand, extra_filters=filters)
    table = pq.read_table(
        str(parquet_path),
        columns=columns,
        filters=combined,
    )
    df = table.to_pandas()

    if as_pyranges:
        out = df.copy()
        out["Start"] = out["Start"] - 1
        return pr.PyRanges(out)

    return df
