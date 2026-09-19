"""Parquet reading utilities."""

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import pyranges1 as pr

from .convert import METADATA_SOURCE_FORMAT


def read_source_format(parquet_path: str | Path) -> str | None:
    """Read the source format from Parquet file-level metadata.

    Returns the value of the ``gxf2parquet.source_format`` key (e.g. ``"gtf"``
    or ``"gff3"``), or ``None`` if the key is absent (e.g. files built before
    this metadata was added).

    Args:
        parquet_path: Path to a Parquet file or partitioned dataset directory.

    Returns:
        Source format string, or ``None`` if metadata key is not present.
    """
    meta = pq.read_metadata(str(parquet_path))
    raw = meta.metadata.get(METADATA_SOURCE_FORMAT)
    if raw is None:
        return None
    return raw.decode() if isinstance(raw, bytes) else raw


def read_gxf_parquet(
    parquet_path: str | Path,
    *,
    columns: list[str] | None = None,
    filters: list[tuple] | list[list[tuple]] | None = None,
    as_pyranges: bool = True,
) -> pr.PyRanges | pd.DataFrame:
    """Read a GXF (GTF/GFF) Parquet file into a PyRanges object or pandas DataFrame.

    Args:
        parquet_path: Path to Parquet file or partitioned dataset directory.
        columns: List of columns to read. If None, reads all columns.
        filters: Row group filters for predicate pushdown.
            Can be a list of tuples like [("Chromosome", "==", "chr1")]
            or a list of lists for OR conditions.
        as_pyranges: If True (default), returns a PyRanges object with 0-based coordinates.
            If False, returns a pandas DataFrame with 1-based coordinates.

    Returns:
        PyRanges object (if as_pyranges=True) or pandas DataFrame (if as_pyranges=False).

    Examples:
        # Read as PyRanges (default)
        gr = read_gxf_parquet("annotations.parquet")

        # Read as DataFrame with 1-based coordinates
        df = read_gxf_parquet("annotations.parquet", as_pyranges=False)

        # Read specific columns
        gr = read_gxf_parquet("annotations.parquet", columns=["Chromosome", "Start", "End", "gene_name"])

        # Read with filters (predicate pushdown)
        gr = read_gxf_parquet(
            "annotations.parquet",
            filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")]
        )
    """
    table = pq.read_table(
        str(parquet_path),
        columns=columns,
        filters=filters,
    )
    df = table.to_pandas()

    if as_pyranges:
        # Convert from 1-based GTF coordinates (stored in Parquet) to 0-based PyRanges coordinates
        # GTF/GFF are 1-based, but PyRanges uses 0-based coordinates internally
        df.loc[:, "Start"] = df.Start - 1
        return pr.PyRanges(df)

    return df
