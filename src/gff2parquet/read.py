"""Parquet reading utilities."""

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


def read_gtf_parquet(
    parquet_path: str | Path,
    *,
    columns: list[str] | None = None,
    filters: list[tuple] | list[list[tuple]] | None = None,
) -> pd.DataFrame:
    """Read a GTF Parquet file into a pandas DataFrame.

    Args:
        parquet_path: Path to Parquet file or partitioned dataset directory.
        columns: List of columns to read. If None, reads all columns.
        filters: Row group filters for predicate pushdown.
            Can be a list of tuples like [("Chromosome", "==", "chr1")]
            or a list of lists for OR conditions.

    Returns:
        pandas DataFrame with the requested data.

    Examples:
        # Read all data
        df = read_gtf_parquet("annotations.parquet")

        # Read specific columns
        df = read_gtf_parquet("annotations.parquet", columns=["Chromosome", "Start", "End", "gene_name"])

        # Read with filters (predicate pushdown)
        df = read_gtf_parquet(
            "annotations.parquet",
            filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")]
        )
    """
    table = pq.read_table(
        str(parquet_path),
        columns=columns,
        filters=filters,
    )
    return table.to_pandas()
