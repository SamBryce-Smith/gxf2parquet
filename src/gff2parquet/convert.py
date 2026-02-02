"""GTF to Parquet conversion."""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyranges1 as pr

from .schema import GENCODE_PRESET, SchemaPreset


def gtf_to_parquet(
    gtf_path: str | Path,
    parquet_path: str | Path,
    *,
    preset: SchemaPreset | None = None,
    partition_cols: list[str] | None = None,
    compression: str = "zstd",
) -> None:
    """Convert a GTF file to Parquet format.

    Args:
        gtf_path: Path to input GTF file (may be gzipped).
        parquet_path: Path for output Parquet file or directory (if partitioned).
        preset: Schema preset defining categorical and list columns.
            Defaults to GENCODE_PRESET.
        partition_cols: Columns to partition by (e.g., ["Chromosome", "Feature"]).
            If provided, output will be a directory with Hive-style partitioning.
        compression: Compression codec ('zstd', 'snappy', 'gzip', 'none').
    """
    if preset is None:
        preset = GENCODE_PRESET

    # Parse GTF using pyranges (pyranges1 returns a DataFrame subclass)
    gr = pr.read_gtf(str(gtf_path))
    df = pd.DataFrame(gr)

    # Convert pyranges coordinates (0-based, half-open) to GTF coordinates (1-based, closed)
    # pyranges Start is 0-based, GTF is 1-based
    df["Start"] = df["Start"] + 1
    # pyranges End is already correct (half-open end equals closed end in 1-based)

    # Apply categorical dtypes to specified columns
    for col in preset.categorical_columns:
        if col in df.columns:
            df[col] = df[col].astype("category")

    # Handle list columns - these are already parsed as lists by pyranges
    # Just ensure they're properly typed for Arrow
    # (pyranges stores multi-value attributes as comma-separated strings or lists)

    # Convert to Arrow Table
    table = pa.Table.from_pandas(df, preserve_index=False)

    # Write to Parquet
    if partition_cols:
        # Validate partition columns exist
        missing = [c for c in partition_cols if c not in df.columns]
        if missing:
            raise ValueError(f"Partition columns not found: {missing}")

        pq.write_to_dataset(
            table,
            root_path=str(parquet_path),
            partition_cols=partition_cols,
            compression=compression if compression != "none" else None,
        )
    else:
        pq.write_table(
            table,
            str(parquet_path),
            compression=compression if compression != "none" else None,
        )
