"""GTF/GFF to Parquet conversion."""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyranges1 as pr

from .schema import BASE_PRESET, SchemaPreset


def _convert_to_parquet(
    df: pd.DataFrame,
    parquet_path: str | Path,
    *,
    preset: SchemaPreset,
    partition_cols: list[str] | None = None,
    compression: str = "zstd",
) -> None:
    """Internal conversion function shared by gtf_to_parquet and gff_to_parquet.

    Args:
        df: DataFrame with pyranges coordinates (0-based, half-open).
        parquet_path: Path for output Parquet file or directory (if partitioned).
        preset: Schema preset defining categorical and list columns.
        partition_cols: Columns to partition by (e.g., ["Chromosome", "Feature"]).
        compression: Compression codec ('zstd', 'snappy', 'gzip', 'none').
    """
    # Convert pyranges coordinates (0-based, half-open) to GTF/GFF coordinates (1-based, closed)
    # pyranges Start is 0-based, GTF/GFF is 1-based
    df["Start"] = df["Start"] + 1
    # pyranges End is already correct (half-open end equals closed end in 1-based)

    # Apply categorical dtypes to specified columns
    for col in preset.categorical_columns:
        if col in df.columns:
            df[col] = df[col].astype("category")

    # Handle list columns: pyranges returns duplicate attributes as comma-separated
    # strings (e.g. "basic,Ensembl_canonical"). Split these into actual Python lists
    # so they are stored as Arrow list<string> in Parquet.
    for col in preset.list_columns:
        if col in df.columns:
            df[col] = df[col].apply(
                lambda x: x.split(",") if isinstance(x, str) else None
            )

    # Convert to Arrow Table
    table = pa.Table.from_pandas(df, preserve_index=False)

    # Explicitly cast list columns to list<string> to guarantee the Arrow type
    # regardless of whether the column was all-null or had mixed content.
    for col in preset.list_columns:
        if col in table.schema.names:
            col_idx = table.schema.get_field_index(col)
            list_array = pa.array(df[col].tolist(), type=pa.list_(pa.string()))
            table = table.set_column(col_idx, col, list_array)

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


def gtf_to_parquet(
    gtf_path: str | Path,
    parquet_path: str | Path,
    *,
    preset: SchemaPreset = BASE_PRESET,
    partition_cols: list[str] | None = None,
    compression: str = "zstd",
) -> None:
    """Convert a GTF file to Parquet format.

    Args:
        gtf_path: Path to input GTF file (may be gzipped).
        parquet_path: Path for output Parquet file or directory (if partitioned).
        preset: Schema preset defining categorical and list columns.
            Defaults to BASE_PRESET.
        partition_cols: Columns to partition by (e.g., ["Chromosome", "Feature"]).
            If provided, output will be a directory with Hive-style partitioning.
        compression: Compression codec ('zstd', 'snappy', 'gzip', 'none').
    """
    # Parse GTF using pyranges with duplicate_attr=True so that repeated attribute
    # keys (e.g. tag "basic"; tag "Ensembl_canonical") are joined into a single
    # comma-separated string rather than silently dropping all but the last value.
    gr = pr.read_gtf(str(gtf_path), duplicate_attr=True)
    df = pd.DataFrame(gr)

    _convert_to_parquet(
        df,
        parquet_path,
        preset=preset,
        partition_cols=partition_cols,
        compression=compression,
    )


def gff_to_parquet(
    gff_path: str | Path,
    parquet_path: str | Path,
    *,
    preset: SchemaPreset = BASE_PRESET,
    partition_cols: list[str] | None = None,
    compression: str = "zstd",
) -> None:
    """Convert a GFF3 file to Parquet format.

    Args:
        gff_path: Path to input GFF3 file (may be gzipped).
        parquet_path: Path for output Parquet file or directory (if partitioned).
        preset: Schema preset defining categorical and list columns.
            Defaults to BASE_PRESET.
        partition_cols: Columns to partition by (e.g., ["Chromosome", "Feature"]).
            If provided, output will be a directory with Hive-style partitioning.
        compression: Compression codec ('zstd', 'snappy', 'gzip', 'none').
    """
    # Parse GFF3 using pyranges (pyranges1 returns a DataFrame subclass)
    gr = pr.read_gff3(str(gff_path))
    df = pd.DataFrame(gr)

    _convert_to_parquet(
        df,
        parquet_path,
        preset=preset,
        partition_cols=partition_cols,
        compression=compression,
    )
