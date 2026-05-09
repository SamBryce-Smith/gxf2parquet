"""Output writers for GFF/GTF and Parquet formats."""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pandas as pd
import pyranges1 as pr

# Columns required for a valid GTF/GFF3 file (the 8 fixed fields)
GTF_REQUIRED_COLUMNS: frozenset[str] = frozenset(
    {"Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame"}
)


def _check_columns_for_text_output(df: pd.DataFrame, fmt: str) -> None:
    """Raise ValueError if any core GTF/GFF3 columns are absent from *df*."""
    missing = GTF_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Cannot write {fmt.upper()} output: the following required columns are "
            f"missing from the data: {sorted(missing)}. "
            "Ensure all core GTF columns (Chromosome, Source, Feature, Start, End, "
            "Score, Strand, Frame) are included in --columns, or write to a Parquet "
            "file instead."
        )


def _warn_columns_for_parquet_output(df: pd.DataFrame) -> None:
    """Emit a warning if any core GTF/GFF3 columns are absent from *df*."""
    missing = GTF_REQUIRED_COLUMNS - set(df.columns)
    if missing:
        warnings.warn(
            f"Writing Parquet output without the following core GTF columns: "
            f"{sorted(missing)}. The resulting file will not be convertible back to a "
            "valid GTF/GFF3.",
            UserWarning,
            stacklevel=3,
        )


def _df_to_pyranges(df: pd.DataFrame) -> pr.PyRanges:
    """Convert a Parquet-sourced DataFrame (1-based coords) to a PyRanges object."""
    out = df.copy()
    out["Start"] = out["Start"] - 1  # 1-based → 0-based for PyRanges
    return pr.PyRanges(out)


def write_gtf(df: pd.DataFrame, output: Path | None = None) -> None:
    """Write a DataFrame with 1-based coordinates to GTF format.

    Args:
        df: DataFrame with 1-based Start/End coordinates (as stored in Parquet).
        output: Output file path, or ``None`` to write to stdout.

    Raises:
        ValueError: If any of the 8 required GTF columns are missing from *df*.
    """
    _check_columns_for_text_output(df, "gtf")
    gr = _df_to_pyranges(df)
    gr.to_gtf(sys.stdout if output is None else str(output))


_GFF3_HEADER = "##gff-version 3\n"


def write_gff3(df: pd.DataFrame, output: Path | None = None) -> None:
    """Write a DataFrame with 1-based coordinates to GFF3 format.

    Args:
        df: DataFrame with 1-based Start/End coordinates (as stored in Parquet).
        output: Output file path, or ``None`` to write to stdout.

    Raises:
        ValueError: If any of the 8 required GFF3 columns are missing from *df*.
    """
    _check_columns_for_text_output(df, "gff3")
    gr = _df_to_pyranges(df)
    content = _GFF3_HEADER + gr.to_gff3()
    if output is None:
        sys.stdout.write(content)
    else:
        output.write_text(content)


def write_parquet(
    df: pd.DataFrame,
    output: Path,
    *,
    compression: str = "zstd",
) -> None:
    """Write a DataFrame to a single Parquet file.

    Args:
        df: DataFrame to write.
        output: Output file path.
        compression: Parquet compression codec (``'zstd'``, ``'snappy'``,
            ``'gzip'``, or ``'none'``).
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    _warn_columns_for_parquet_output(df)
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(
        table,
        str(output),
        compression=compression if compression != "none" else None,
    )


def write_output(
    df: pd.DataFrame,
    output: Path | None,
    fmt: str,
    *,
    compression: str = "zstd",
) -> None:
    """Write *df* to *output* in the requested format.

    Args:
        df: DataFrame with 1-based Start/End coordinates (as stored in Parquet).
        output: Output file path, or ``None`` to write to stdout (text formats only).
        fmt: One of ``'gtf'``, ``'gff3'``, or ``'parquet'``.
        compression: Compression codec for Parquet output (default: ``'zstd'``).

    Raises:
        ValueError: If *fmt* is not recognised, or if Parquet output is requested
            without an output path.
    """
    if fmt == "gtf":
        write_gtf(df, output)
    elif fmt == "gff3":
        write_gff3(df, output)
    elif fmt == "parquet":
        if output is None:
            raise ValueError("An output path is required for Parquet format.")
        write_parquet(df, output, compression=compression)
    else:
        raise ValueError(f"Unknown output format {fmt!r}. Expected 'gtf', 'gff3', or 'parquet'.")


def detect_output_format(output: Path | None, *, default: str = "gtf") -> str:
    """Infer the output format from a file extension.

    Args:
        output: Output path, or ``None`` for stdout.
        default: Fallback format when ``output`` is ``None`` or the extension is
            not recognised. Must be one of ``'gtf'``, ``'gff3'``, ``'parquet'``.

    Returns:
        One of ``'gtf'``, ``'gff3'``, or ``'parquet'``.
    """
    if output is None:
        return default

    name = output.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]

    if name.endswith(".gtf"):
        return "gtf"
    if name.endswith(".gff3") or name.endswith(".gff"):
        return "gff3"
    if name.endswith(".parquet"):
        return "parquet"

    return default
