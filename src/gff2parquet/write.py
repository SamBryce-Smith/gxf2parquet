"""Output writers for GFF/GTF and Parquet formats."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pandas as pd
import pyranges1 as pr


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
    """
    gr = _df_to_pyranges(df)
    if output is None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / "out.gtf"
            gr.to_gtf(str(tmp_path))
            print(tmp_path.read_text(), end="")
    else:
        gr.to_gtf(str(output))


_GFF3_HEADER = "##gff-version 3\n"


def write_gff3(df: pd.DataFrame, output: Path | None = None) -> None:
    """Write a DataFrame with 1-based coordinates to GFF3 format.

    Args:
        df: DataFrame with 1-based Start/End coordinates (as stored in Parquet).
        output: Output file path, or ``None`` to write to stdout.
    """
    gr = _df_to_pyranges(df)
    if output is None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / "out.gff3"
            gr.to_gff3(str(tmp_path))
            print(_GFF3_HEADER + tmp_path.read_text(), end="")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir) / "out.gff3"
            gr.to_gff3(str(tmp_path))
            output.write_text(_GFF3_HEADER + tmp_path.read_text())


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
