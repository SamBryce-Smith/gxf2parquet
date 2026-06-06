"""Output format detection for GFF/GTF and Parquet formats."""

from __future__ import annotations

from pathlib import Path


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
