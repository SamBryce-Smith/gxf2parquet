"""Output format detection for GFF/GTF and Parquet formats."""

from __future__ import annotations

from pathlib import Path

# The eight core GTF/GFF columns. Always emitted for gtf/gff3 output regardless of
# any ``--columns`` selection; also used to warn when writing Parquet without them.
CORE_GXF_COLUMNS = (
    "Chromosome",
    "Source",
    "Feature",
    "Start",
    "End",
    "Score",
    "Strand",
    "Frame",
)

# Columns preserved for BED output. ``Name`` is absent in GXF data and is filled with
# ``"."`` by ``PyRanges.to_bed``; ``Score`` and ``Strand`` carry real values worth keeping.
CORE_BED_COLUMNS = (
    "Chromosome",
    "Start",
    "End",
    "Score",
    "Strand",
)


def detect_output_format(output: Path | None) -> str | None:
    """Infer the output format from a file extension.

    Args:
        output: Output path, or ``None`` for stdout.

    Returns:
        One of ``'gtf'``, ``'gff3'``, ``'bed'``, ``'tsv'``, ``'csv'``, or
        ``'parquet'``, or ``None`` when ``output`` is ``None`` or the extension is
        not recognised (the caller applies its own fallback).
    """
    if output is None:
        return None

    name = output.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]

    if name.endswith(".gtf"):
        return "gtf"
    if name.endswith(".gff3") or name.endswith(".gff"):
        return "gff3"
    if name.endswith(".bed"):
        return "bed"
    if name.endswith(".tsv"):
        return "tsv"
    if name.endswith(".csv"):
        return "csv"
    if name.endswith(".parquet"):
        return "parquet"

    return None
