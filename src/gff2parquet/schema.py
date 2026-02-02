"""Schema presets for GTF to Parquet conversion."""

from dataclasses import dataclass, field


@dataclass
class SchemaPreset:
    """Schema configuration for GTF conversion.

    Attributes:
        categorical_columns: Columns to convert to pandas Categorical dtype.
        list_columns: Columns that contain multiple values (stored as list<string>).
    """

    categorical_columns: list[str] = field(default_factory=list)
    list_columns: list[str] = field(default_factory=list)


GENCODE_PRESET = SchemaPreset(
    categorical_columns=[
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
        "gene_type",
        "transcript_type",
    ],
    list_columns=["tag", "ont"],
)

ENSEMBL_PRESET = SchemaPreset(
    categorical_columns=[
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
        "gene_biotype",
        "transcript_biotype",
    ],
    list_columns=["tag"],
)

_PRESETS = {
    "gencode": GENCODE_PRESET,
    "ensembl": ENSEMBL_PRESET,
}


def get_preset(name: str) -> SchemaPreset:
    """Get a schema preset by name.

    Args:
        name: Preset name ('gencode' or 'ensembl').

    Returns:
        The requested SchemaPreset.

    Raises:
        ValueError: If preset name is not recognized.
    """
    name_lower = name.lower()
    if name_lower not in _PRESETS:
        available = ", ".join(_PRESETS.keys())
        raise ValueError(f"Unknown preset '{name}'. Available presets: {available}")
    return _PRESETS[name_lower]
