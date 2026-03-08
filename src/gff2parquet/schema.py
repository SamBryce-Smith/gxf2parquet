"""Schema presets for GTF to Parquet conversion."""

from dataclasses import dataclass, field


@dataclass
class SchemaPreset:
    """Schema configuration for GTF conversion.

    Attributes:
        categorical_columns: Columns to convert to pandas Categorical dtype.
        list_columns: Columns that contain multiple values (stored as list<string>).
        int16_columns: Columns to cast to pandas nullable Int16 dtype.
        int32_columns: Columns to cast to pandas nullable Int32 dtype.
        int64_columns: Columns to cast to pandas nullable Int64 dtype.
    """

    categorical_columns: list[str] = field(default_factory=list)
    list_columns: list[str] = field(default_factory=list)
    int16_columns: list[str] = field(default_factory=list)
    int32_columns: list[str] = field(default_factory=list)
    int64_columns: list[str] = field(default_factory=list)


BASE_PRESET = SchemaPreset(
    categorical_columns=[
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
    ],
    list_columns=[],
)

GENCODE_PRESET = SchemaPreset(
    categorical_columns=[
        # Standard GTF columns
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
        "Score",   # always "." in GENCODE; 1-element dict is near-zero overhead
        "Frame",   # 4 values: ".", "0", "1", "2"
        # GENCODE attribute columns (low-cardinality, used for filtering)
        "gene_type",
        "transcript_type",
        "level",                    # annotation confidence: "1", "2", "3"
        "transcript_support_level", # TSL: "1"–"5", "NA"
    ],
    list_columns=["tag", "ont"],
    int16_columns=[
        # exon_number is an ordinal integer (max ~363 in human GENCODE); nullable
        # because gene/transcript rows carry no exon_number value.
        "exon_number",
    ],
)

ENSEMBL_PRESET = SchemaPreset(
    categorical_columns=[
        # Standard GTF columns
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
        "Score",   # always "." in Ensembl releases
        "Frame",   # 4 values: ".", "0", "1", "2"
        # Ensembl attribute columns (low-cardinality, used for filtering)
        "gene_biotype",
        "transcript_biotype",
        "gene_source",              # e.g. "havana", "ensembl_havana", "ensembl"
        "transcript_source",        # same small vocabulary
        "transcript_support_level", # TSL: "1"–"5", "NA"
    ],
    list_columns=["tag"],
    int16_columns=[
        # Integer-valued string fields; nullable because only some feature rows carry them.
        "exon_number",
        "gene_version",
        "transcript_version",
        "exon_version",
    ],
)

_PRESETS = {
    "base": BASE_PRESET,
    "gencode": GENCODE_PRESET,
    "ensembl": ENSEMBL_PRESET,
}


def get_preset(name: str) -> SchemaPreset:
    """Get a schema preset by name.

    Args:
        name: Preset name ('base', 'gencode', or 'ensembl').

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
