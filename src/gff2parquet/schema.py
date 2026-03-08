"""Schema presets for GTF to Parquet conversion."""

from dataclasses import dataclass, field


@dataclass
class SchemaPreset:
    """Schema configuration for GTF conversion.

    Attributes:
        categorical_columns: Columns to convert to pandas Categorical dtype.
        list_columns: Columns that contain multiple values (stored as list<string>).
        column_dtypes: Explicit dtype overrides for individual columns, expressed as
            pandas dtype strings (e.g. ``{"exon_number": "Int16"}``). Applied after
            categorical conversion. Use pandas nullable integer types (capital-I
            ``"Int8"``, ``"Int16"``, ``"Int32"``) for integer-valued string columns
            that may contain NULLs.
    """

    categorical_columns: list[str] = field(default_factory=list)
    list_columns: list[str] = field(default_factory=list)
    column_dtypes: dict[str, str] = field(default_factory=dict)


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
        "Score",                    # always "." in GENCODE
        "Frame",                    # 4 values: ".", "0", "1", "2"
        # GENCODE attribute columns (low-cardinality, used for filtering)
        "gene_type",
        "transcript_type",
        "level",                    # annotation confidence: "1", "2", "3"
        "transcript_support_level", # TSL: "1"–"5", "NA"
    ],
    list_columns=["tag", "ont"],
    column_dtypes={
        # exon_number is an ordinal integer; nullable because gene/transcript rows
        # carry no exon_number value.
        "exon_number": "Int16",
    },
)

ENSEMBL_PRESET = SchemaPreset(
    categorical_columns=[
        # Standard GTF columns
        "Chromosome",
        "Source",
        "Feature",
        "Strand",
        "Score",                    # always "." in Ensembl releases
        "Frame",                    # 4 values: ".", "0", "1", "2"
        # Ensembl attribute columns (low-cardinality, used for filtering)
        "gene_biotype",
        "transcript_biotype",
        "gene_source",              # e.g. "havana", "ensembl_havana", "ensembl"
        "transcript_source",        # same small vocabulary
        "transcript_support_level", # TSL: "1"–"5", "NA"
    ],
    list_columns=["tag"],
    column_dtypes={
        # Integer-valued string fields; nullable because only some feature rows carry them.
        "exon_number":        "Int16",
        "gene_version":       "Int16",
        "transcript_version": "Int16",
        "exon_version":       "Int16",
    },
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
