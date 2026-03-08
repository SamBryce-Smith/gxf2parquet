"""GTF to Parquet conversion for genomics workflows."""

from .convert import detect_format, gff_to_parquet, gtf_to_parquet
from .query import query_gff_parquet
from .read import read_gtf_parquet
from .schema import (
    BASE_PRESET,
    ENSEMBL_PRESET,
    GENCODE_PRESET,
    SchemaPreset,
    get_preset,
)

__all__ = [
    "gtf_to_parquet",
    "gff_to_parquet",
    "query_gff_parquet",
    "read_gtf_parquet",
    "detect_format",
    "SchemaPreset",
    "BASE_PRESET",
    "GENCODE_PRESET",
    "ENSEMBL_PRESET",
    "get_preset",
]
