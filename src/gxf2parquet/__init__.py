"""GTF to Parquet conversion for genomics workflows."""

from .convert import detect_format, gff_to_parquet, gtf_to_parquet
from .query import query_gxf_parquet
from .read import read_gxf_parquet, read_source_format
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
    "query_gxf_parquet",
    "read_gxf_parquet",
    "read_source_format",
    "detect_format",
    "SchemaPreset",
    "BASE_PRESET",
    "GENCODE_PRESET",
    "ENSEMBL_PRESET",
    "get_preset",
]
