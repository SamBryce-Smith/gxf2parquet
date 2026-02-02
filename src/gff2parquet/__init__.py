"""GTF to Parquet conversion for genomics workflows."""

from .convert import gtf_to_parquet
from .read import read_gtf_parquet
from .schema import ENSEMBL_PRESET, GENCODE_PRESET, SchemaPreset, get_preset

__all__ = [
    "gtf_to_parquet",
    "read_gtf_parquet",
    "SchemaPreset",
    "GENCODE_PRESET",
    "ENSEMBL_PRESET",
    "get_preset",
]
