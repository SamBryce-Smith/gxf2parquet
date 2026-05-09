"""Tests for source-format metadata stored in Parquet files."""

import tempfile
from pathlib import Path

import pyarrow.parquet as pq
import pytest
import pyranges1 as pr

from gff2parquet import gff_to_parquet, gtf_to_parquet, read_source_format
from gff2parquet.convert import METADATA_SOURCE_FORMAT


@pytest.fixture
def ensembl_gtf_path():
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = Path(tmpdir) / "test.gtf"
        pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
        yield gtf_path


@pytest.fixture
def temp_parquet_path():
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir) / "output.parquet"


class TestSourceFormatMetadata:
    def test_gtf_stores_source_format(self, ensembl_gtf_path, temp_parquet_path):
        gtf_to_parquet(ensembl_gtf_path, temp_parquet_path)
        meta = pq.read_metadata(str(temp_parquet_path))
        assert METADATA_SOURCE_FORMAT in meta.metadata
        assert meta.metadata[METADATA_SOURCE_FORMAT] == b"gtf"

    def test_gff_stores_source_format(self, tmp_path):
        """GFF3 conversion stores 'gff3' as source_format."""
        gff_path = tmp_path / "test.gff3"
        parquet_path = tmp_path / "output.parquet"

        # Write a minimal valid GFF3 file
        gff_path.write_text(
            "##gff-version 3\n"
            "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\t"
            "ID=gene1;gene_id=ENSG00000223972;gene_name=DDX11L1\n"
        )

        gff_to_parquet(gff_path, parquet_path)
        meta = pq.read_metadata(str(parquet_path))
        assert METADATA_SOURCE_FORMAT in meta.metadata
        assert meta.metadata[METADATA_SOURCE_FORMAT] == b"gff3"


class TestReadSourceFormat:
    def test_gtf_source_format(self, ensembl_gtf_path, temp_parquet_path):
        gtf_to_parquet(ensembl_gtf_path, temp_parquet_path)
        assert read_source_format(temp_parquet_path) == "gtf"

    def test_gff3_source_format(self, tmp_path):
        gff_path = tmp_path / "test.gff3"
        parquet_path = tmp_path / "output.parquet"
        gff_path.write_text(
            "##gff-version 3\n"
            "chr1\tHAVANA\tgene\t11869\t14409\t.\t+\t.\t"
            "ID=gene1;gene_id=ENSG00000223972;gene_name=DDX11L1\n"
        )
        gff_to_parquet(gff_path, parquet_path)
        assert read_source_format(parquet_path) == "gff3"

    def test_missing_metadata_returns_none(self, tmp_path):
        import pyarrow as pa
        import pyarrow.parquet as pq

        parquet_path = tmp_path / "no_meta.parquet"
        table = pa.table({"a": [1, 2, 3]})
        pq.write_table(table, str(parquet_path))
        assert read_source_format(parquet_path) is None
