"""Tests for write.py helpers."""

from gxf2parquet.write import detect_output_format


class TestDetectOutputFormat:
    def test_gtf_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gtf") == "gtf"

    def test_gff3_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gff3") == "gff3"

    def test_gff_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gff") == "gff3"

    def test_bed_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.bed") == "bed"

    def test_tsv_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.tsv") == "tsv"

    def test_csv_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.csv") == "csv"

    def test_parquet_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.parquet") == "parquet"

    def test_gz_stripped(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gtf.gz") == "gtf"

    def test_none_returns_none(self):
        assert detect_output_format(None) is None

    def test_unknown_extension_returns_none(self, tmp_path):
        assert detect_output_format(tmp_path / "out.xyz") is None
