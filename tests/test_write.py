"""Tests for write.py helpers."""

from pathlib import Path

import pytest

from gff2parquet.write import detect_output_format


class TestDetectOutputFormat:
    def test_gtf_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gtf") == "gtf"

    def test_gff3_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gff3") == "gff3"

    def test_gff_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gff") == "gff3"

    def test_parquet_extension(self, tmp_path):
        assert detect_output_format(tmp_path / "out.parquet") == "parquet"

    def test_gz_stripped(self, tmp_path):
        assert detect_output_format(tmp_path / "out.gtf.gz") == "gtf"

    def test_none_returns_default(self):
        assert detect_output_format(None) == "gtf"

    def test_none_custom_default(self):
        assert detect_output_format(None, default="gff3") == "gff3"

    def test_unknown_extension_returns_default(self, tmp_path):
        assert detect_output_format(tmp_path / "out.bed") == "gtf"
