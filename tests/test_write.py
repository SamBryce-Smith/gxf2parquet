"""Round-trip and output tests for write.py helpers."""

import sys
from io import StringIO

import pandas as pd
import pytest
import pyarrow.parquet as pq
import pyranges1 as pr

from gff2parquet import gtf_to_parquet, query_gxf_parquet, write_gff3, write_gtf
from gff2parquet.schema import ENSEMBL_PRESET
from gff2parquet.write import write_output, write_parquet


@pytest.fixture
def ensembl_gtf_path(tmp_path):
    gtf_path = tmp_path / "ensembl.gtf"
    pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
    return gtf_path


@pytest.fixture
def parquet_from_ensembl(tmp_path, ensembl_gtf_path):
    parquet_path = tmp_path / "ensembl.parquet"
    gtf_to_parquet(ensembl_gtf_path, parquet_path, preset=ENSEMBL_PRESET)
    return parquet_path


@pytest.fixture
def sample_df(parquet_from_ensembl):
    return query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)


class TestWriteGtf:
    def test_write_to_file(self, sample_df, tmp_path):
        out = tmp_path / "out.gtf"
        write_gtf(sample_df, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_written_gtf_is_readable_by_pyranges(self, sample_df, tmp_path):
        out = tmp_path / "out.gtf"
        write_gtf(sample_df, out)
        gr = pr.read_gtf(str(out))
        assert len(gr) == len(sample_df)

    def test_write_to_stdout(self, sample_df, capsys):
        write_gtf(sample_df, None)
        captured = capsys.readouterr()
        assert len(captured.out) > 0

    def test_row_count_preserved(self, sample_df, tmp_path):
        out = tmp_path / "out.gtf"
        write_gtf(sample_df, out)
        gr = pr.read_gtf(str(out))
        assert len(gr) == len(sample_df)


class TestWriteGff3:
    def test_write_to_file(self, sample_df, tmp_path):
        out = tmp_path / "out.gff3"
        write_gff3(sample_df, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_gff3_header_present(self, sample_df, tmp_path):
        out = tmp_path / "out.gff3"
        write_gff3(sample_df, out)
        content = out.read_text()
        assert "##gff-version 3" in content

    def test_write_to_stdout(self, sample_df, capsys):
        write_gff3(sample_df, None)
        captured = capsys.readouterr()
        assert len(captured.out) > 0
        assert "##gff-version 3" in captured.out


class TestWriteParquet:
    def test_write_produces_valid_parquet(self, sample_df, tmp_path):
        out = tmp_path / "out.parquet"
        write_parquet(sample_df, out)
        assert out.exists()
        table = pq.read_table(str(out))
        assert len(table) == len(sample_df)

    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_compression_options(self, sample_df, tmp_path, compression):
        out = tmp_path / f"out_{compression}.parquet"
        write_parquet(sample_df, out, compression=compression)
        assert out.exists()
        table = pq.read_table(str(out))
        assert len(table) == len(sample_df)


class TestWriteOutput:
    def test_dispatch_gtf(self, sample_df, tmp_path):
        out = tmp_path / "out.gtf"
        write_output(sample_df, out, "gtf")
        assert out.exists()

    def test_dispatch_gff3(self, sample_df, tmp_path):
        out = tmp_path / "out.gff3"
        write_output(sample_df, out, "gff3")
        assert out.exists()
        assert "##gff-version 3" in out.read_text()

    def test_dispatch_parquet(self, sample_df, tmp_path):
        out = tmp_path / "out.parquet"
        write_output(sample_df, out, "parquet")
        assert out.exists()
        assert pq.read_table(str(out)).num_rows == len(sample_df)

    def test_parquet_without_output_raises(self, sample_df):
        with pytest.raises(ValueError, match="output path"):
            write_output(sample_df, None, "parquet")

    def test_unknown_format_raises(self, sample_df, tmp_path):
        with pytest.raises(ValueError, match="Unknown output format"):
            write_output(sample_df, tmp_path / "out.bed", "bed")
