"""Tests for the gff2parquet CLI."""

import pytest
import pyranges1 as pr

from gff2parquet.cli import main


@pytest.fixture
def ensembl_gtf_path(tmp_path):
    gtf_path = tmp_path / "ensembl.gtf"
    pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
    return gtf_path


@pytest.fixture
def ensembl_parquet_path(ensembl_gtf_path, tmp_path):
    parquet_path = tmp_path / "ensembl.parquet"
    assert main(["build", str(ensembl_gtf_path), str(parquet_path)]) == 0
    return parquet_path


class TestBuildSubcommand:
    def test_no_args_prints_help(self):
        assert main([]) == 0

    def test_build_gtf_succeeds(self, ensembl_gtf_path, tmp_path):
        out = tmp_path / "out.parquet"
        rc = main(["build", str(ensembl_gtf_path), str(out)])
        assert rc == 0
        assert out.exists()

    def test_build_auto_detects_gtf(self, ensembl_gtf_path, tmp_path):
        """build dispatches to gtf_to_parquet based on .gtf extension."""
        out = tmp_path / "out.parquet"
        rc = main(["build", str(ensembl_gtf_path), str(out)])
        assert rc == 0

    def test_build_missing_input_returns_error(self, tmp_path):
        rc = main(
            ["build", str(tmp_path / "nonexistent.gtf"), str(tmp_path / "out.parquet")]
        )
        assert rc == 1

    def test_build_unknown_extension_returns_error(self, tmp_path):
        bad_input = tmp_path / "annotations.xyz"
        bad_input.write_text("data")
        rc = main(["build", str(bad_input), str(tmp_path / "out.parquet")])
        assert rc == 1

    def test_build_with_preset(self, ensembl_gtf_path, tmp_path):
        out = tmp_path / "out.parquet"
        rc = main(["build", str(ensembl_gtf_path), str(out), "--preset", "ensembl"])
        assert rc == 0

    def test_build_with_compression(self, ensembl_gtf_path, tmp_path):
        out = tmp_path / "out.parquet"
        rc = main(["build", str(ensembl_gtf_path), str(out), "--compression", "snappy"])
        assert rc == 0


class TestQuerySubcommand:
    def test_query_no_filters_exits_zero(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.gtf"
        rc = main(["query", str(ensembl_parquet_path), "--output", str(out)])
        assert rc == 0
        assert out.exists()

    def test_query_with_region(self, ensembl_parquet_path, tmp_path):
        import pyarrow.parquet as pq

        df = pq.read_table(str(ensembl_parquet_path)).to_pandas()
        chrom = str(df["Chromosome"].iloc[0])

        out = tmp_path / "result.gtf"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--region",
                chrom,
                "--output",
                str(out),
            ]
        )
        assert rc == 0

    def test_query_with_feature_filter(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.gtf"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--filter",
                "Feature",
                "eq",
                "gene",
                "--output",
                str(out),
            ]
        )
        assert rc == 0

    def test_query_parquet_output(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.parquet"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--output",
                str(out),
                "--format",
                "parquet",
            ]
        )
        assert rc == 0
        assert out.exists()

    def test_query_parquet_format_without_output_returns_error(
        self, ensembl_parquet_path
    ):
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--format",
                "parquet",
            ]
        )
        assert rc == 1

    def test_query_missing_input_returns_error(self, tmp_path):
        rc = main(["query", str(tmp_path / "nonexistent.parquet")])
        assert rc == 1

    def test_query_invalid_filter_op_returns_error(self, ensembl_parquet_path):
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--filter",
                "Feature",
                "like",
                "gene%",
            ]
        )
        assert rc == 1

    def test_query_strand_filter(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.parquet"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--strand",
                "plus",
                "--output",
                str(out),
                "--format",
                "parquet",
            ]
        )
        assert rc == 0
        import pyarrow.parquet as pq

        df = pq.read_table(str(out)).to_pandas()
        assert (df["Strand"] == "+").all()

    def test_query_format_auto_detected_from_extension(
        self, ensembl_parquet_path, tmp_path
    ):
        out = tmp_path / "result.parquet"
        rc = main(["query", str(ensembl_parquet_path), "--output", str(out)])
        assert rc == 0
        import pyarrow.parquet as pq

        pq.read_table(str(out))  # should not raise
