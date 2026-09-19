"""Tests for the gxf2parquet CLI."""

import pytest
import pyranges1 as pr

from gxf2parquet.cli import main


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
                "-of",
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
                "-of",
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
                "-of",
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

    def test_query_gtf_preserves_core_columns_with_column_subset(
        self, ensembl_parquet_path, tmp_path
    ):
        """A narrow --columns must not drop real core column values (issue #30)."""
        out = tmp_path / "result.gtf"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--filter",
                "Feature",
                "eq",
                "gene",
                "--columns",
                "gene_name",
                "-of",
                "gtf",
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        lines = [
            ln for ln in out.read_text().splitlines() if ln and not ln.startswith("#")
        ]
        assert lines, "expected at least one gene row"
        for ln in lines:
            fields = ln.split("\t")
            # Feature (col 3) and Source (col 2) carry real values, not "."
            assert fields[2] == "gene"
            assert fields[1] != "."
            assert 'gene_name "' in fields[8]

    def test_query_bed_output_keeps_extra_column(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.bed"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--filter",
                "Feature",
                "eq",
                "gene",
                "--columns",
                "gene_name",
                "-of",
                "bed",
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        lines = [ln for ln in out.read_text().splitlines() if ln]
        assert lines
        # Standard 6 BED fields plus one extra (gene_name) = 7 tab-separated fields
        assert len(lines[0].split("\t")) == 7

    def test_query_bed_format_auto_detected_from_extension(
        self, ensembl_parquet_path, tmp_path
    ):
        out = tmp_path / "result.bed"
        rc = main(["query", str(ensembl_parquet_path), "--output", str(out)])
        assert rc == 0
        assert out.exists()
        assert len(out.read_text().splitlines()[0].split("\t")) >= 6

    def test_query_tsv_output_has_header_and_1based_start(
        self, ensembl_parquet_path, tmp_path
    ):
        out = tmp_path / "result.tsv"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--columns",
                "Chromosome",
                "Start",
                "End",
                "gene_name",
                "-of",
                "tsv",
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        lines = out.read_text().splitlines()
        # Header present, tab-separated, no leading index column
        assert lines[0].split("\t") == ["Chromosome", "Start", "End", "gene_name"]

    def test_query_tsv_zero_based_shifts_start(self, ensembl_parquet_path, tmp_path):
        cols = ["Chromosome", "Start", "End", "gene_name"]
        one_based = tmp_path / "one.tsv"
        zero_based = tmp_path / "zero.tsv"
        base_args = [
            "query",
            str(ensembl_parquet_path),
            "--columns",
            *cols,
            "-of",
            "tsv",
        ]
        assert main([*base_args, "--output", str(one_based)]) == 0
        assert main([*base_args, "--xsv-zero-based", "--output", str(zero_based)]) == 0

        import pandas as pd

        df1 = pd.read_csv(one_based, sep="\t")
        df0 = pd.read_csv(zero_based, sep="\t")
        assert (df1["Start"] - df0["Start"] == 1).all()

    def test_query_csv_output(self, ensembl_parquet_path, tmp_path):
        out = tmp_path / "result.csv"
        rc = main(
            [
                "query",
                str(ensembl_parquet_path),
                "--columns",
                "Chromosome",
                "Start",
                "End",
                "--output-format",
                "csv",
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        assert out.read_text().splitlines()[0] == "Chromosome,Start,End"
