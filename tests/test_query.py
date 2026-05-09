"""Tests for query_gxf_parquet()."""

import pandas as pd
import pytest
import pyranges1 as pr

from gff2parquet import gtf_to_parquet, query_gxf_parquet
from gff2parquet.schema import ENSEMBL_PRESET


@pytest.fixture
def parquet_from_ensembl(tmp_path):
    """Build a Parquet file from the pyranges ensembl_gtf example data."""
    gtf_path = tmp_path / "ensembl.gtf"
    pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
    parquet_path = tmp_path / "ensembl.parquet"
    gtf_to_parquet(gtf_path, parquet_path, preset=ENSEMBL_PRESET)
    return parquet_path


class TestQueryGxfParquetDataFrame:
    """Tests with as_pyranges=False (returns DataFrame with 1-based coords)."""

    def test_returns_dataframe(self, parquet_from_ensembl):
        df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_coordinates_are_one_based(self, parquet_from_ensembl):
        df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        assert (df["Start"] >= 1).all()

    def test_column_selection(self, parquet_from_ensembl):
        cols = ["Chromosome", "Start", "End", "Feature"]
        df = query_gxf_parquet(parquet_from_ensembl, columns=cols, as_pyranges=False)
        assert list(df.columns) == cols

    def test_filter_by_feature(self, parquet_from_ensembl):
        df = query_gxf_parquet(
            parquet_from_ensembl,
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )
        assert len(df) > 0
        assert (df["Feature"] == "gene").all()

    def test_region_chromosome_only(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chrom = all_df["Chromosome"].iloc[0]

        df = query_gxf_parquet(parquet_from_ensembl, regions=[chrom], as_pyranges=False)
        assert len(df) > 0
        assert (df["Chromosome"] == chrom).all()

    def test_region_with_interval(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chrom = all_df["Chromosome"].iloc[0]
        start = int(all_df["Start"].min())
        end = int(all_df["End"].max())

        df = query_gxf_parquet(
            parquet_from_ensembl,
            regions=[f"{chrom}:{start}-{end}"],
            as_pyranges=False,
        )
        assert len(df) > 0
        assert (df["Chromosome"] == chrom).all()
        assert (df["Start"] >= start).all()
        assert (df["End"] <= end).all()

    def test_strand_filter(self, parquet_from_ensembl):
        df_plus = query_gxf_parquet(parquet_from_ensembl, strand="plus", as_pyranges=False)
        df_minus = query_gxf_parquet(parquet_from_ensembl, strand="minus", as_pyranges=False)
        assert (df_plus["Strand"] == "+").all()
        assert (df_minus["Strand"] == "-").all()

    def test_multiple_regions_or_combined(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chroms = all_df["Chromosome"].unique()
        if len(chroms) < 2:
            pytest.skip("Need at least 2 chromosomes for this test")

        c1, c2 = chroms[:2]
        df = query_gxf_parquet(parquet_from_ensembl, regions=[c1, c2], as_pyranges=False)
        assert set(df["Chromosome"].unique()) <= {c1, c2}

    def test_no_results_returns_empty_dataframe(self, parquet_from_ensembl):
        df = query_gxf_parquet(
            parquet_from_ensembl,
            regions=["chrNONEXISTENT"],
            as_pyranges=False,
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


class TestQueryGxfParquetPyRanges:
    """Tests with as_pyranges=True (default, returns PyRanges with 0-based coords)."""

    def test_returns_pyranges_by_default(self, parquet_from_ensembl):
        result = query_gxf_parquet(parquet_from_ensembl)
        assert isinstance(result, pr.PyRanges)

    def test_coordinates_are_zero_based(self, parquet_from_ensembl):
        df_1based = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        gr = query_gxf_parquet(parquet_from_ensembl, as_pyranges=True)
        assert (gr.Start.values == (df_1based["Start"] - 1).values).all()
        assert (gr.End.values == df_1based["End"].values).all()

    def test_region_filter_returns_pyranges(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chrom = all_df["Chromosome"].iloc[0]
        gr = query_gxf_parquet(parquet_from_ensembl, regions=[chrom])
        assert isinstance(gr, pr.PyRanges)
        assert len(gr) > 0


class TestQueryGxfParquetOutput:
    """Tests for the output= parameter (write to file, return None)."""

    def test_output_to_gtf_returns_none(self, parquet_from_ensembl, tmp_path):
        out = tmp_path / "out.gtf"
        result = query_gxf_parquet(parquet_from_ensembl, output=out)
        assert result is None
        assert out.exists()
        assert out.stat().st_size > 0

    def test_output_to_gff3_returns_none(self, parquet_from_ensembl, tmp_path):
        out = tmp_path / "out.gff3"
        result = query_gxf_parquet(parquet_from_ensembl, output=out)
        assert result is None
        assert out.exists()
        content = out.read_text()
        assert "##gff-version 3" in content

    def test_output_format_override(self, parquet_from_ensembl, tmp_path):
        out = tmp_path / "out.gtf"
        query_gxf_parquet(parquet_from_ensembl, output=out, output_format="gtf")
        assert out.exists()


class TestQueryGxfParquetPerRegionStrand:
    """Tests for per-region strand pairing."""

    def test_single_strand_broadcasts(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chroms = all_df["Chromosome"].unique()
        if len(chroms) < 2:
            pytest.skip("Need at least 2 chromosomes")

        c1, c2 = chroms[:2]
        df = query_gxf_parquet(
            parquet_from_ensembl,
            regions=[c1, c2],
            strand="+",
            as_pyranges=False,
        )
        assert (df["Strand"] == "+").all()

    def test_per_region_strand_list(self, parquet_from_ensembl):
        all_df = query_gxf_parquet(parquet_from_ensembl, as_pyranges=False)
        chroms = all_df["Chromosome"].unique()
        if len(chroms) < 2:
            pytest.skip("Need at least 2 chromosomes")

        c1, c2 = chroms[:2]
        df = query_gxf_parquet(
            parquet_from_ensembl,
            regions=[c1, c2],
            strand=["+", "-"],
            as_pyranges=False,
        )
        c1_rows = df[df["Chromosome"] == c1]
        c2_rows = df[df["Chromosome"] == c2]
        if len(c1_rows):
            assert (c1_rows["Strand"] == "+").all()
        if len(c2_rows):
            assert (c2_rows["Strand"] == "-").all()

    def test_strand_region_count_mismatch_raises(self, parquet_from_ensembl):
        with pytest.raises(ValueError, match="must match"):
            query_gxf_parquet(
                parquet_from_ensembl,
                regions=["chr1"],
                strand=["+", "-"],
            )
