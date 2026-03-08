"""Tests for query_gff_parquet()."""

import pandas as pd
import pytest
import pyranges1 as pr

from gff2parquet import gtf_to_parquet, query_gff_parquet
from gff2parquet.schema import ENSEMBL_PRESET


@pytest.fixture
def parquet_from_ensembl(tmp_path):
    """Build a Parquet file from the pyranges ensembl_gtf example data."""
    gtf_path = tmp_path / "ensembl.gtf"
    pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
    parquet_path = tmp_path / "ensembl.parquet"
    gtf_to_parquet(gtf_path, parquet_path, preset=ENSEMBL_PRESET)
    return parquet_path


class TestQueryGffParquet:
    def test_returns_dataframe(self, parquet_from_ensembl):
        df = query_gff_parquet(parquet_from_ensembl)
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

    def test_coordinates_are_one_based(self, parquet_from_ensembl):
        df = query_gff_parquet(parquet_from_ensembl)
        # Start must be >= 1 (1-based storage)
        assert (df["Start"] >= 1).all()

    def test_column_selection(self, parquet_from_ensembl):
        cols = ["Chromosome", "Start", "End", "Feature"]
        df = query_gff_parquet(parquet_from_ensembl, columns=cols)
        assert list(df.columns) == cols

    def test_filter_by_feature(self, parquet_from_ensembl):
        df = query_gff_parquet(
            parquet_from_ensembl,
            filters=[("Feature", "==", "gene")],
        )
        assert len(df) > 0
        assert (df["Feature"] == "gene").all()

    def test_region_chromosome_only(self, parquet_from_ensembl):
        all_df = query_gff_parquet(parquet_from_ensembl)
        chrom = all_df["Chromosome"].iloc[0]

        df = query_gff_parquet(parquet_from_ensembl, regions=[chrom])
        assert len(df) > 0
        assert (df["Chromosome"] == chrom).all()

    def test_region_with_interval(self, parquet_from_ensembl):
        all_df = query_gff_parquet(parquet_from_ensembl)
        chrom = all_df["Chromosome"].iloc[0]
        # Use a very wide interval that should match everything on this chrom
        start = int(all_df["Start"].min())
        end = int(all_df["End"].max())

        df = query_gff_parquet(
            parquet_from_ensembl,
            regions=[f"{chrom}:{start}-{end}"],
        )
        assert len(df) > 0
        assert (df["Chromosome"] == chrom).all()
        assert (df["Start"] >= start).all()
        assert (df["End"] <= end).all()

    def test_strand_filter(self, parquet_from_ensembl):
        df_plus = query_gff_parquet(parquet_from_ensembl, strand="plus")
        df_minus = query_gff_parquet(parquet_from_ensembl, strand="minus")
        assert (df_plus["Strand"] == "+").all()
        assert (df_minus["Strand"] == "-").all()

    def test_multiple_regions_or_combined(self, parquet_from_ensembl):
        all_df = query_gff_parquet(parquet_from_ensembl)
        chroms = all_df["Chromosome"].unique()
        if len(chroms) < 2:
            pytest.skip("Need at least 2 chromosomes for this test")

        c1, c2 = chroms[:2]
        df = query_gff_parquet(parquet_from_ensembl, regions=[c1, c2])
        assert set(df["Chromosome"].unique()) <= {c1, c2}

    def test_no_results_returns_empty_dataframe(self, parquet_from_ensembl):
        df = query_gff_parquet(
            parquet_from_ensembl,
            regions=["chrNONEXISTENT"],
        )
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0
