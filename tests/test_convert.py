"""Tests for GTF to Parquet conversion."""

import tempfile
from pathlib import Path

import pandas as pd
import pyranges1 as pr
import pytest

from gff2parquet import gtf_to_parquet, read_gtf_parquet, ENSEMBL_PRESET


@pytest.fixture
def ensembl_gtf_path():
    """Create a temporary GTF file from pyranges example data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = Path(tmpdir) / "test.gtf"
        gr = pr.example_data.ensembl_gtf
        gr.to_gtf(str(gtf_path))
        yield gtf_path


@pytest.fixture
def temp_parquet_path():
    """Create a temporary path for Parquet output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir) / "output.parquet"


class TestRoundtrip:
    """Test GTF -> Parquet -> DataFrame roundtrip."""

    def test_roundtrip(self, ensembl_gtf_path, temp_parquet_path):
        """Test basic conversion and reading as DataFrame."""
        # Convert GTF to Parquet
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read back as DataFrame
        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # Verify basic structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in df.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 1-based)
        assert df["Start"].min() >= 1, "Start coordinates should be 1-based"

    def test_roundtrip_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test basic conversion and reading as PyRanges (default)."""
        # Convert GTF to Parquet
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read back as PyRanges (default behavior)
        gr = read_gtf_parquet(temp_parquet_path)

        # Verify it's a PyRanges object
        assert isinstance(gr, pr.PyRanges)
        assert len(gr) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in gr.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 0-based for PyRanges)
        assert gr.Start.min() >= 0, "Start coordinates should be 0-based in PyRanges"

    def test_roundtrip_preserves_data(self, ensembl_gtf_path, temp_parquet_path):
        """Test that data is preserved through roundtrip when reading as DataFrame."""
        # Read original GTF with pyranges (pyranges1 is a DataFrame subclass)
        original_gr = pr.read_gtf(str(ensembl_gtf_path))
        original_df = pd.DataFrame(original_gr)

        # Convert to Parquet and read back as DataFrame
        gtf_to_parquet(ensembl_gtf_path, temp_parquet_path, preset=ENSEMBL_PRESET)
        roundtrip_df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # Same number of rows
        assert len(roundtrip_df) == len(original_df)

        # Check coordinate conversion
        # Original pyranges Start is 0-based, Parquet should be 1-based
        assert (roundtrip_df["Start"] == original_df["Start"] + 1).all()
        # End should be the same (half-open end equals closed end)
        assert (roundtrip_df["End"] == original_df["End"]).all()

    def test_roundtrip_pyranges_equals_original(self, ensembl_gtf_path, temp_parquet_path):
        """Test that PyRanges from parquet equals original pr.read_gtf()."""
        # Read original GTF with pyranges
        original_gr = pr.read_gtf(str(ensembl_gtf_path))

        # Convert to Parquet and read back as PyRanges
        gtf_to_parquet(ensembl_gtf_path, temp_parquet_path, preset=ENSEMBL_PRESET)
        roundtrip_gr = read_gtf_parquet(temp_parquet_path, as_pyranges=True)

        # Same number of rows
        assert len(roundtrip_gr) == len(original_gr)

        # Convert both to DataFrames for comparison
        original_df = pd.DataFrame(original_gr).sort_values(
            by=["Chromosome", "Start", "End"]
        ).reset_index(drop=True)
        roundtrip_df = pd.DataFrame(roundtrip_gr).sort_values(
            by=["Chromosome", "Start", "End"]
        ).reset_index(drop=True)

        # Coordinates should match exactly (both 0-based)
        assert (roundtrip_df["Start"] == original_df["Start"]).all()
        assert (roundtrip_df["End"] == original_df["End"]).all()

        # Check other key columns
        for col in ["Chromosome", "Strand", "Feature"]:
            if col in original_df.columns:
                assert (roundtrip_df[col] == original_df[col]).all(), f"{col} mismatch"


class TestFilteredRead:
    """Test filtered reading with predicate pushdown."""

    def test_filtered_read_by_chromosome(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading with chromosome filter as DataFrame."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read full data to find available chromosomes
        full_df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)
        chromosomes = full_df["Chromosome"].unique()

        if len(chromosomes) > 0:
            target_chrom = chromosomes[0]
            filtered_df = read_gtf_parquet(
                temp_parquet_path,
                filters=[("Chromosome", "==", target_chrom)],
                as_pyranges=False,
            )

            assert len(filtered_df) > 0
            assert (filtered_df["Chromosome"] == target_chrom).all()
            assert len(filtered_df) <= len(full_df)

    def test_filtered_read_by_chromosome_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading with chromosome filter as PyRanges."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read full data to find available chromosomes
        full_gr = read_gtf_parquet(temp_parquet_path)
        chromosomes = full_gr.Chromosome.unique()

        if len(chromosomes) > 0:
            target_chrom = chromosomes[0]
            filtered_gr = read_gtf_parquet(
                temp_parquet_path,
                filters=[("Chromosome", "==", target_chrom)],
            )

            assert isinstance(filtered_gr, pr.PyRanges)
            assert len(filtered_gr) > 0
            assert (filtered_gr.Chromosome == target_chrom).all()
            assert len(filtered_gr) <= len(full_gr)

    def test_filtered_read_by_feature(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading with feature filter as DataFrame."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read only gene features
        filtered_df = read_gtf_parquet(
            temp_parquet_path,
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )

        if len(filtered_df) > 0:
            assert (filtered_df["Feature"] == "gene").all()

    def test_filtered_read_by_feature_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading with feature filter as PyRanges."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        # Read only gene features
        filtered_gr = read_gtf_parquet(
            temp_parquet_path,
            filters=[("Feature", "==", "gene")],
        )

        if len(filtered_gr) > 0:
            assert isinstance(filtered_gr, pr.PyRanges)
            assert (filtered_gr.Feature == "gene").all()


class TestColumnSelection:
    """Test selective column loading."""

    def test_column_selection(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading specific columns as DataFrame."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        selected_cols = ["Chromosome", "Start", "End", "gene_id"]
        df = read_gtf_parquet(temp_parquet_path, columns=selected_cols, as_pyranges=False)

        assert list(df.columns) == selected_cols

    def test_column_selection_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test reading specific columns as PyRanges."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        selected_cols = ["Chromosome", "Start", "End", "gene_id"]
        gr = read_gtf_parquet(temp_parquet_path, columns=selected_cols)

        assert isinstance(gr, pr.PyRanges)
        assert set(gr.columns) == set(selected_cols)

    def test_column_selection_with_filter(self, ensembl_gtf_path, temp_parquet_path):
        """Test combined column selection and filtering as DataFrame."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        df = read_gtf_parquet(
            temp_parquet_path,
            columns=["Chromosome", "Start", "End"],
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )

        assert "Chromosome" in df.columns
        assert "Feature" not in df.columns  # Not selected

    def test_column_selection_with_filter_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test combined column selection and filtering as PyRanges."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        gr = read_gtf_parquet(
            temp_parquet_path,
            columns=["Chromosome", "Start", "End"],
            filters=[("Feature", "==", "gene")],
        )

        assert isinstance(gr, pr.PyRanges)
        assert "Chromosome" in gr.columns
        assert "Feature" not in gr.columns  # Not selected


class TestPartitioning:
    """Test partitioned writes."""

    def test_partitioned_write(self, ensembl_gtf_path):
        """Test writing with partitioning and reading as DataFrame."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_dir = Path(tmpdir) / "partitioned"

            gtf_to_parquet(
                ensembl_gtf_path,
                parquet_dir,
                preset=ENSEMBL_PRESET,
                partition_cols=["Chromosome"],
            )

            # Should create directory structure
            assert parquet_dir.is_dir()

            # Read back partitioned data
            df = read_gtf_parquet(parquet_dir, as_pyranges=False)
            assert len(df) > 0

    def test_partitioned_write_as_pyranges(self, ensembl_gtf_path):
        """Test writing with partitioning and reading as PyRanges."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_dir = Path(tmpdir) / "partitioned"

            gtf_to_parquet(
                ensembl_gtf_path,
                parquet_dir,
                preset=ENSEMBL_PRESET,
                partition_cols=["Chromosome"],
            )

            # Should create directory structure
            assert parquet_dir.is_dir()

            # Read back partitioned data as PyRanges
            gr = read_gtf_parquet(parquet_dir)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0


class TestPresets:
    """Test schema presets."""

    def test_ensembl_preset(self, ensembl_gtf_path, temp_parquet_path):
        """Test ENSEMBL_PRESET works with example data."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # Categorical columns should be category dtype
        for col in ["Chromosome", "Source", "Feature", "Strand"]:
            if col in df.columns:
                assert df[col].dtype.name == "category", f"{col} should be categorical"

    def test_ensembl_preset_as_pyranges(self, ensembl_gtf_path, temp_parquet_path):
        """Test ENSEMBL_PRESET works with example data when returning PyRanges."""
        gtf_to_parquet(
            ensembl_gtf_path,
            temp_parquet_path,
            preset=ENSEMBL_PRESET,
        )

        gr = read_gtf_parquet(temp_parquet_path)

        assert isinstance(gr, pr.PyRanges)
        # Categorical columns should be category dtype
        for col in ["Chromosome", "Source", "Feature", "Strand"]:
            if col in gr.columns:
                assert gr[col].dtype.name == "category", f"{col} should be categorical"


class TestCompression:
    """Test compression options."""

    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_compression_options(self, ensembl_gtf_path, compression):
        """Test different compression codecs with DataFrame output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_path = Path(tmpdir) / f"test_{compression}.parquet"

            gtf_to_parquet(
                ensembl_gtf_path,
                parquet_path,
                preset=ENSEMBL_PRESET,
                compression=compression,
            )

            df = read_gtf_parquet(parquet_path, as_pyranges=False)
            assert len(df) > 0

    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_compression_options_as_pyranges(self, ensembl_gtf_path, compression):
        """Test different compression codecs with PyRanges output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_path = Path(tmpdir) / f"test_{compression}.parquet"

            gtf_to_parquet(
                ensembl_gtf_path,
                parquet_path,
                preset=ENSEMBL_PRESET,
                compression=compression,
            )

            gr = read_gtf_parquet(parquet_path)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0
