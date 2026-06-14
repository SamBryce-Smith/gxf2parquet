"""Tests for GTF/GFF to Parquet conversion."""

import tempfile
from pathlib import Path

import pandas as pd
import pyranges1 as pr
import pytest

from gxf2parquet import (
    ENSEMBL_PRESET,
    GENCODE_PRESET,
    gff_to_parquet,
    gtf_to_parquet,
    read_gxf_parquet,
)


@pytest.fixture
def ensembl_gtf_path():
    """Create a temporary GTF file from pyranges example data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        gtf_path = Path(tmpdir) / "test.gtf"
        gr = pr.example_data.ensembl_gtf
        gr.to_gtf(str(gtf_path))
        yield gtf_path


@pytest.fixture
def gencode_gtf_path():
    """Path to GENCODE GTF test file (gzipped)."""
    return Path(__file__).parent / "pyranges_data.gencode.gtf.gz"


@pytest.fixture
def gencode_gff_path():
    """Path to GENCODE GFF test file (gzipped)."""
    return Path(__file__).parent / "pyranges_data.gencode.gff.gz"


@pytest.fixture
def temp_parquet_path():
    """Create a temporary path for Parquet output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir) / "output.parquet"


# Parametrize GTF test data sources
GTF_FIXTURES = [
    ("ensembl_gtf_path", ENSEMBL_PRESET),
    ("gencode_gtf_path", GENCODE_PRESET),
]


class TestRoundtrip:
    """Test GTF -> Parquet -> DataFrame roundtrip."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_roundtrip(self, gtf_fixture, preset, temp_parquet_path, request):
        """Test basic conversion and reading as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        # Convert GTF to Parquet
        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read back as DataFrame
        df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)

        # Verify basic structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in df.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 1-based)
        assert df["Start"].min() >= 1, "Start coordinates should be 1-based"

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_roundtrip_as_pyranges(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test basic conversion and reading as PyRanges (default)."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        # Convert GTF to Parquet
        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read back as PyRanges (default behavior)
        gr = read_gxf_parquet(temp_parquet_path)

        # Verify it's a PyRanges object
        assert isinstance(gr, pr.PyRanges)
        assert len(gr) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in gr.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 0-based for PyRanges)
        assert gr.Start.min() >= 0, "Start coordinates should be 0-based in PyRanges"

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_roundtrip_preserves_data(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test that data is preserved through roundtrip when reading as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        # Read original GTF with pyranges (pyranges1 is a DataFrame subclass)
        original_gr = pr.read_gtf(str(gtf_path))
        original_df = pd.DataFrame(original_gr)

        # Convert to Parquet and read back as DataFrame
        gtf_to_parquet(gtf_path, temp_parquet_path, preset=preset)
        roundtrip_df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)

        # Same number of rows
        assert len(roundtrip_df) == len(original_df)

        # Check coordinate conversion
        # Original pyranges Start is 0-based, Parquet should be 1-based
        assert (roundtrip_df["Start"] == original_df["Start"] + 1).all()
        # End should be the same (half-open end equals closed end)
        assert (roundtrip_df["End"] == original_df["End"]).all()

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_roundtrip_pyranges_equals_original(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test that PyRanges from parquet equals original pr.read_gtf()."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        # Read original GTF with pyranges
        original_gr = pr.read_gtf(str(gtf_path))

        # Convert to Parquet and read back as PyRanges
        gtf_to_parquet(gtf_path, temp_parquet_path, preset=preset)
        roundtrip_gr = read_gxf_parquet(temp_parquet_path, as_pyranges=True)

        # Same number of rows
        assert len(roundtrip_gr) == len(original_gr)

        # Convert both to DataFrames for comparison
        original_df = (
            pd.DataFrame(original_gr)
            .sort_values(by=["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )
        roundtrip_df = (
            pd.DataFrame(roundtrip_gr)
            .sort_values(by=["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )

        # Coordinates should match exactly (both 0-based)
        assert (roundtrip_df["Start"] == original_df["Start"]).all()
        assert (roundtrip_df["End"] == original_df["End"]).all()

        # Check other key columns
        for col in ["Chromosome", "Strand", "Feature"]:
            if col in original_df.columns:
                assert (roundtrip_df[col] == original_df[col]).all(), f"{col} mismatch"


class TestFilteredRead:
    """Test filtered reading with predicate pushdown."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_filtered_read_by_chromosome(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test reading with chromosome filter as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read full data to find available chromosomes
        full_df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)
        chromosomes = full_df["Chromosome"].unique()

        if len(chromosomes) > 0:
            target_chrom = chromosomes[0]
            filtered_df = read_gxf_parquet(
                temp_parquet_path,
                filters=[("Chromosome", "==", target_chrom)],
                as_pyranges=False,
            )

            assert len(filtered_df) > 0
            assert (filtered_df["Chromosome"] == target_chrom).all()
            assert len(filtered_df) <= len(full_df)

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_filtered_read_by_chromosome_as_pyranges(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test reading with chromosome filter as PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read full data to find available chromosomes
        full_gr = read_gxf_parquet(temp_parquet_path)
        chromosomes = full_gr.Chromosome.unique()

        if len(chromosomes) > 0:
            target_chrom = chromosomes[0]
            filtered_gr = read_gxf_parquet(
                temp_parquet_path,
                filters=[("Chromosome", "==", target_chrom)],
            )

            assert isinstance(filtered_gr, pr.PyRanges)
            assert len(filtered_gr) > 0
            assert (filtered_gr.Chromosome == target_chrom).all()
            assert len(filtered_gr) <= len(full_gr)

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_filtered_read_by_feature(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test reading with feature filter as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read only gene features
        filtered_df = read_gxf_parquet(
            temp_parquet_path,
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )

        if len(filtered_df) > 0:
            assert (filtered_df["Feature"] == "gene").all()

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_filtered_read_by_feature_as_pyranges(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test reading with feature filter as PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        # Read only gene features
        filtered_gr = read_gxf_parquet(
            temp_parquet_path,
            filters=[("Feature", "==", "gene")],
        )

        if len(filtered_gr) > 0:
            assert isinstance(filtered_gr, pr.PyRanges)
            assert (filtered_gr.Feature == "gene").all()


class TestColumnSelection:
    """Test selective column loading."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_column_selection(self, gtf_fixture, preset, temp_parquet_path, request):
        """Test reading specific columns as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        selected_cols = ["Chromosome", "Start", "End", "gene_id"]
        df = read_gxf_parquet(
            temp_parquet_path, columns=selected_cols, as_pyranges=False
        )

        assert list(df.columns) == selected_cols

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_column_selection_as_pyranges(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test reading specific columns as PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        selected_cols = ["Chromosome", "Start", "End", "gene_id"]
        gr = read_gxf_parquet(temp_parquet_path, columns=selected_cols)

        assert isinstance(gr, pr.PyRanges)
        assert set(gr.columns) == set(selected_cols)

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_column_selection_with_filter(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test combined column selection and filtering as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        df = read_gxf_parquet(
            temp_parquet_path,
            columns=["Chromosome", "Start", "End"],
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )

        assert "Chromosome" in df.columns
        assert "Feature" not in df.columns  # Not selected

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_column_selection_with_filter_as_pyranges(
        self, gtf_fixture, preset, temp_parquet_path, request
    ):
        """Test combined column selection and filtering as PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        gr = read_gxf_parquet(
            temp_parquet_path,
            columns=["Chromosome", "Start", "End"],
            filters=[("Feature", "==", "gene")],
        )

        assert isinstance(gr, pr.PyRanges)
        assert "Chromosome" in gr.columns
        assert "Feature" not in gr.columns  # Not selected


class TestPartitioning:
    """Test partitioned writes."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_partitioned_write(self, gtf_fixture, preset, request):
        """Test writing with partitioning and reading as DataFrame."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_dir = Path(tmpdir) / "partitioned"

            gtf_to_parquet(
                gtf_path,
                parquet_dir,
                preset=preset,
                partition_cols=["Chromosome"],
            )

            # Should create directory structure
            assert parquet_dir.is_dir()

            # Read back partitioned data
            df = read_gxf_parquet(parquet_dir, as_pyranges=False)
            assert len(df) > 0

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_partitioned_write_as_pyranges(self, gtf_fixture, preset, request):
        """Test writing with partitioning and reading as PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_dir = Path(tmpdir) / "partitioned"

            gtf_to_parquet(
                gtf_path,
                parquet_dir,
                preset=preset,
                partition_cols=["Chromosome"],
            )

            # Should create directory structure
            assert parquet_dir.is_dir()

            # Read back partitioned data as PyRanges
            gr = read_gxf_parquet(parquet_dir)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0


class TestPresets:
    """Test schema presets."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_preset(self, gtf_fixture, preset, temp_parquet_path, request):
        """Test preset works with corresponding data."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)

        # All categorical columns declared in the preset should round-trip as category
        for col in preset.categorical_columns:
            if col in df.columns:
                assert df[col].dtype.name == "category", f"{col} should be categorical"

        # All integer columns declared in the preset should round-trip as integer dtype
        for col in [
            *preset.int16_columns,
            *preset.int32_columns,
            *preset.int64_columns,
        ]:
            if col in df.columns:
                assert pd.api.types.is_integer_dtype(df[col]), (
                    f"{col} should be integer dtype"
                )

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_preset_as_pyranges(self, gtf_fixture, preset, temp_parquet_path, request):
        """Test preset works with corresponding data when returning PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        gr = read_gxf_parquet(temp_parquet_path)

        assert isinstance(gr, pr.PyRanges)
        # All categorical columns declared in the preset should round-trip as category
        for col in preset.categorical_columns:
            if col in gr.columns:
                assert gr[col].dtype.name == "category", f"{col} should be categorical"


class TestCompression:
    """Test compression options."""

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_compression_options(self, gtf_fixture, preset, compression, request):
        """Test different compression codecs with DataFrame output."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_path = Path(tmpdir) / f"test_{compression}.parquet"

            gtf_to_parquet(
                gtf_path,
                parquet_path,
                preset=preset,
                compression=compression,
            )

            df = read_gxf_parquet(parquet_path, as_pyranges=False)
            assert len(df) > 0

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_compression_options_as_pyranges(
        self, gtf_fixture, preset, compression, request
    ):
        """Test different compression codecs with PyRanges output."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_path = Path(tmpdir) / f"test_{compression}.parquet"

            gtf_to_parquet(
                gtf_path,
                parquet_path,
                preset=preset,
                compression=compression,
            )

            gr = read_gxf_parquet(parquet_path)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0


# GFF-specific tests
class TestGFFRoundtrip:
    """Test GFF -> Parquet -> DataFrame roundtrip."""

    def test_gff_roundtrip(self, gencode_gff_path, temp_parquet_path):
        """Test basic GFF conversion and reading as DataFrame."""
        # Convert GFF to Parquet
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        # Read back as DataFrame
        df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)

        # Verify basic structure
        assert isinstance(df, pd.DataFrame)
        assert len(df) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in df.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 1-based)
        assert df["Start"].min() >= 1, "Start coordinates should be 1-based"

    def test_gff_roundtrip_as_pyranges(self, gencode_gff_path, temp_parquet_path):
        """Test basic GFF conversion and reading as PyRanges (default)."""
        # Convert GFF to Parquet
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        # Read back as PyRanges (default behavior)
        gr = read_gxf_parquet(temp_parquet_path)

        # Verify it's a PyRanges object
        assert isinstance(gr, pr.PyRanges)
        assert len(gr) > 0

        # Verify required columns exist
        required_cols = ["Chromosome", "Start", "End", "Strand", "Feature"]
        for col in required_cols:
            assert col in gr.columns, f"Missing column: {col}"

        # Verify coordinate conversion (Start should be 0-based for PyRanges)
        assert gr.Start.min() >= 0, "Start coordinates should be 0-based in PyRanges"

    def test_gff_roundtrip_preserves_data(self, gencode_gff_path, temp_parquet_path):
        """Test that GFF data is preserved through roundtrip when reading as DataFrame."""
        # Read original GFF with pyranges (pyranges1 is a DataFrame subclass)
        original_gr = pr.read_gff3(str(gencode_gff_path))
        original_df = pd.DataFrame(original_gr)

        # Convert to Parquet and read back as DataFrame
        gff_to_parquet(gencode_gff_path, temp_parquet_path, preset=GENCODE_PRESET)
        roundtrip_df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)

        # Same number of rows
        assert len(roundtrip_df) == len(original_df)

        # Check coordinate conversion
        # Original pyranges Start is 0-based, Parquet should be 1-based
        assert (roundtrip_df["Start"] == original_df["Start"] + 1).all()
        # End should be the same (half-open end equals closed end)
        assert (roundtrip_df["End"] == original_df["End"]).all()

    def test_gff_roundtrip_pyranges_equals_original(
        self, gencode_gff_path, temp_parquet_path
    ):
        """Test that PyRanges from parquet equals original pr.read_gff3()."""
        # Read original GFF with pyranges
        original_gr = pr.read_gff3(str(gencode_gff_path))

        # Convert to Parquet and read back as PyRanges
        gff_to_parquet(gencode_gff_path, temp_parquet_path, preset=GENCODE_PRESET)
        roundtrip_gr = read_gxf_parquet(temp_parquet_path, as_pyranges=True)

        # Same number of rows
        assert len(roundtrip_gr) == len(original_gr)

        # Convert both to DataFrames for comparison
        original_df = (
            pd.DataFrame(original_gr)
            .sort_values(by=["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )
        roundtrip_df = (
            pd.DataFrame(roundtrip_gr)
            .sort_values(by=["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )

        # Coordinates should match exactly (both 0-based)
        assert (roundtrip_df["Start"] == original_df["Start"]).all()
        assert (roundtrip_df["End"] == original_df["End"]).all()

        # Check other key columns
        for col in ["Chromosome", "Strand", "Feature"]:
            if col in original_df.columns:
                assert (roundtrip_df[col] == original_df[col]).all(), f"{col} mismatch"


class TestGFFFilteredRead:
    """Test filtered reading of GFF files with predicate pushdown."""

    def test_gff_filtered_read_by_chromosome(self, gencode_gff_path, temp_parquet_path):
        """Test reading GFF with chromosome filter as DataFrame."""
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        # Read full data to find available chromosomes
        full_df = read_gxf_parquet(temp_parquet_path, as_pyranges=False)
        chromosomes = full_df["Chromosome"].unique()

        if len(chromosomes) > 0:
            target_chrom = chromosomes[0]
            filtered_df = read_gxf_parquet(
                temp_parquet_path,
                filters=[("Chromosome", "==", target_chrom)],
                as_pyranges=False,
            )

            assert len(filtered_df) > 0
            assert (filtered_df["Chromosome"] == target_chrom).all()
            assert len(filtered_df) <= len(full_df)

    def test_gff_filtered_read_by_feature(self, gencode_gff_path, temp_parquet_path):
        """Test reading GFF with feature filter as PyRanges."""
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        # Read only gene features
        filtered_gr = read_gxf_parquet(
            temp_parquet_path,
            filters=[("Feature", "==", "gene")],
        )

        if len(filtered_gr) > 0:
            assert isinstance(filtered_gr, pr.PyRanges)
            assert (filtered_gr.Feature == "gene").all()


class TestGFFColumnSelection:
    """Test selective column loading for GFF files."""

    def test_gff_column_selection(self, gencode_gff_path, temp_parquet_path):
        """Test reading specific columns from GFF as DataFrame."""
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        selected_cols = ["Chromosome", "Start", "End", "gene_id"]
        df = read_gxf_parquet(
            temp_parquet_path, columns=selected_cols, as_pyranges=False
        )

        assert list(df.columns) == selected_cols

    def test_gff_column_selection_with_filter(
        self, gencode_gff_path, temp_parquet_path
    ):
        """Test combined column selection and filtering for GFF as PyRanges."""
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        gr = read_gxf_parquet(
            temp_parquet_path,
            columns=["Chromosome", "Start", "End"],
            filters=[("Feature", "==", "gene")],
        )

        assert isinstance(gr, pr.PyRanges)
        assert "Chromosome" in gr.columns
        assert "Feature" not in gr.columns  # Not selected


class TestGFFPartitioning:
    """Test partitioned writes for GFF files."""

    def test_gff_partitioned_write(self, gencode_gff_path):
        """Test writing GFF with partitioning and reading as PyRanges."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_dir = Path(tmpdir) / "partitioned"

            gff_to_parquet(
                gencode_gff_path,
                parquet_dir,
                preset=GENCODE_PRESET,
                partition_cols=["Chromosome"],
            )

            # Should create directory structure
            assert parquet_dir.is_dir()

            # Read back partitioned data as PyRanges
            gr = read_gxf_parquet(parquet_dir)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0


class TestGFFCompression:
    """Test compression options for GFF files."""

    @pytest.mark.parametrize("compression", ["zstd", "snappy", "gzip", "none"])
    def test_gff_compression_options(self, gencode_gff_path, compression):
        """Test different compression codecs with GFF files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            parquet_path = Path(tmpdir) / f"test_{compression}.parquet"

            gff_to_parquet(
                gencode_gff_path,
                parquet_path,
                preset=GENCODE_PRESET,
                compression=compression,
            )

            gr = read_gxf_parquet(parquet_path)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0
