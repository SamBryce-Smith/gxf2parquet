"""Tests for GTF/GFF to Parquet conversion."""

import tempfile
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyranges1 as pr
import pytest

from gff2parquet import (
    ENSEMBL_PRESET,
    GENCODE_PRESET,
    gff_to_parquet,
    gtf_to_parquet,
    read_gtf_parquet,
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
def duplicate_tags_gtf_path():
    """Path to GTF file with duplicate tag/ont attributes."""
    return Path(__file__).parent / "gencode.example-duplicate-tags.gtf"


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
        roundtrip_df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

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
        roundtrip_gr = read_gtf_parquet(temp_parquet_path, as_pyranges=True)

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
        filtered_df = read_gtf_parquet(
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
        filtered_gr = read_gtf_parquet(
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
        df = read_gtf_parquet(
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
        gr = read_gtf_parquet(temp_parquet_path, columns=selected_cols)

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

        df = read_gtf_parquet(
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
            df = read_gtf_parquet(parquet_dir, as_pyranges=False)
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
            gr = read_gtf_parquet(parquet_dir)
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

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # Categorical columns should be category dtype
        for col in ["Chromosome", "Source", "Feature", "Strand"]:
            if col in df.columns:
                assert df[col].dtype.name == "category", f"{col} should be categorical"

    @pytest.mark.parametrize("gtf_fixture,preset", GTF_FIXTURES)
    def test_preset_as_pyranges(self, gtf_fixture, preset, temp_parquet_path, request):
        """Test preset works with corresponding data when returning PyRanges."""
        gtf_path = request.getfixturevalue(gtf_fixture)

        gtf_to_parquet(
            gtf_path,
            temp_parquet_path,
            preset=preset,
        )

        gr = read_gtf_parquet(temp_parquet_path)

        assert isinstance(gr, pr.PyRanges)
        # Categorical columns should be category dtype
        for col in ["Chromosome", "Source", "Feature", "Strand"]:
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

            df = read_gtf_parquet(parquet_path, as_pyranges=False)
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

            gr = read_gtf_parquet(parquet_path)
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

    def test_gff_roundtrip_as_pyranges(self, gencode_gff_path, temp_parquet_path):
        """Test basic GFF conversion and reading as PyRanges (default)."""
        # Convert GFF to Parquet
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
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

    def test_gff_roundtrip_preserves_data(self, gencode_gff_path, temp_parquet_path):
        """Test that GFF data is preserved through roundtrip when reading as DataFrame."""
        # Read original GFF with pyranges (pyranges1 is a DataFrame subclass)
        original_gr = pr.read_gff3(str(gencode_gff_path))
        original_df = pd.DataFrame(original_gr)

        # Convert to Parquet and read back as DataFrame
        gff_to_parquet(gencode_gff_path, temp_parquet_path, preset=GENCODE_PRESET)
        roundtrip_df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

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
        roundtrip_gr = read_gtf_parquet(temp_parquet_path, as_pyranges=True)

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

    def test_gff_filtered_read_by_feature(self, gencode_gff_path, temp_parquet_path):
        """Test reading GFF with feature filter as PyRanges."""
        gff_to_parquet(
            gencode_gff_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        # Read only gene features
        filtered_gr = read_gtf_parquet(
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
        df = read_gtf_parquet(
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

        gr = read_gtf_parquet(
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
            gr = read_gtf_parquet(parquet_dir)
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

            gr = read_gtf_parquet(parquet_path)
            assert isinstance(gr, pr.PyRanges)
            assert len(gr) > 0


class TestDuplicateTags:
    """Test that duplicate GTF attributes are stored as list<string> in Parquet."""

    def test_list_columns_are_array_like(
        self, duplicate_tags_gtf_path, temp_parquet_path
    ):
        """Tag and ont columns should contain array-like sequences, not plain strings.

        PyArrow stores list<string> columns and pandas materialises them as
        numpy.ndarray objects (not plain Python lists), so we check for the
        sequence type rather than str.
        """
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        for col in ("tag", "ont"):
            if col in df.columns:
                non_null = df[col].dropna()
                for v in non_null:
                    assert not isinstance(v, str), (
                        f"'{col}' values should not be plain strings; got {v!r}"
                    )
                    assert hasattr(v, "__len__"), (
                        f"'{col}' values should be sequence-like; got {type(v)!r}"
                    )

    def test_single_tag_stored_as_single_element_list(
        self, duplicate_tags_gtf_path, temp_parquet_path
    ):
        """Rows with one tag value should be stored as a one-element sequence."""
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # The first transcript row (ENST00000456328.2) has only tag "basic"
        transcripts = df[df["Feature"] == "transcript"]
        first_transcript = transcripts.iloc[0]
        assert list(first_transcript["tag"]) == ["basic"]

    def test_multiple_tags_stored_as_multi_element_list(
        self, duplicate_tags_gtf_path, temp_parquet_path
    ):
        """Rows with multiple tag values should be stored as a multi-element sequence."""
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # The second transcript (ENST00000450305.2) has tag "basic" and "Ensembl_canonical"
        multi_tag_rows = df[df["tag"].apply(lambda x: x is not None and len(x) > 1)]
        assert len(multi_tag_rows) > 0, "Expected rows with multiple tags"
        assert any(
            list(v) == ["basic", "Ensembl_canonical"] for v in multi_tag_rows["tag"]
        )

    def test_missing_tag_is_null(self, duplicate_tags_gtf_path, temp_parquet_path):
        """Rows without a tag attribute should have null (None) in the tag column."""
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        # Gene rows in this file have no tag attribute
        gene_rows = df[df["Feature"] == "gene"]
        assert gene_rows["tag"].isna().all(), "Gene rows should have null tag"

    def test_ont_multiple_values(self, duplicate_tags_gtf_path, temp_parquet_path):
        """Rows with multiple ont values should be stored as a multi-element sequence."""
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        df = read_gtf_parquet(temp_parquet_path, as_pyranges=False)

        if "ont" in df.columns:
            multi_ont_rows = df[df["ont"].apply(lambda x: x is not None and len(x) > 1)]
            assert len(multi_ont_rows) > 0, "Expected rows with multiple ont values"
            assert any(
                list(v) == ["PGO:0000005", "PGO:0000019"] for v in multi_ont_rows["ont"]
            )

    def test_parquet_schema_uses_list_type(
        self, duplicate_tags_gtf_path, temp_parquet_path
    ):
        """The Parquet schema should encode tag and ont as list<string>."""
        gtf_to_parquet(
            duplicate_tags_gtf_path,
            temp_parquet_path,
            preset=GENCODE_PRESET,
        )

        schema = pq.read_schema(temp_parquet_path)
        for col in ("tag", "ont"):
            if col in schema.names:
                field = schema.field(col)
                assert pa.types.is_list(field.type), (
                    f"'{col}' should be list type in Parquet schema, got {field.type}"
                )
