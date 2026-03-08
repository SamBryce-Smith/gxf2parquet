"""Tests for region/strand/filter parsing utilities."""

import pytest

from gff2parquet.filters import (
    build_filters,
    parse_filter,
    parse_region,
    parse_strand,
)


class TestParseRegion:
    def test_chromosome_only(self):
        assert parse_region("chr1") == [("Chromosome", "==", "chr1")]

    def test_chromosome_with_interval(self):
        assert parse_region("chr1:1000-2000") == [
            ("Chromosome", "==", "chr1"),
            ("Start", ">=", 1000),
            ("End", "<=", 2000),
        ]

    def test_chromosome_with_large_coords(self):
        result = parse_region("chr1:1-248956422")
        assert result[1] == ("Start", ">=", 1)
        assert result[2] == ("End", "<=", 248956422)

    def test_non_standard_chrom_name(self):
        assert parse_region("scaffold_1234:100-200")[0] == (
            "Chromosome",
            "==",
            "scaffold_1234",
        )

    def test_missing_end_raises(self):
        with pytest.raises(ValueError, match="Expected 'CHROM' or 'CHROM:START-END'"):
            parse_region("chr1:1000")

    def test_non_integer_coords_raises(self):
        with pytest.raises(ValueError, match="must be integers"):
            parse_region("chr1:abc-200")

    def test_start_less_than_one_raises(self):
        with pytest.raises(ValueError, match="start must be >= 1"):
            parse_region("chr1:0-100")

    def test_end_less_than_start_raises(self):
        with pytest.raises(ValueError, match="end must be >= start"):
            parse_region("chr1:500-100")


class TestParseStrand:
    @pytest.mark.parametrize(
        "alias,expected",
        [
            ("plus", "+"),
            ("minus", "-"),
            ("+", "+"),
            ("-", "-"),
            ("PLUS", "+"),
            ("MINUS", "-"),
        ],
    )
    def test_valid_strands(self, alias, expected):
        assert parse_strand(alias) == expected

    def test_invalid_strand_raises(self):
        with pytest.raises(ValueError, match="Unknown strand"):
            parse_strand("forward")


class TestParseFilter:
    @pytest.mark.parametrize(
        "op,arrow_op",
        [
            ("eq", "=="),
            ("ne", "!="),
            ("gt", ">"),
            ("lt", "<"),
            ("ge", ">="),
            ("le", "<="),
        ],
    )
    def test_scalar_operators(self, op, arrow_op):
        col, result_op, val = parse_filter("gene_type", op, "protein_coding")
        assert col == "gene_type"
        assert result_op == arrow_op
        assert val == "protein_coding"

    def test_isin_splits_on_comma(self):
        col, op, val = parse_filter("Chromosome", "isin", "chr1,chr2,chr3")
        assert op == "in"
        assert val == ["chr1", "chr2", "chr3"]

    def test_notin_splits_on_comma(self):
        col, op, val = parse_filter("Feature", "notin", "CDS,UTR")
        assert op == "not in"
        assert val == ["CDS", "UTR"]

    def test_isin_strips_whitespace(self):
        _, _, val = parse_filter("Chromosome", "isin", "chr1, chr2 , chr3")
        assert val == ["chr1", "chr2", "chr3"]

    def test_unknown_op_raises(self):
        with pytest.raises(ValueError, match="Unknown operator"):
            parse_filter("gene_type", "like", "protein%")


class TestBuildFilters:
    def test_no_args_returns_none(self):
        assert build_filters() is None

    def test_single_region(self):
        result = build_filters(regions=["chr1"])
        assert result == [("Chromosome", "==", "chr1")]

    def test_multiple_regions_returns_list_of_lists(self):
        result = build_filters(regions=["chr1", "chr2"])
        assert isinstance(result, list)
        assert isinstance(result[0], list)
        assert len(result) == 2

    def test_strand_only(self):
        result = build_filters(strand="plus")
        assert result == [("Strand", "==", "+")]

    def test_region_and_strand_combined(self):
        result = build_filters(regions=["chr1"], strand="plus")
        assert ("Strand", "==", "+") in result
        assert ("Chromosome", "==", "chr1") in result

    def test_extra_filters(self):
        extra = [("Feature", "==", "gene")]
        result = build_filters(extra_filters=extra)
        assert result == [("Feature", "==", "gene")]

    def test_region_strand_and_extra_combined(self):
        result = build_filters(
            regions=["chr1:1000-2000"],
            strand="minus",
            extra_filters=[("Feature", "==", "exon")],
        )
        assert ("Chromosome", "==", "chr1") in result
        assert ("Strand", "==", "-") in result
        assert ("Feature", "==", "exon") in result

    def test_multi_region_with_strand(self):
        result = build_filters(regions=["chr1", "chr2"], strand="plus")
        # Each group should contain the strand filter
        assert all(("Strand", "==", "+") in group for group in result)
