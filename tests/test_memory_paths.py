"""Equivalence tests for the memory-saving conversion and read paths.

The column-by-column Arrow conversion, the self-destructing ``to_pandas`` and
the streamed CLI writers must give exactly the same results as the plain
pandas/pyarrow calls they replace.
"""

from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyranges1 as pr
import pytest

from gxf2parquet import gtf_to_parquet, query_gxf_parquet, read_gxf_parquet
from gxf2parquet.cli import main
from gxf2parquet.convert import _frame_to_table, _prepare_frame
from gxf2parquet.schema import BASE_PRESET, ENSEMBL_PRESET, GENCODE_PRESET

TESTS_DIR = Path(__file__).parent


def _ensembl_frame(tmp_path: Path) -> pd.DataFrame:
    gtf_path = tmp_path / "ensembl.gtf"
    pr.example_data.ensembl_gtf.to_gtf(str(gtf_path))
    return pd.DataFrame(pr.read_gtf(str(gtf_path), duplicate_attr=True))


FRAME_CASES = {
    "gencode_gtf": (
        lambda tmp: pd.DataFrame(
            pr.read_gtf(
                str(TESTS_DIR / "pyranges_data.gencode.gtf.gz"), duplicate_attr=True
            )
        ),
        GENCODE_PRESET,
    ),
    "gencode_gtf_base": (
        lambda tmp: pd.DataFrame(
            pr.read_gtf(
                str(TESTS_DIR / "pyranges_data.gencode.gtf.gz"), duplicate_attr=True
            )
        ),
        BASE_PRESET,
    ),
    "gencode_gff3": (
        lambda tmp: pd.DataFrame(
            pr.read_gff3(str(TESTS_DIR / "pyranges_data.gencode.gff.gz"))
        ),
        GENCODE_PRESET,
    ),
    "duplicate_tags": (
        lambda tmp: pd.DataFrame(
            pr.read_gtf(
                str(TESTS_DIR / "gencode.exampleduplicatetags.gtf"),
                duplicate_attr=True,
            )
        ),
        GENCODE_PRESET,
    ),
    "ensembl_gtf": (_ensembl_frame, ENSEMBL_PRESET),
}


class TestFrameToTable:
    @pytest.mark.parametrize("case", FRAME_CASES)
    def test_matches_from_pandas(self, case, tmp_path):
        make_frame, preset = FRAME_CASES[case]
        df = make_frame(tmp_path)
        _prepare_frame(df, preset)
        expected = pa.Table.from_pandas(df.copy(), preserve_index=False)

        actual = _frame_to_table(df)

        assert actual.schema.equals(expected.schema, check_metadata=False)
        assert actual.schema.pandas_metadata == expected.schema.pandas_metadata
        assert actual.equals(expected)
        # Restores the same dtypes on read.
        pd.testing.assert_frame_equal(actual.to_pandas(), expected.to_pandas())

    def test_consumes_frame(self, tmp_path):
        make_frame, preset = FRAME_CASES["gencode_gtf"]
        df = make_frame(tmp_path)
        _prepare_frame(df, preset)
        n_cols = df.shape[1]

        table = _frame_to_table(df)

        assert table.num_columns == n_cols
        assert df.columns.empty


@pytest.fixture(params=[None, ["Chromosome", "Feature"]], ids=["single", "partitioned"])
def gencode_parquet(request, tmp_path):
    path = tmp_path / "gencode.parquet"
    gtf_to_parquet(
        TESTS_DIR / "pyranges_data.gencode.gtf.gz",
        path,
        preset=GENCODE_PRESET,
        partition_cols=request.param,
    )
    return path


READ_ARGS = [
    {},
    {"columns": ["Chromosome", "Start", "End", "Strand", "gene_name"]},
    {"filters": [("Feature", "==", "exon")]},
    {"filters": [("gene_name", "in", ["DDX11L1", "WASH7P"])]},
]


class TestReadEquivalence:
    @pytest.mark.parametrize("kwargs", READ_ARGS)
    def test_dataframe_matches_plain_to_pandas(self, gencode_parquet, kwargs):
        expected = pq.read_table(str(gencode_parquet), **kwargs).to_pandas()
        actual = read_gxf_parquet(gencode_parquet, as_pyranges=False, **kwargs)
        pd.testing.assert_frame_equal(actual, expected)

    @pytest.mark.parametrize("kwargs", READ_ARGS)
    def test_pyranges_matches_plain_to_pandas(self, gencode_parquet, kwargs):
        expected = pq.read_table(str(gencode_parquet), **kwargs).to_pandas()
        expected["Start"] = expected["Start"] - 1
        actual = read_gxf_parquet(gencode_parquet, **kwargs)
        assert isinstance(actual, pr.PyRanges)
        pd.testing.assert_frame_equal(
            pd.DataFrame(actual), pd.DataFrame(pr.PyRanges(expected))
        )

    def test_query_matches_plain_to_pandas(self, gencode_parquet):
        filters = [("Chromosome", "==", "chr1"), ("Feature", "==", "gene")]
        expected = pq.read_table(str(gencode_parquet), filters=filters).to_pandas()
        actual = query_gxf_parquet(
            gencode_parquet,
            regions=["chr1"],
            filters=[("Feature", "==", "gene")],
            as_pyranges=False,
        )
        pd.testing.assert_frame_equal(
            actual.reset_index(drop=True), expected.reset_index(drop=True)
        )

        expected["Start"] = expected["Start"] - 1
        actual_gr = query_gxf_parquet(
            gencode_parquet, regions=["chr1"], filters=[("Feature", "==", "gene")]
        )
        pd.testing.assert_frame_equal(
            pd.DataFrame(actual_gr).reset_index(drop=True),
            pd.DataFrame(pr.PyRanges(expected)).reset_index(drop=True),
        )


class TestCliStreamedOutput:
    """With --output, gtf/gff3/bed are streamed to the file; content must match stdout."""

    @pytest.mark.parametrize("fmt", ["gtf", "gff3", "bed"])
    def test_file_matches_stdout(self, gencode_parquet, tmp_path, capsys, fmt):
        args = ["query", str(gencode_parquet), "--region", "chr1", "-of", fmt]
        assert main(args) == 0
        stdout = capsys.readouterr().out

        out_file = tmp_path / f"out.{fmt}"
        assert main([*args, "--output", str(out_file)]) == 0

        assert out_file.read_text() == stdout
        if fmt == "gff3":
            assert stdout.startswith("##gff-version 3\n")

    def test_gz_named_output_stays_plain_text(self, gencode_parquet, tmp_path):
        # Unchanged behaviour: the format comes from the name, the file is not gzipped.
        out_file = tmp_path / "out.gtf.gz"
        assert main(["query", str(gencode_parquet), "--output", str(out_file)]) == 0
        assert not out_file.read_bytes().startswith(b"\x1f\x8b")

    def test_parquet_output_matches_from_pandas(self, gencode_parquet, tmp_path):
        out_file = tmp_path / "out.parquet"
        assert main(["query", str(gencode_parquet), "--output", str(out_file)]) == 0

        expected = query_gxf_parquet(gencode_parquet)
        expected_df = pd.DataFrame(expected).copy()
        expected_df["Start"] = expected_df["Start"] + 1
        expected_table = pa.Table.from_pandas(expected_df, preserve_index=False)

        actual = pq.read_table(str(out_file))
        assert actual.schema.equals(expected_table.schema, check_metadata=False)
        assert actual.equals(expected_table)
