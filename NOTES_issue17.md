# Issue #17 — CLI Refactor: `build` and `query` Subcommands

Tracking notes for implementing https://github.com/smsaladi/gff2parquet/issues/17
(or whichever the repo URL resolves to — see `gh issue view 17`).

## Summary

Refactor `gff2parquet` CLI into two subcommands (`build` and `query`) and add a
`query_gff_parquet()` library function. See the issue body for full design spec.

---

## Implementation Plan Progress

| Step | Status | Description |
|------|--------|-------------|
| **1** | ✅ Done (cc0abb5) | Refactor CLI to subcommand structure (`build` + placeholder `query`) |
| **2** | ✅ Done | Auto-detect input format in `build` (`.gtf` / `.gff3` / `.gz`) |
| **3** | ✅ Done | Store source format in Parquet metadata (`gff2parquet.source_format`) |
| **4** | ✅ Done | Region and filter parsing utilities (`filters.py`) |
| **5** | ✅ Done | Output writers (`write.py`: `write_gtf`, `write_gff3`, `write_parquet`) |
| **6** | ✅ Done | `query_gff_parquet()` library function (`query.py`) |
| **7** | ✅ Done | Wire up full `query` CLI subcommand |
| **8** | ✅ Done | Update exports (`__init__.py`) |
| **9** | ✅ Done | Tests (`test_filters.py`, `test_query.py`, `test_cli.py`, `test_metadata.py`) |

---

## Step 1 — Done

**Files changed:** `src/gff2parquet/cli.py`

- `main()` now creates argparse subparsers with `build` and `query`
- Existing conversion logic moved into `_cmd_build(args)`
- `_cmd_query(args)` is a placeholder: prints "not yet implemented", exits 1
- `gff2parquet` with no subcommand prints help and exits 0

Acceptance criteria confirmed:
- `gff2parquet build annotations.gtf out.parquet` works identically to old flat CLI ✓
- `gff2parquet query` prints placeholder message ✓
- `gff2parquet` (no args) prints help with both subcommands ✓
- All 61 existing tests still pass ✓

---

## Step 2 — Next

**Files to change:** `src/gff2parquet/cli.py`, `src/gff2parquet/convert.py`

Key tasks:
- Add `detect_format(path: Path) -> str` helper (strips `.gz`, matches `.gtf` → `"gtf"`,
  `.gff`/`.gff3` → `"gff3"`, else raises `ValueError`)
- Update `_cmd_build` to dispatch to `gtf_to_parquet()` or `gff_to_parquet()` based on result
- Remove hardcoded `gtf_to_parquet` import in favour of conditional dispatch

---

## Key Design Notes (from issue)

### Metadata key
`gff2parquet.source_format` → `"gtf"` or `"gff3"` stored in Parquet file-level metadata.

### Region format (step 4)
- `chr1` → `[("Chromosome", "==", "chr1")]`
- `chr1:1000-2000` → `[("Chromosome", "==", "chr1"), ("Start", ">=", 1000), ("End", "<=", 2000)]`
- Coordinates are **1-based closed** (GTF convention), same as what's stored in Parquet.
- Multiple regions → OR-groups in pyarrow filter format.

### Filter operations (step 4)
| CLI op | pyarrow op |
|--------|-----------|
| `eq`   | `==`      |
| `ne`   | `!=`      |
| `gt`   | `>`       |
| `lt`   | `<`       |
| `ge`   | `>=`      |
| `le`   | `<=`      |
| `isin` | `in`      |
| `notin`| `not in`  |

### Strand mapping (step 4/7)
- `plus` → `"+"`
- `minus` → `"-"`
- Single strand value → broadcast to all regions
- List of strands → must match region count exactly

### Output writers (step 5)
- Convert 1-based DataFrame back to PyRanges (subtract 1 from Start) before calling
  `gr.to_gtf()` / `gr.to_gff3()`
- For STDOUT: write to temp file then stream, or pass file-like if pyranges1 supports it

### New files to create
- `src/gff2parquet/filters.py` — region/filter parsing
- `src/gff2parquet/write.py` — output writers
- `src/gff2parquet/query.py` — `query_gff_parquet()`
- `tests/test_filters.py`, `tests/test_write.py`, `tests/test_query.py`,
  `tests/test_cli.py`, `tests/test_metadata.py`
