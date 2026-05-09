# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a GTF to Parquet conversion utility for genomics workflows. It converts GTF (Gene Transfer Format) annotation files to Apache Parquet format to enable:
- One-time parsing with persistent column structure
- Predicate pushdown for efficient filtered reads
- Columnar storage with compression
- Explicit schema with optimized dtypes (categoricals)

## Quick Reference

```
src/gxf2parquet/
├── __init__.py       # Public API exports
├── convert.py:13     # gtf_to_parquet() function
├── read.py:10        # read_gxf_parquet() function
├── schema.py:19-41   # GENCODE_PRESET, ENSEMBL_PRESET, get_preset()
└── cli.py            # CLI entry point: gxf2parquet

tests/test_convert.py # 6 test classes, multiple test methods
benchmarks/benchmark.py # Performance comparison script
```

## API Reference

### `gtf_to_parquet()` (src/gxf2parquet/convert.py:13)

```python
def gtf_to_parquet(
    gtf_path: str | Path,
    parquet_path: str | Path,
    *,
    preset: SchemaPreset | None = None,      # Defaults to GENCODE_PRESET
    partition_cols: list[str] | None = None,  # e.g., ["Chromosome", "Feature"]
    compression: str = "zstd",                # Options: "zstd", "snappy", "gzip", "none"
) -> None
```

Converts a GTF file to Parquet format with optional partitioning.

### `read_gxf_parquet()` (src/gxf2parquet/read.py:10)

```python
def read_gxf_parquet(
    parquet_path: str | Path,
    *,
    columns: list[str] | None = None,
    filters: list[tuple] | list[list[tuple]] | None = None,
    as_pyranges: bool = True,  # False returns DataFrame with 1-based coords
) -> pr.PyRanges | pd.DataFrame
```

Reads a GTF Parquet file as a PyRanges object or pandas DataFrame. When `as_pyranges=True` (default), coordinates are converted to 0-based. When `as_pyranges=False`, coordinates remain 1-based as stored.

## Build & Development Commands

```bash
# Install in development mode
uv sync

# Install without dev dependencies
uv sync --no-dev

# Run all tests
uv run pytest

# Run specific test
uv run pytest tests/test_convert.py::TestRoundtrip::test_roundtrip

# Run benchmarks
uv run benchmarks/benchmark.py /path/to/annotations.gtf
```

## Pre-commit Hooks

```bash
# Install pre-commit hooks (after uv sync)
uv run prek install

# Run hooks manually on all files
uv run prek run --all-files
```

## Schema Presets

### GENCODE_PRESET (src/gxf2parquet/schema.py:19)

**Categorical columns:** Chromosome, Source, Feature, Strand, gene_type, transcript_type
**List columns:** tag, ont
**Use for:** GENCODE annotation files

### ENSEMBL_PRESET (src/gxf2parquet/schema.py:31)

**Categorical columns:** Chromosome, Source, Feature, Strand, gene_biotype, transcript_biotype
**List columns:** tag
**Use for:** Ensembl annotation files

### Usage in code

```python
from gxf2parquet import GENCODE_PRESET, ENSEMBL_PRESET, get_preset

# Direct import
gtf_to_parquet("file.gtf", "file.parquet", preset=GENCODE_PRESET)

# Or by name
preset = get_preset("ensembl")
gtf_to_parquet("file.gtf", "file.parquet", preset=preset)
```

## Coordinate System

The package stores GTF-native coordinates (1-based, closed intervals) in the Parquet file.

```python
# Stored in Parquet: GTF-native 1-based coordinates
# Row in file: Start=1001, End=2000

# Reading as PyRanges (default, as_pyranges=True)
gr = read_gxf_parquet("file.parquet")
gr.Start  # 1000 (converted to 0-based)
gr.End    # 2000 (unchanged, half-open end)

# Reading as DataFrame (as_pyranges=False)
df = read_gxf_parquet("file.parquet", as_pyranges=False)
df["Start"]  # 1001 (1-based, as stored)
df["End"]    # 2000 (1-based, as stored)
```

**Key conversion (pyranges to GTF):**
```python
parquet_start = pyranges_start + 1  # pyranges is 0-based
parquet_end = pyranges_end          # already correct (half-open end equals closed end)
```

## Common Query Patterns

```python
from gxf2parquet import read_gxf_parquet

# Filter by chromosome and feature type (predicate pushdown)
gr = read_gxf_parquet(
    "annotations.parquet",
    filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")]
)

# Load only specific columns
gr = read_gxf_parquet(
    "annotations.parquet",
    columns=["Chromosome", "Start", "End", "gene_name", "gene_type"]
)

# Get raw DataFrame for non-interval analysis
df = read_gxf_parquet("annotations.parquet", as_pyranges=False)

# Combine filters and column selection
gr = read_gxf_parquet(
    "annotations.parquet",
    columns=["Chromosome", "Start", "End", "gene_name"],
    filters=[("Chromosome", "in", ["chr1", "chr2"]), ("Feature", "==", "exon")]
)
```

## Key Design Decisions

- **Partitioning**: Default partitioning by `Chromosome` and `Feature` for common query patterns
- **Coordinates**: Store GTF-native coordinates (1-based, closed intervals), not pyranges coordinates (0-based, half-open)
- **Categoricals**: Columns like Chromosome, Source, Feature, Strand, gene_type, transcript_type are stored as categoricals
- **Multi-value attributes**: Columns like `tag` and `ont` are stored as `list<string>`

## Dependencies

- `pyranges1>=1.2.0` - GTF parsing
- `pandas<3.0.0` - DataFrame handling
- `pyarrow>=14.0.0` - Parquet read/write

## pyranges1 vs pyranges & accessing pyranges1 documentation

**Important:** pyranges1 has significant structural differences from the original pyranges and has renamed many methods. Always consult the official documentation rather than relying on training data.

### Accessing Documentation

Get a primer prompt for pyranges1 tasks:
```python
import pyranges1 as pr
pr.assistant.prompt()
```

Export complete pyranges1 documentation as a reference file:
```python
pr.assistant.export_docs("pr_docs.txt")
```

## CLI Reference

See README.md for comprehensive CLI usage examples. The CLI command is `gxf2parquet`.
