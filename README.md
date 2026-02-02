# gff2parquet

Convert GTF annotation files to Apache Parquet format for efficient genomics workflows.

## Features

- **One-time parsing**: Parse GTF once, read Parquet many times
- **Predicate pushdown**: Efficient filtered reads by chromosome, feature type, etc.
- **Columnar storage**: Compression and selective column loading
- **Optimized dtypes**: Categorical columns for memory efficiency

## Installation

Requires Python 3.12+.

```bash
# Create virtual environment with uv
uv venv --python 3.12 .venv
source .venv/bin/activate

# Install package
uv pip install -e .

# Or with dev dependencies for testing
uv pip install -e ".[dev]"
```

## Quick Start

### CLI Usage

```bash
# Basic conversion
gtf-to-parquet annotations.gtf annotations.parquet

# With partitioning for faster filtered reads
gtf-to-parquet annotations.gtf annotations.parquet --partition-cols Chromosome Feature

# Specify schema preset (gencode or ensembl)
gtf-to-parquet annotations.gtf annotations.parquet --preset ensembl

# Choose compression (zstd, snappy, gzip, none)
gtf-to-parquet annotations.gtf annotations.parquet --compression zstd
```

### Python API

```python
from gff2parquet import gtf_to_parquet, read_gtf_parquet, GENCODE_PRESET

# Convert GTF to Parquet
gtf_to_parquet(
    "annotations.gtf",
    "annotations.parquet",
    preset=GENCODE_PRESET,
    partition_cols=["Chromosome", "Feature"],
)

# Read with filters (predicate pushdown)
df = read_gtf_parquet(
    "annotations.parquet",
    columns=["Chromosome", "Start", "End", "gene_name"],
    filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")],
)
```

## Testing

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run a specific test
pytest tests/test_convert.py::TestRoundtrip::test_roundtrip
```

## Benchmarking

```bash
# Run benchmarks with a GTF file
python benchmarks/benchmark.py /path/to/gencode.v44.annotation.gtf
```

## Schema Presets

Two presets are available for common GTF sources:

- **GENCODE_PRESET**: Categorical columns for `gene_type`, `transcript_type`
- **ENSEMBL_PRESET**: Categorical columns for `gene_biotype`, `transcript_biotype`

## Coordinate System

The package stores GTF-native coordinates (1-based, closed intervals). When reading GTF files with pyranges (which uses 0-based, half-open intervals), coordinates are automatically converted.

## License

MIT
