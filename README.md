# gff2parquet

Parse and transform Gene Transfer Format (GTF) annotation files to Apache Parquet format for more efficient and powerful downstream analysis.

## Motivation

I heavily use the pyranges1 library during my day-to-day analysis when working with transcriptome annotations and intervals. Loading a full GTF file into memory is fairly time intensive (~ 1 min), mainly due to the requirement to perform complex parsing of the attribute field to extract key-value pairs, which is inconvenient in interactive/exploratory analysis. I also regularly find myself only needing a subset of the intervals and metadata for analysis (e.g. exon intervals, protein coding genes), which is currently only possible by first reading and parsing the complete GTF file into a pyranges object (pandas dataframe) into memory before subsetting.

## Features

- **One-time parsing**: Parse the key-value pairs from the GTF attribute field into individual columns once for a given reference file, speedily read Parquet many times
- **Load what you need**: Leverage predicate pushdown to pre-filter for intervals of interest (e.g. chromosome, strand, exon/CDS entries etc.) and load only the columns (e.g. attribute keys) you need for analysis
- **Optimized dtypes**: Efficiently encode datatypes, reducing in-memory object size and avoiding per-run inference
- **Compatible with pyranges1**: Returns pyranges1 objects by default for downstream analysis, relying on the same core dependencies (pandas, pyarrow)
- **Reduced disk space usage with Parquet vs uncompressed/gzipped TSV**

## Installation

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

### General users

```bash
uv pip install gff2parquet
```

### Developers

From the repository root, sync the project environment (includes dev dependencies by default):

```bash
uv sync
```

To install without dev dependencies:

```bash
uv sync --no-dev
```

#### Setting up pre-commit hooks

The project uses [prek](https://prek.j178.dev/) for pre-commit hooks. To enable automatic checks on every commit:

```bash
# Install the git hooks
uv run prek install

# Run hooks manually on all files (optional)
uv run prek run --all-files
```

The configured hooks will automatically run before each commit to:
- Validate TOML file syntax
- Detect accidentally committed private keys
- Ensure executable scripts have proper shebangs

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

# Read as PyRanges object (default) with 0-based coordinates
gr = read_gtf_parquet(
    "annotations.parquet",
    columns=["Chromosome", "Start", "End", "gene_name"],
    filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")],
)

# Or read as pandas DataFrame with 1-based coordinates
df = read_gtf_parquet(
    "annotations.parquet",
    columns=["Chromosome", "Start", "End", "gene_name"],
    filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")],
    as_pyranges=False,
)
```

## Testing

```bash
# Run all tests
uv run pytest

# Run with verbose output
uv run pytest -v

# Run a specific test
uv run pytest tests/test_convert.py::TestRoundtrip::test_roundtrip
```

## Benchmarking

The benchmark script compares conversion time, file size, read performance, and memory usage between GTF and Parquet formats.

### Example: GENCODE v40 chr2 & chr20-22 Subset

```bash
# Create a smaller subset of multiple chromosomes for testing
grep 'chr2' gencode.v40.annotation.sorted.gtf > gencode.v40.chr2s.annotation.sorted.gtf

# Run benchmark with filtered reads
uv run benchmarks/benchmark.py gencode.v40.chr2s.annotation.sorted.gtf --filter-chrom chr2
```

**Results** (435,497 rows, 25 columns):

```
Benchmarking: gencode.v40.chr2s.annotation.sorted.gtf
Runs per benchmark: 3
============================================================

Converting GTF to Parquet...
Conversion time: 8.55s

Converting GTF to Parquet (partitioned by Chromosome, Feature)...
Conversion time (partitioned): 8.44s

============================================================
FILE SIZE COMPARISON
============================================================
GTF file:              189.99 MB
Parquet file:          5.22 MB
Parquet (partitioned): 8.48 MB
Compression ratio:     2.7%

============================================================
FULL READ TIME COMPARISON
============================================================
GTF read time (pyranges):  7.241s
Parquet read time:         0.258s
Speedup:                   28.0x

============================================================
IN-MEMORY SIZE COMPARISON
============================================================
GTF DataFrame memory:     511.06 MB
Parquet DataFrame memory: 384.33 MB

============================================================
FILTERED READ COMPARISON (chr2, gene)
============================================================
Naive approach (GTF + pandas filter):  7.411s
Parquet filtered read:                  0.012s
Speedup:                                635.0x

Filtered rows:         4,267
Filtered memory:       594.93 KB
Memory reduction:      99.8%

============================================================
SUMMARY
============================================================
Rows in dataset:       435,497
Columns:               25
File size reduction:   97.3%
Read speedup:          28.0x
```

**Key Takeaways:**
- **97.3%** smaller file size with Parquet
- **27.3x** faster full reads compared to parsing GTF
- **99.8%** memory reduction with filtered reads using predicate pushdown
- Conversion time: ~8.5s for 435K rows

## Schema Presets

Two presets are available for common GTF sources:

- **GENCODE_PRESET**: Categorical columns for `gene_type`, `transcript_type`
- **ENSEMBL_PRESET**: Categorical columns for `gene_biotype`, `transcript_biotype`

## Coordinate System

The package stores GTF-native coordinates (1-based, closed intervals) in the Parquet file.

When reading with `read_gtf_parquet()`:
- **Default behavior** (`as_pyranges=True`): Returns a PyRanges object with 0-based, half-open coordinates (Start coordinate is automatically converted by subtracting 1)
- **DataFrame mode** (`as_pyranges=False`): Returns a pandas DataFrame with 1-based coordinates as stored in the Parquet file

This ensures compatibility with both the GTF standard and PyRanges conventions.

## License

MIT
