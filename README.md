# gff2parquet

Parse and transform Gene Transfer Format (GTF) annotation files to Apache Parquet format for more efficient and powerful downstream analysis.

## Motivation

I heavily use the pyranges1 library during my day-to-day analysis when working with transcriptome annotations and intervals. Loading a full GTF file into memory is fairly time intensive (~ 1 min), mainly due to the requirement to perform complex parsing of the attribute field to extract key-value pairs, which is inconvenient in interactive/exploratory analysis. I also regularly find myself only needing a subset of the intervals and metadata for analysis (e.g. exon intervals, protein coding genes), which is currently only possible by first reading and parsing the complete GTF file into a pyranges object (pandas dataframe) into memory before subsetting.

## Features

- **One-time parsing**: Parse the key-value pairs from the GTF attribute field into individual columns once for a given reference file, speedily read Parquet many times
- **Load what you need**: Leverage predicate pushdown to pre-filter for intervals of interest (e.g. chromosome, strand, exon/CDS entries etc.) and load only the columns (e.g. attribute keys) you need for analysis
- **Optimized dtypes**: Efficiently encode datatypes, reducing in-memory object size and avoiding per-run inference
- **Compatible with pyranges1**: import directly as a pyranges1 object for downstream analysis, relying on the same core dependencies (pandas, pyarrow) (TODO!)
- **Reduced disk space usage with Parquet vs uncompressed/gzipped TSV**

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

The benchmark script compares conversion time, file size, read performance, and memory usage between GTF and Parquet formats.

### Example: GENCODE v40 chr2 & chr20-22 Subset

```bash
# Create a smaller subset of multiple chromosomes for testing
grep 'chr2' gencode.v40.annotation.sorted.gtf > gencode.v40.chr2s.annotation.sorted.gtf

# Run benchmark with filtered reads
python benchmarks/benchmark.py gencode.v40.chr2s.annotation.sorted.gtf --filter-chrom chr2
```

**Results** (435,497 rows, 25 columns):

```
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
GTF read time (pyranges):  6.913s
Parquet read time:         0.253s
Speedup:                   27.3x

============================================================
IN-MEMORY SIZE COMPARISON
============================================================
GTF DataFrame memory:     511.06 MB
Parquet DataFrame memory: 384.33 MB

============================================================
FILTERED READ COMPARISON (chr2, gene)
============================================================
Filtered read time:    0.009s
Filtered rows:         4,267
Filtered memory:       594.93 KB
Memory reduction:      99.8%

============================================================
SUMMARY
============================================================
File size reduction:   97.3%
Read speedup:          27.3x
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

The package stores GTF-native coordinates (1-based, closed intervals). When reading GTF files with pyranges (which uses 0-based, half-open intervals), coordinates are automatically converted.

## License

MIT
