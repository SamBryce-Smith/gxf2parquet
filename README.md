# gxf2parquet

Parse and transform Gene Transfer Format (GTF) annotation files to Apache Parquet format for more efficient and powerful downstream analysis.

## Contents

- [Motivation](#motivation)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
  - [CLI Usage](#cli-usage)
  - [Python API](#python-api)
- [Schema Presets](#schema-presets)
- [Coordinate System](#coordinate-system)
- [Benchmarks](#benchmarks)
- [Developer Guide](#developer-guide)
  - [Setting Up](#setting-up)
  - [Pre-commit Hooks](#pre-commit-hooks)
  - [Testing](#testing)
  - [Benchmarking](#benchmarking)
- [License](#license)

---

## Motivation

I heavily use the pyranges1 library during my day-to-day analysis when working with transcriptome annotations and intervals. Loading a full GTF file into memory is fairly time intensive (~ 1 min), mainly due to the requirement to perform complex parsing of the attribute field to extract key-value pairs, which is inconvenient in interactive/exploratory analysis. I also regularly find myself only needing a subset of the intervals and metadata for analysis (e.g. exon intervals, protein coding genes), which is currently only possible by first reading and parsing the complete GTF file into a pyranges object (pandas dataframe) into memory before subsetting.

[↑ Back to contents](#contents)

---

## Features

- **One-time parsing**: Parse the key-value pairs from the GTF attribute field into individual columns once for a given reference file, speedily read Parquet many times
- **Load what you need**: Leverage predicate pushdown to pre-filter for intervals of interest (e.g. chromosome, strand, exon/CDS entries etc.) and load only the columns (e.g. attribute keys) you need for analysis
- **Optimized dtypes**: Efficiently encode datatypes, reducing in-memory object size and avoiding per-run inference
- **Compatible with pyranges1**: Returns pyranges1 objects by default for downstream analysis, relying on the same core dependencies (pandas, pyarrow)
- **Reduced disk space usage with Parquet vs uncompressed/gzipped TSV**

[↑ Back to contents](#contents)

---

## Installation

Requires Python 3.12+.

```bash
uv pip install gxf2parquet
```

[↑ Back to contents](#contents)

---

## Quick Start

### CLI Usage

```bash
# Build: convert GENCODE GTF to Parquet (default gencode preset)
gxf2parquet build gencode.v47.annotation.gtf.gz gencode.parquet

# Build: partition by Chromosome and Feature for faster region/feature queries
gxf2parquet build gencode.v47.annotation.gtf.gz gencode.parquet \
    --partition-cols Chromosome Feature

# Build: Ensembl GTF with matching preset
gxf2parquet build Homo_sapiens.GRCh38.113.gtf.gz ensembl.parquet --preset ensembl

# Build: GFF3 input (auto-detected from extension), snappy compression
gxf2parquet build Homo_sapiens.GRCh38.113.gff3.gz ensembl.parquet \
    --preset ensembl --compression snappy

# Query: all genes on chr1, written as GTF to stdout
gxf2parquet query gencode.parquet --region chr1 --filter Feature eq gene

# Query: region with strand filter, keep a couple of attribute columns, write to GTF file.
# For gtf/gff3/bed output the core columns are ALWAYS included; --columns filters
# attribute/optional columns only (core columns listed here are ignored).
# Only intervals fully contained within the region are returned.
# Coordinates are assumed to follow GFF/GTF convention
gxf2parquet query gencode.parquet \
    --region chr1:11869-14409 --strand plus \
    --columns gene_name transcript_id \
    --output chr1_region.gtf

# Query: BED output — standard 6 columns plus gene_name as an extra field
gxf2parquet query gencode.parquet \
    --region chr1:11869-14409 --strand plus \
    --columns gene_name --output chr1_region.bed

# Query: TSV output with a header row (e.g. a tx2gene table).
# Coordinates are 1-based as stored; add --xsv-zero-based for BED-like 0-based Start.
gxf2parquet query gencode.parquet \
    --filter Feature eq transcript \
    --columns transcript_id gene_id gene_name -of tsv \
    --output tx2gene.tsv

# Query: region selecting non-standard columns — write to Parquet (breaks the GTF layout)
gxf2parquet query gencode.parquet \
    --region chr1:11869-14409 --strand plus \
    --columns Chromosome Start End Strand Feature gene_name \
    --output chr1_region.parquet

# Query: filter by feature set and gene type, save as Parquet for downstream use
gxf2parquet query gencode.parquet \
    --filter Feature isin exon,CDS \
    --filter gene_type eq protein_coding \
    --output coding_exons.parquet

# Query: whole chromosome, GFF3 output
gxf2parquet query gencode.parquet --region chr2 --output chr2.gff3
```

### Python API

```python
from gxf2parquet import gtf_to_parquet, read_gxf_parquet, GENCODE_PRESET

# Convert GTF to Parquet (partitioned for efficient filtered reads)
gtf_to_parquet(
    "gencode.v47.annotation.gtf.gz",
    "gencode.parquet",
    preset=GENCODE_PRESET,
    partition_cols=["Chromosome", "Feature"],
)

# Read as PyRanges object (default, 0-based coordinates) with predicate pushdown
gr = read_gtf_parquet(
    "gencode.parquet",
    columns=["Chromosome", "Start", "End", "gene_name", "gene_type"],
    filters=[("Chromosome", "==", "chr1"), ("Feature", "==", "gene")],
)

# Region query with strand — combines coordinate range and column filters
gr = read_gtf_parquet(
    "gencode.parquet",
    columns=["Chromosome", "Start", "End", "Strand", "Feature", "gene_name"],
    filters=[
        ("Chromosome", "==", "chr1"),
        ("Start", ">=", 11869),
        ("End", "<=", 14409),
        ("Strand", "==", "+"),
    ],
)

# Read as pandas DataFrame with 1-based coordinates (as stored)
df = read_gtf_parquet("gencode.parquet", as_pyranges=False)
```

[↑ Back to contents](#contents)

---

## Schema Presets

Two presets are available for common GTF sources:

- **GENCODE_PRESET**: Categorical columns for `gene_type`, `transcript_type`
- **ENSEMBL_PRESET**: Categorical columns for `gene_biotype`, `transcript_biotype`

[↑ Back to contents](#contents)

---

## Coordinate System

The package stores GTF-native coordinates (1-based, closed intervals) in the Parquet file.

When reading with `read_gxf_parquet()`:
- **Default behavior** (`as_pyranges=True`): Returns a PyRanges object with 0-based, half-open coordinates (Start coordinate is automatically converted by subtracting 1)
- **DataFrame mode** (`as_pyranges=False`): Returns a pandas DataFrame with 1-based coordinates as stored in the Parquet file

This ensures compatibility with both the GTF standard and PyRanges conventions.

[↑ Back to contents](#contents)

---

## Benchmarks

The benchmark below compares read performance and memory usage between GTF and Parquet formats.

Note(SBS): I'm not fully satisfied with the benchmarking workflow as presented here, it is likely to change significantly. Take the figures with a pinch of salt. At the moment, comparisons are deliberately kept 'biased' by comparing to the naive pyranges1 workflow. More appropriate benchmarks would include comparing against other competitors with comparable functionality e.g. tabix, GFFx, gffutils, polars-bio, gff2parquet (UriNeli) etc.).

### Example: GENCODE v40 chr2 & chr20-22 Subset

**Results** (435,497 rows, 25 columns):

```
Benchmarking: gencode.v40.chr2s.annotation.sorted.gtf
Runs per benchmark: 3
============================================================

Converting GTF to Parquet...
Conversion time: 6.46s

Converting GTF to Parquet (partitioned by Chromosome, Feature)...
Conversion time (partitioned): 5.93s

============================================================
FILE SIZE COMPARISON
============================================================
GTF file:              189.99 MB
Parquet file:          5.22 MB
Parquet (partitioned): 8.49 MB
Compression ratio:     2.7%

============================================================
FULL READ TIME COMPARISON
============================================================
GTF read time (pyranges):  13.962s
GTF peak memory:           614.87 MB
Parquet read time:         0.262s
Parquet peak memory:       60.95 MB
Speedup (time):            53.3x
Memory reduction (peak):   90.1%

============================================================
IN-MEMORY SIZE COMPARISON
============================================================
GTF DataFrame memory:     458.95 MB
Parquet DataFrame memory: 287.09 MB

============================================================
FILTERED READ COMPARISON (chr2, gene)
============================================================
Naive approach (GTF + pandas filter):  14.415s
Naive peak memory:                      614.87 MB
Parquet filtered read:                  0.017s
Parquet peak memory:                    610.07 KB
Speedup (time):                         835.5x
Memory reduction (peak):                99.9%

Filtered rows:         4,267
Filtered memory:       594.93 KB
Memory reduction:      99.8%

============================================================
REGION QUERY COMPARISON (chr2:30000-500000, +)
============================================================
Naive approach (GTF + gr.loci):  14.676s
Naive peak memory:               614.87 MB
Parquet region query:            0.016s
Parquet peak memory:             36.33 KB
Speedup (time):                  940.0x
Memory reduction (peak):         100.0%

Region query rows:     212
Region query memory:   30.70 KB
Memory reduction:      100.0%

============================================================
SUMMARY
============================================================
Rows in dataset:       435,497
Columns:               25
File size reduction:   97.3%
Read speedup:          53.3x
```

[↑ Back to contents](#contents)

---

## Developer Guide

### Setting Up

Requires Python 3.12+ and [uv](https://docs.astral.sh/uv/).

From the repository root, sync the project environment (includes dev dependencies by default):

```bash
uv sync
```

To install without dev dependencies:

```bash
uv sync --no-dev
```

### Pre-commit Hooks

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
- Run ruff for linting and python code formatting

### Testing

```bash
# Run all tests
uv run pytest

# Run with verbose output
uv run pytest -v

# Run a specific test
uv run pytest tests/test_convert.py::TestRoundtrip::test_roundtrip
```

### Benchmarking

```bash
# Create a smaller subset of multiple chromosomes for testing
grep 'chr2' gencode.v40.annotation.sorted.gtf > gencode.v40.chr2s.annotation.sorted.gtf

# Run benchmark with filtered reads
uv run benchmarks/benchmark.py gencode.v40.chr2s.annotation.sorted.gtf --filter-chrom chr2

# Run benchmark with region query (e.g., chr2:30000-500000 on + strand)
uv run benchmarks/benchmark.py gencode.v40.chr2s.annotation.sorted.gtf \
    --filter-chrom chr2 \
    --region-chrom chr2 \
    --region-start 30000 \
    --region-end 500000 \
    --region-strand +
```

[↑ Back to contents](#contents)

---

## License

MIT
