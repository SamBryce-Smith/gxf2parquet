# gxf2parquet

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
uv pip install gxf2parquet
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
- Run ruff for linting and python code formatting

## Quick Start

### CLI Usage

```bash
# Build: convert GENCODE GTF to Parquet (default gencode preset)
gff2parquet build gencode.v47.annotation.gtf.gz gencode.parquet

# Build: partition by Chromosome and Feature for faster region/feature queries
gff2parquet build gencode.v47.annotation.gtf.gz gencode.parquet \
    --partition-cols Chromosome Feature

# Build: Ensembl GTF with matching preset
gff2parquet build Homo_sapiens.GRCh38.113.gtf.gz ensembl.parquet --preset ensembl

# Build: GFF3 input (auto-detected from extension), snappy compression
gff2parquet build Homo_sapiens.GRCh38.113.gff3.gz ensembl.parquet \
    --preset ensembl --compression snappy

# Query: all genes on chr1, written as GTF to stdout
gff2parquet query gencode.parquet --region chr1 --filter Feature eq gene

# Query: region with strand filter, select subset of attribute columns, write to GTF file
# (all 8 core GTF columns must be included; extra attribute columns are optional)
gff2parquet query gencode.parquet \
    --region chr1:11869-14409 --strand plus \
    --columns Chromosome Source Feature Start End Score Strand Frame gene_name transcript_id \
    --output chr1_region.gtf

# Query: region selecting non-standard columns — write to Parquet (not GTF/GFF3)
gff2parquet query gencode.parquet \
    --region chr1:11869-14409 --strand plus \
    --columns Chromosome Start End Strand Feature gene_name \
    --output chr1_region.parquet

# Query: filter by feature set and gene type, save as Parquet for downstream use
gff2parquet query gencode.parquet \
    --filter Feature isin exon,CDS \
    --filter gene_type eq protein_coding \
    --output coding_exons.parquet

# Query: whole chromosome, GFF3 output
gff2parquet query gencode.parquet --region chr2 --output chr2.gff3
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

# Run benchmark with region query (e.g., chr2:30000-500000 on + strand)
uv run benchmarks/benchmark.py gencode.v40.chr2s.annotation.sorted.gtf \
    --filter-chrom chr2 \
    --region-chrom chr2 \
    --region-start 30000 \
    --region-end 500000 \
    --region-strand +
```

**Results** (435,497 rows, 25 columns):

```
Benchmarking: gencode.v40.chr2s.annotation.sorted.gtf
Runs per benchmark: 3
============================================================

Converting GTF to Parquet...
Conversion time: 8.97s

Converting GTF to Parquet (partitioned by Chromosome, Feature)...
Conversion time (partitioned): 8.62s

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
GTF read time (pyranges):  31.324s
GTF peak memory:           489.63 MB
Parquet read time:         0.341s
Parquet peak memory:       60.95 MB
Speedup (time):            91.8x
Memory reduction (peak):   87.6%

============================================================
IN-MEMORY SIZE COMPARISON
============================================================
GTF DataFrame memory:     511.06 MB
Parquet DataFrame memory: 287.09 MB

============================================================
FILTERED READ COMPARISON (chr2, gene)
============================================================
Naive approach (GTF + pandas filter):  31.605s
Naive peak memory:                      489.62 MB
Parquet filtered read:                  0.022s
Parquet peak memory:                    609.87 KB
Speedup (time):                         1423.5x
Memory reduction (peak):                99.9%

Filtered rows:         4,267
Filtered memory:       594.93 KB
Memory reduction:      99.8%

============================================================
REGION QUERY COMPARISON (chr2:30000-500000, +)
============================================================
Naive approach (GTF + gr.loci):  30.804s
Naive peak memory:               489.62 MB
Parquet region query:            0.019s
Parquet peak memory:             36.17 KB
Speedup (time):                  1589.5x
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
Read speedup:          91.8x
```

## Schema Presets

Two presets are available for common GTF sources:

- **GENCODE_PRESET**: Categorical columns for `gene_type`, `transcript_type`
- **ENSEMBL_PRESET**: Categorical columns for `gene_biotype`, `transcript_biotype`

## Coordinate System

The package stores GTF-native coordinates (1-based, closed intervals) in the Parquet file.

When reading with `read_gxf_parquet()`:
- **Default behavior** (`as_pyranges=True`): Returns a PyRanges object with 0-based, half-open coordinates (Start coordinate is automatically converted by subtracting 1)
- **DataFrame mode** (`as_pyranges=False`): Returns a pandas DataFrame with 1-based coordinates as stored in the Parquet file

This ensures compatibility with both the GTF standard and PyRanges conventions.

## License

MIT
