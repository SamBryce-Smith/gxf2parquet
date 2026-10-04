# gxf2parquet benchmark pipeline

A Snakemake pipeline, run in its own [pixi](https://pixi.sh) environment, that measures
**runtime**, **peak memory** and **disk usage** of gxf2parquet against the naive
pyranges1 workflow (parse the whole GTF with `pr.read_gtf`, then filter in memory).

The structure is adapted from the
[riker benchmark pipeline](https://github.com/fulcrumgenomics/riker/tree/main/benchmark-pipeline)
(MIT License, Copyright (c) 2026 Fulcrum Genomics LLC). Every timed cell is one process
wrapped in GNU `time -v`, and a `bench=100` resource lock makes sure only one timed job
runs at a time.

This directory is not part of the gxf2parquet package or its distributions.

## Quick start

```bash
cd benchmarks
pixi install                       # conda-forge/bioconda env + editable gxf2parquet from ..

pixi run smoke                     # whole DAG on tests/pyranges_data.gencode.gtf.gz (~5k rows)

pixi run download-gencode-v50      # download + md5-validate GENCODE v50 into data/
pixi run bench                     # full benchmark (config/performance.config.yaml)
```

`run.sh` wraps `pixi run snakemake` and passes on everything after `--` unchanged
(targets, `--config` overrides, other Snakemake flags):

```bash
./run.sh config/performance.config.yaml --dry-run
./run.sh config/performance.config.yaml --cores 4 -- --config cache_mode=evict replicates=5
./run.sh config/performance.config.yaml -- build_all     # build subworkflow only
./run.sh config/performance.config.yaml -- query_all     # query subworkflow only
```

The download tasks need a pixi version that supports task arguments. `pixi run download
<name>` works for any entry in `config/annotations.yaml`.

## What is measured

### Build subworkflow (`build_all`)

Each named parameter set under `builds:` in the config is passed to `gxf2parquet build`.
Each one is timed per annotation and replicate. Add entries under `builds:` to benchmark other
partitioning or compression choices; no code changes are needed:

```yaml
builds:
  single_zstd:        {partition_cols: [],                    compression: zstd}
  part_chrom_feature: {partition_cols: [Chromosome, Feature], compression: zstd}
```

The on-disk size of the decompressed GTF, the `.gtf.gz` and every build is recorded in
`disk.tsv`.

### Query subworkflow (`query_all`)

Engines:

| engine | input | how |
|---|---|---|
| `pyranges1` | decompressed GTF (`naive_input: gtf`) | `pr.read_gtf(path, duplicate_attr=True)`, then filter |
| `gxf2parquet-<build>` | rep1 output of each build in `query_builds` | `read_gxf_parquet(path, columns=..., filters=...)` |

Queries (from `queries:` in the config; each value is one timed cell). Both engines end
with an in-memory `pyranges1.PyRanges`:

| query | parameter | pyranges1 (naive) | gxf2parquet |
|---|---|---|---|
| `full` | `all` | read everything | read everything |
| `chrom` | chromosome | `gr[gr.Chromosome == c]` | `filters=[("Chromosome", "==", c)]` |
| `columns` | name in `column_sets:` | keep core GTF columns + set | `columns=CORE + set` |
| `gene_set` | N | `gr[gr.gene_name.isin(names)]` | `filters=[("gene_name", "in", names)]` |

For `gene_set`, N gene names are drawn once per annotation with a fixed seed
(`gene_set_seed`), so every engine and replicate queries the same set.

`workflow/scripts/run_query.py` records two timings. The GNU `time` wall clock includes
interpreter start-up and imports. `result.json` adds `import_s`, and `op_wall_s` for the
read + filter alone. `n_rows` must agree across engines for each (query, parameter);
`bench_query_summary.tsv` has an `n_rows_consistent` column, and the summary step prints
a warning when engines disagree.

## Cache mode

Linux keeps recently read file data in RAM (the *page cache*). The first read of a file
comes from disk; later reads of the same file come from memory. Snakemake runs the timed
cells one after another, so without a fixed policy whichever engine or replicate runs first
pays for the disk read and the rest get it free. The results would then depend on run
order rather than on the tool. This matters most for gxf2parquet, whose advantage partly
comes from reading fewer bytes (column pruning, predicate pushdown).

riker drops the whole OS cache before every timed run (needs root). Here
`cache_mode` picks the policy, applied by `workflow/scripts/cache_prep.py`, untimed, just before
each timed run:

| mode | what happens | measures | needs |
|---|---|---|---|
| `warm` (default) | every input file is read once | parse/filter cost with inputs in RAM | nothing; Linux + macOS |
| `evict` | `posix_fadvise(POSIX_FADV_DONTNEED)` on every input file | cold read of the files under test | Linux; no root |
| `drop` | `sync; echo 3 > /proc/sys/vm/drop_caches` via `sudo -n` | whole-system cold cache (riker) | passwordless sudo |

`evict` is the closest match to riker's cold-cache runs on a local machine without root.
Only the benchmarked files are evicted; Python and shared libraries stay cached. For `drop`,
allow just that one command, e.g. with `sudo visudo -f /etc/sudoers.d/drop-caches`:

```
<user> ALL=(root) NOPASSWD: /usr/bin/tee /proc/sys/vm/drop_caches
```

`sudo -n` makes a missing rule fail straight away instead of hanging on a password prompt
mid-run. The mode used is recorded in the `cache_mode` column of every output table.

## Annotations and checksums

`config/annotations.yaml` lists each annotation as either a `url` or a local `path`
(relative to this directory). `workflow/scripts/download_annotation.py` is used both by
`pixi run download <name>` and by the pipeline's staging rule. It:

1. downloads with `curl` to `data/<name>/<file>.partial`;
2. checks the md5 (`md5sum`) against the pinned `md5`, or, if none is pinned, against the
   line for the file in `md5sums_url` (GENCODE publishes an `MD5SUMS` per release);
3. runs `gzip -t`, then moves the file into place and writes a `<file>.md5` sidecar.

A failed check deletes the partial download. A file already present whose md5 still matches
is not downloaded again. After the first download, pin the printed checksum as `md5:` in
`annotations.yaml` so later runs don't rely on `MD5SUMS`.

## Outputs

Under `results_dir` (default `results/`; the smoke config uses `results-smoke/`):

| file | contents |
|---|---|
| `bench_query.tsv` | one row per (annotation, engine, query, param, rep): `wall_s`, `user_s`, `sys_s`, `cpu_percent`, `max_rss_kb`/`max_rss_gb`, `exit_status` (GNU time); `import_s`, `op_wall_s`, `n_rows`, `n_cols`, `df_mem_bytes` (run_query.py); `input_artifact`, `input_bytes`; `cache_mode`; package versions; `host_*` |
| `bench_build.tsv` | one row per (annotation, build, rep): `partition_cols`, `compression`, GNU time columns, `output_bytes`, `output_n_files`, versions, `host_*` |
| `disk.tsv` | one row per (annotation, artifact): `artifact` is `gtf`, `gtf_gz` or a build name; `bytes`, `n_files` |
| `bench_{query,build}_summary.tsv` | replicates collapsed to `<metric>_{median,min,max}` plus `rep_count`, `any_failed` (and `n_rows_consistent` for queries) |
| `plots/*.pdf` | query wall time, query peak RSS, ratio vs pyranges1, disk size, build cost |
| `run_config.json`, `host.json` | the resolved config and host description |

Per-run files are in `results/build/<annotation>/<build>/rep<N>/` and
`results/query/<annotation>/<engine>/<query>/<param>/rep<N>/`. Each run directory has
`time.txt`, `cmdline.txt` and `tool.log`; query runs also have `result.json`. Read
`tool.log` when a run fails.

## Layout

```
config/           annotations.yaml, performance + smoke configs
workflow/rules/   stage.smk, build.smk (build subworkflow), query.smk (query subworkflow), aggregate.smk
workflow/scripts/ download_annotation.py, cache_prep.py, run_query.py, sample_genes.py,
                  disk_usage.py, parse_gnu_time.py, host_info.py, merge_results.py,
                  summarize.py, plot.R
```
