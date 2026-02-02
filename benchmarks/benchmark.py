#!/usr/bin/env python3
"""Benchmark GTF vs Parquet read performance and file sizes."""

import argparse
import gc
import os
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd
import pyranges1 as pr

from gff2parquet import gtf_to_parquet, read_gtf_parquet, GENCODE_PRESET


def get_file_size(path: Path) -> int:
    """Get file or directory size in bytes."""
    if path.is_file():
        return path.stat().st_size
    elif path.is_dir():
        return sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    return 0


def format_size(size_bytes: int) -> str:
    """Format bytes as human-readable string."""
    for unit in ["B", "KB", "MB", "GB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.2f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.2f} TB"


def get_memory_usage(df: pd.DataFrame) -> int:
    """Get DataFrame memory usage in bytes."""
    return df.memory_usage(deep=True).sum()


def benchmark_read_gtf(gtf_path: Path, n_runs: int = 3) -> tuple[float, pd.DataFrame]:
    """Benchmark pyranges GTF read time."""
    times = []
    df = None

    for _ in range(n_runs):
        gc.collect()
        start = time.perf_counter()
        gr = pr.read_gtf(str(gtf_path))
        df = pd.DataFrame(gr)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    avg_time = sum(times) / len(times)
    return avg_time, df


def benchmark_read_parquet(
    parquet_path: Path, n_runs: int = 3
) -> tuple[float, pd.DataFrame]:
    """Benchmark Parquet read time."""
    times = []
    df = None

    for _ in range(n_runs):
        gc.collect()
        start = time.perf_counter()
        df = read_gtf_parquet(parquet_path)
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    avg_time = sum(times) / len(times)
    return avg_time, df


def benchmark_filtered_read(
    parquet_path: Path,
    chromosome: str,
    feature: str = "gene",
    n_runs: int = 3,
) -> tuple[float, pd.DataFrame]:
    """Benchmark filtered Parquet read time."""
    times = []
    df = None

    for _ in range(n_runs):
        gc.collect()
        start = time.perf_counter()
        df = read_gtf_parquet(
            parquet_path,
            columns=["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"],
            filters=[("Chromosome", "==", chromosome), ("Feature", "==", feature)],
        )
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    avg_time = sum(times) / len(times)
    return avg_time, df


def main():
    parser = argparse.ArgumentParser(description="Benchmark GTF vs Parquet performance")
    parser.add_argument("gtf_path", type=Path, help="Path to GTF file")
    parser.add_argument(
        "--n-runs", type=int, default=3, help="Number of benchmark runs (default: 3)"
    )
    parser.add_argument(
        "--filter-chrom",
        default="chr1",
        help="Chromosome for filtered read benchmark (default: chr1)",
    )
    args = parser.parse_args()

    if not args.gtf_path.exists():
        print(f"Error: GTF file not found: {args.gtf_path}", file=sys.stderr)
        return 1

    print(f"Benchmarking: {args.gtf_path}")
    print(f"Runs per benchmark: {args.n_runs}")
    print("=" * 60)

    with tempfile.TemporaryDirectory() as tmpdir:
        parquet_path = Path(tmpdir) / "benchmark.parquet"
        parquet_partitioned_path = Path(tmpdir) / "benchmark_partitioned"

        # Convert to Parquet (non-partitioned)
        print("\nConverting GTF to Parquet...")
        start = time.perf_counter()
        gtf_to_parquet(args.gtf_path, parquet_path, preset=GENCODE_PRESET)
        convert_time = time.perf_counter() - start
        print(f"Conversion time: {convert_time:.2f}s")

        # Convert to Parquet (partitioned)
        print("\nConverting GTF to Parquet (partitioned by Chromosome, Feature)...")
        start = time.perf_counter()
        gtf_to_parquet(
            args.gtf_path,
            parquet_partitioned_path,
            preset=GENCODE_PRESET,
            partition_cols=["Chromosome", "Feature"],
        )
        convert_partitioned_time = time.perf_counter() - start
        print(f"Conversion time (partitioned): {convert_partitioned_time:.2f}s")

        # File size comparison
        print("\n" + "=" * 60)
        print("FILE SIZE COMPARISON")
        print("=" * 60)

        gtf_size = get_file_size(args.gtf_path)
        parquet_size = get_file_size(parquet_path)
        parquet_partitioned_size = get_file_size(parquet_partitioned_path)

        print(f"GTF file:              {format_size(gtf_size)}")
        print(f"Parquet file:          {format_size(parquet_size)}")
        print(f"Parquet (partitioned): {format_size(parquet_partitioned_size)}")
        print(f"Compression ratio:     {parquet_size / gtf_size * 100:.1f}%")

        # Read time comparison
        print("\n" + "=" * 60)
        print("FULL READ TIME COMPARISON")
        print("=" * 60)

        gtf_read_time, gtf_df = benchmark_read_gtf(args.gtf_path, args.n_runs)
        parquet_read_time, parquet_df = benchmark_read_parquet(
            parquet_path, args.n_runs
        )

        print(f"GTF read time (pyranges):  {gtf_read_time:.3f}s")
        print(f"Parquet read time:         {parquet_read_time:.3f}s")
        print(f"Speedup:                   {gtf_read_time / parquet_read_time:.1f}x")

        # Memory usage comparison
        print("\n" + "=" * 60)
        print("IN-MEMORY SIZE COMPARISON")
        print("=" * 60)

        gtf_memory = get_memory_usage(gtf_df)
        parquet_memory = get_memory_usage(parquet_df)

        print(f"GTF DataFrame memory:     {format_size(gtf_memory)}")
        print(f"Parquet DataFrame memory: {format_size(parquet_memory)}")

        # Filtered read comparison
        print("\n" + "=" * 60)
        print(f"FILTERED READ COMPARISON ({args.filter_chrom}, gene)")
        print("=" * 60)

        # Check if the filter chromosome exists
        available_chroms = parquet_df["Chromosome"].unique()
        filter_chrom = args.filter_chrom
        if filter_chrom not in available_chroms:
            filter_chrom = available_chroms[0]
            print(f"Note: Using {filter_chrom} (requested {args.filter_chrom} not found)")

        filtered_read_time, filtered_df = benchmark_filtered_read(
            parquet_partitioned_path,
            filter_chrom,
            n_runs=args.n_runs,
        )

        filtered_memory = get_memory_usage(filtered_df)

        print(f"Filtered read time:    {filtered_read_time:.3f}s")
        print(f"Filtered rows:         {len(filtered_df):,}")
        print(f"Filtered memory:       {format_size(filtered_memory)}")
        print(f"Memory reduction:      {(1 - filtered_memory / parquet_memory) * 100:.1f}%")

        # Summary
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Rows in dataset:       {len(parquet_df):,}")
        print(f"Columns:               {len(parquet_df.columns)}")
        print(f"File size reduction:   {(1 - parquet_size / gtf_size) * 100:.1f}%")
        print(f"Read speedup:          {gtf_read_time / parquet_read_time:.1f}x")

    return 0


if __name__ == "__main__":
    sys.exit(main())
