#!/usr/bin/env python3
"""Benchmark GTF vs Parquet read performance and file sizes."""

import argparse
import gc
import sys
import tempfile
import time
import tracemalloc
from contextlib import contextmanager
from pathlib import Path
from typing import Generator

import pandas as pd
import pyranges1 as pr

from gxf2parquet import gtf_to_parquet, read_gxf_parquet, GENCODE_PRESET


@contextmanager
def track_memory() -> Generator[dict, None, None]:
    """Context manager to track peak memory usage during execution.

    Usage:
        with track_memory() as mem_stats:
            # ... do work ...
        print(f"Peak memory: {mem_stats['peak'] / 1024 / 1024:.2f} MiB")

    Yields:
        dict with keys 'current' and 'peak' (both in bytes)
    """
    tracemalloc.start()
    stats = {"peak": 0, "current": 0}
    try:
        yield stats
    finally:
        current, peak = tracemalloc.get_traced_memory()
        stats["current"] = current
        stats["peak"] = peak
        tracemalloc.stop()


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


def benchmark_read_gtf(
    gtf_path: Path, n_runs: int = 3
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark pyranges GTF read time.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start = time.perf_counter()
            gr = pr.read_gtf(str(gtf_path))
            df = pd.DataFrame(gr)
            elapsed = time.perf_counter() - start

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


def benchmark_read_parquet(
    parquet_path: Path, n_runs: int = 3
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark Parquet read time.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start = time.perf_counter()
            df = read_gxf_parquet(parquet_path)
            elapsed = time.perf_counter() - start

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


def benchmark_filtered_read(
    parquet_path: Path,
    chromosome: str,
    feature: str = "gene",
    n_runs: int = 3,
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark filtered Parquet read time.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start = time.perf_counter()
            df = read_gxf_parquet(
                parquet_path,
                columns=[
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "gene_id",
                    "gene_name",
                ],
                filters=[("Chromosome", "==", chromosome), ("Feature", "==", feature)],
            )
            elapsed = time.perf_counter() - start

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


def benchmark_naive_filtered_read(
    gtf_path: Path,
    chromosome: str,
    feature: str = "gene",
    n_runs: int = 3,
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark naive approach: read full GTF then filter with pandas.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start = time.perf_counter()
            gr = pr.read_gtf(str(gtf_path))
            full_df = pd.DataFrame(gr)
            df = full_df[
                (full_df["Chromosome"] == chromosome) & (full_df["Feature"] == feature)
            ][["Chromosome", "Start", "End", "Strand", "gene_id", "gene_name"]]
            elapsed = time.perf_counter() - start

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


def benchmark_region_query(
    parquet_path: Path,
    chromosome: str,
    start: int,
    end: int,
    strand: str | None = None,
    n_runs: int = 3,
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark region query using Parquet with Start/End range filters.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start_time = time.perf_counter()

            # Build filters for region query
            filters = [
                ("Chromosome", "==", chromosome),
                ("Start", "<=", end),
                ("End", ">=", start),
            ]

            if strand is not None:
                filters.append(("Strand", "==", strand))

            df = read_gxf_parquet(
                parquet_path,
                columns=[
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "Feature",
                    "gene_id",
                    "gene_name",
                ],
                filters=filters,
            )
            elapsed = time.perf_counter() - start_time

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


def benchmark_naive_region_query(
    gtf_path: Path,
    chromosome: str,
    start: int,
    end: int,
    strand: str | None = None,
    n_runs: int = 3,
) -> tuple[float, pd.DataFrame, int]:
    """Benchmark naive region query: read full GTF then filter with gr.loci.

    Returns:
        tuple of (avg_time_seconds, resulting_dataframe, avg_peak_memory_bytes)
    """
    times = []
    peak_memories = []
    df = None

    for _ in range(n_runs):
        gc.collect()

        with track_memory() as mem_stats:
            start_time = time.perf_counter()
            gr = pr.read_gtf(str(gtf_path))

            # Use gr.loci for region filtering
            if strand is not None:
                filtered_gr = gr.loci[chromosome, strand, start:end]
            else:
                # For unstranded queries, filter both strands
                filtered_gr = gr.loci[chromosome, start:end]

            df = pd.DataFrame(filtered_gr)[
                [
                    "Chromosome",
                    "Start",
                    "End",
                    "Strand",
                    "Feature",
                    "gene_id",
                    "gene_name",
                ]
            ]
            elapsed = time.perf_counter() - start_time

        times.append(elapsed)
        peak_memories.append(mem_stats["peak"])

    avg_time = sum(times) / len(times)
    avg_peak_memory = sum(peak_memories) / len(peak_memories)
    return avg_time, df, avg_peak_memory


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
    parser.add_argument(
        "--region-chrom",
        help="Chromosome for region query benchmark (e.g., chr2)",
    )
    parser.add_argument(
        "--region-start",
        type=int,
        help="Start position for region query benchmark",
    )
    parser.add_argument(
        "--region-end",
        type=int,
        help="End position for region query benchmark",
    )
    parser.add_argument(
        "--region-strand",
        choices=["+", "-"],
        help="Strand for region query benchmark (optional, unstranded if not provided)",
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

        gtf_read_time, gtf_df, gtf_peak_memory = benchmark_read_gtf(
            args.gtf_path, args.n_runs
        )
        parquet_read_time, parquet_df, parquet_peak_memory = benchmark_read_parquet(
            parquet_path, args.n_runs
        )

        print(f"GTF read time (pyranges):  {gtf_read_time:.3f}s")
        print(f"GTF peak memory:           {format_size(gtf_peak_memory)}")
        print(f"Parquet read time:         {parquet_read_time:.3f}s")
        print(f"Parquet peak memory:       {format_size(parquet_peak_memory)}")
        print(f"Speedup (time):            {gtf_read_time / parquet_read_time:.1f}x")
        print(
            f"Memory reduction (peak):   {(1 - parquet_peak_memory / gtf_peak_memory) * 100:.1f}%"
        )

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
            print(
                f"Note: Using {filter_chrom} (requested {args.filter_chrom} not found)"
            )

        # Benchmark naive approach: read full GTF and filter with pandas
        naive_filtered_time, naive_filtered_df, naive_filtered_peak_memory = (
            benchmark_naive_filtered_read(
                args.gtf_path,
                filter_chrom,
                n_runs=args.n_runs,
            )
        )

        # Benchmark parquet filtered read
        filtered_read_time, filtered_df, filtered_peak_memory = benchmark_filtered_read(
            parquet_partitioned_path,
            filter_chrom,
            n_runs=args.n_runs,
        )

        filtered_memory = get_memory_usage(filtered_df)

        print(f"Naive approach (GTF + pandas filter):  {naive_filtered_time:.3f}s")
        print(
            f"Naive peak memory:                      {format_size(naive_filtered_peak_memory)}"
        )
        print(f"Parquet filtered read:                  {filtered_read_time:.3f}s")
        print(
            f"Parquet peak memory:                    {format_size(filtered_peak_memory)}"
        )
        print(
            f"Speedup (time):                         {naive_filtered_time / filtered_read_time:.1f}x"
        )
        print(
            f"Memory reduction (peak):                {(1 - filtered_peak_memory / naive_filtered_peak_memory) * 100:.1f}%"
        )
        print(f"\nFiltered rows:         {len(filtered_df):,}")
        print(f"Filtered memory:       {format_size(filtered_memory)}")
        print(
            f"Memory reduction:      {(1 - filtered_memory / parquet_memory) * 100:.1f}%"
        )

        # Region query comparison (if region parameters provided)
        if (
            args.region_chrom
            and args.region_start is not None
            and args.region_end is not None
        ):
            print("\n" + "=" * 60)
            strand_str = (
                f", {args.region_strand}" if args.region_strand else " (unstranded)"
            )
            print(
                f"REGION QUERY COMPARISON ({args.region_chrom}:{args.region_start}-{args.region_end}{strand_str})"
            )
            print("=" * 60)

            # Benchmark naive approach: read full GTF and filter with gr.loci
            naive_region_time, naive_region_df, naive_region_peak_memory = (
                benchmark_naive_region_query(
                    args.gtf_path,
                    args.region_chrom,
                    args.region_start,
                    args.region_end,
                    args.region_strand,
                    n_runs=args.n_runs,
                )
            )

            # Benchmark parquet region query
            region_time, region_df, region_peak_memory = benchmark_region_query(
                parquet_partitioned_path,
                args.region_chrom,
                args.region_start,
                args.region_end,
                args.region_strand,
                n_runs=args.n_runs,
            )

            region_memory = get_memory_usage(region_df)

            print(f"Naive approach (GTF + gr.loci):  {naive_region_time:.3f}s")
            print(
                f"Naive peak memory:               {format_size(naive_region_peak_memory)}"
            )
            print(f"Parquet region query:            {region_time:.3f}s")
            print(f"Parquet peak memory:             {format_size(region_peak_memory)}")
            print(
                f"Speedup (time):                  {naive_region_time / region_time:.1f}x"
            )
            print(
                f"Memory reduction (peak):         {(1 - region_peak_memory / naive_region_peak_memory) * 100:.1f}%"
            )
            print(f"\nRegion query rows:     {len(region_df):,}")
            print(f"Region query memory:   {format_size(region_memory)}")
            print(
                f"Memory reduction:      {(1 - region_memory / parquet_memory) * 100:.1f}%"
            )

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
