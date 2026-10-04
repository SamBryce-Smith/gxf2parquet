#!/usr/bin/env Rscript
# Diagnostic plots for the gxf2parquet benchmark pipeline.
#
# Usage:
#   Rscript plot.R bench_query_summary.tsv bench_build_summary.tsv disk.tsv <outdir>
#
# Writes:
#   01_query_wall_time.pdf  wall time per engine, faceted by query (median, min-max bars)
#   02_query_peak_rss.pdf   peak RSS per engine, faceted by query
#   03_ratio_vs_naive.pdf   pyranges1 / engine ratio for wall time and peak RSS
#   04_disk_size.pdf        on-disk size of the GTF, its gzip and each build
#   05_build_cost.pdf       build wall time and peak RSS per build

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(forcats)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("usage: plot.R <query_summary.tsv> <build_summary.tsv> <disk.tsv> <outdir>")
}
query_path <- args[[1]]
build_path <- args[[2]]
disk_path <- args[[3]]
outdir <- args[[4]]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

naive_engine <- "pyranges1"

read_tsv_quiet <- function(path, ...) {
  if (file.size(path) == 0) {
    return(tibble())
  }
  read_tsv(path, show_col_types = FALSE, ...)
}

save_plot <- function(plot, name, width = 10, height = 7) {
  ggsave(file.path(outdir, name), plot, width = width, height = height)
}

placeholder <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = msg) +
    theme_void()
}

# param mixes values like "all", "chr1" and "10": keep it as text.
query <- read_tsv_quiet(query_path, col_types = cols(param = col_character(), .default = col_guess()))
build <- read_tsv_quiet(build_path)
disk <- read_tsv_quiet(disk_path)

# Naive engine first, then builds in config order of appearance.
if (nrow(query) > 0) {
  query <- query |>
    mutate(
      engine = fct_relevel(factor(engine, levels = unique(engine)), naive_engine),
      cell = paste(query, param, sep = ": ")
    )
}

# ---- 01 / 02: per-query wall time and peak RSS ----------------------------

metric_plot <- function(df, metric, ylab) {
  med <- paste0(metric, "_median")
  lo <- paste0(metric, "_min")
  hi <- paste0(metric, "_max")
  ggplot(df, aes(x = engine, y = .data[[med]], fill = engine)) +
    geom_col(width = 0.7) +
    geom_errorbar(aes(ymin = .data[[lo]], ymax = .data[[hi]]), width = 0.25) +
    facet_wrap(vars(annotation, cell), scales = "free_y") +
    labs(x = NULL, y = ylab, fill = "engine") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
}

if (nrow(query) > 0) {
  save_plot(
    metric_plot(query, "wall_s", "wall time (s, GNU time; median with min-max)") +
      ggtitle("Query wall time"),
    "01_query_wall_time.pdf"
  )
  save_plot(
    metric_plot(query, "max_rss_gb", "peak RSS (GB; median with min-max)") +
      ggtitle("Query peak memory"),
    "02_query_peak_rss.pdf"
  )
} else {
  walk(c("01_query_wall_time.pdf", "02_query_peak_rss.pdf"),
       \(f) save_plot(placeholder("no query results"), f))
}

# ---- 03: ratio vs the naive pyranges1 workflow -----------------------------

if (nrow(query) > 0 && naive_engine %in% query$engine) {
  ratios <- query |>
    select(annotation, cell, engine, wall_s_median, max_rss_gb_median) |>
    pivot_longer(c(wall_s_median, max_rss_gb_median), names_to = "metric", values_to = "value") |>
    group_by(annotation, cell, metric) |>
    mutate(ratio = value[engine == naive_engine][1] / value) |>
    ungroup() |>
    filter(engine != naive_engine) |>
    mutate(metric = recode(metric,
                           wall_s_median = "wall time",
                           max_rss_gb_median = "peak RSS"))

  p3 <- ggplot(ratios, aes(x = cell, y = ratio, fill = engine)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_y_log10() +
    facet_grid(rows = vars(metric), cols = vars(annotation), scales = "free_y") +
    labs(x = NULL, y = "pyranges1 / engine (log scale; > 1 = better than naive)",
         fill = "engine", title = "Improvement over the naive pyranges1 workflow") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
  save_plot(p3, "03_ratio_vs_naive.pdf")
} else {
  save_plot(placeholder("no naive (pyranges1) results to compare against"), "03_ratio_vs_naive.pdf")
}

# ---- 04: disk size -----------------------------------------------------------

if (nrow(disk) > 0) {
  p4 <- disk |>
    mutate(
      mb = bytes / 1024^2,
      artifact = fct_reorder(artifact, mb, .desc = TRUE),
      label = sprintf("%.1f MB", mb)
    ) |>
    ggplot(aes(x = artifact, y = mb)) +
    geom_col(fill = "grey40", width = 0.7) +
    geom_text(aes(label = label), vjust = -0.3, size = 3) +
    facet_wrap(vars(annotation), scales = "free_y") +
    labs(x = NULL, y = "size on disk (MB)", title = "Disk usage") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p4, "04_disk_size.pdf", height = 5)
} else {
  save_plot(placeholder("no disk usage results"), "04_disk_size.pdf")
}

# ---- 05: build cost ----------------------------------------------------------

if (nrow(build) > 0) {
  p5 <- build |>
    select(annotation, build, wall_s_median, wall_s_min, wall_s_max,
           max_rss_gb_median, max_rss_gb_min, max_rss_gb_max) |>
    pivot_longer(-c(annotation, build),
                 names_to = c("metric", ".value"),
                 names_pattern = "(wall_s|max_rss_gb)_(median|min|max)") |>
    mutate(metric = recode(metric, wall_s = "wall time (s)", max_rss_gb = "peak RSS (GB)")) |>
    ggplot(aes(x = build, y = median)) +
    geom_col(fill = "grey40", width = 0.7) +
    geom_errorbar(aes(ymin = min, ymax = max), width = 0.25) +
    facet_grid(rows = vars(metric), cols = vars(annotation), scales = "free_y") +
    labs(x = NULL, y = NULL, title = "Build cost (gxf2parquet build; median with min-max)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p5, "05_build_cost.pdf", height = 6)
} else {
  save_plot(placeholder("no build results"), "05_build_cost.pdf")
}
