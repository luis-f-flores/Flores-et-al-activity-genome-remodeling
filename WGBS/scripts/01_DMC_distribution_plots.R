#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_DMC_distribution_plots.R
# Package: WGBS / methylation analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for WGBS DMC/DMR/DMG analysis, annotation,
#   overlap analysis, methylation-expression integration, signal tracks, and
#   OLS modeling.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# 01_DMC_distribution_plots.R
#
# Plot DMC distributions (volcano plots, histograms, density plots).
#
# Usage:
#   Rscript 01_DMC_distribution_plots.R \
#     --input-dir DMC_outputs \
#     --pattern "*_analysis/all_DMCs_*.txt" \
#     --outdir DMC_plots \
#     --diff-threshold 5 \
#     --p-threshold 0.001

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

input_dir <- get_arg("--input-dir")
pattern <- get_arg("--pattern", "*_analysis/all_DMCs_*.txt")
outdir <- get_arg("--outdir", "DMC_plots")
diff_threshold <- as.numeric(get_arg("--diff-threshold", "5"))
p_threshold <- as.numeric(get_arg("--p-threshold", "0.001"))

if (is.null(input_dir)) {
  stop("Usage: Rscript 01_DMC_distribution_plots.R --input-dir DMC_outputs")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all DMC files
files <- Sys.glob(file.path(input_dir, pattern))

if (length(files) == 0) {
  stop("No files found matching pattern: ", file.path(input_dir, pattern))
}

for (f in files) {
  contrast <- basename(dirname(f))
  contrast <- sub("DMC_", "", contrast)
  contrast <- sub("_analysis", "", contrast)

  df <- read_tsv(f, show_col_types = FALSE)

  if (!"meth.diff" %in% names(df) || !"pvalue" %in% names(df)) {
    message("Skipping ", f, ": missing required columns")
    next
  }

  df$significant <- abs(df$meth.diff) >= diff_threshold & df$pvalue <= p_threshold
  df$direction <- ifelse(df$meth.diff >= diff_threshold & df$pvalue <= p_threshold, "Hyper",
                  ifelse(df$meth.diff <= -diff_threshold & df$pvalue <= p_threshold, "Hypo", "NS"))

  # Volcano plot
  p_volcano <- ggplot(df, aes(x = meth.diff, y = -log10(pvalue), color = direction)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Hyper" = "red", "Hypo" = "blue", "NS" = "gray")) +
    geom_vline(xintercept = c(-diff_threshold, diff_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    labs(
      title = paste0("DMC Volcano Plot: ", contrast),
      x = "Methylation Difference (%)",
      y = "-log10(p-value)"
    ) +
    theme_bw() +
    theme(legend.position = "top")

  ggsave(
    file.path(outdir, paste0("volcano_", contrast, ".pdf")),
    p_volcano,
    width = 8,
    height = 6
  )

  # Methylation difference histogram
  p_hist <- ggplot(df %>% filter(significant), aes(x = meth.diff, fill = direction)) +
    geom_histogram(bins = 50, color = "black", alpha = 0.7) +
    scale_fill_manual(values = c("Hyper" = "red", "Hypo" = "blue")) +
    labs(
      title = paste0("Significant DMC Distribution: ", contrast),
      x = "Methylation Difference (%)",
      y = "Count"
    ) +
    theme_bw()

  ggsave(
    file.path(outdir, paste0("histogram_", contrast, ".pdf")),
    p_hist,
    width = 8,
    height = 5
  )

  message("Plotted: ", contrast)
}

message("[done] All DMC distribution plots generated.")
