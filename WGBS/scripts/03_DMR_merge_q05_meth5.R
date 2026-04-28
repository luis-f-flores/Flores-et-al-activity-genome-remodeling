#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_DMR_merge_q05_meth5.R
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

# 03_DMR_merge_q05_meth5.R
#
# Merge neighboring DMRs and apply q-value and methylation difference thresholds.
#
# Usage:
#   Rscript 03_DMR_merge_q05_meth5.R \
#     --input-dir methylKit_DMR \
#     --pattern "DMRs_*_hyper_*.txt" \
#     --max-gap 500 \
#     --qvalue 0.05 \
#     --diff 5 \
#     --outdir DMR_merged

suppressPackageStartupMessages({
  library(GenomicRanges)
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
pattern <- get_arg("--pattern", "DMRs_*_q*.txt")
max_gap <- as.integer(get_arg("--max-gap", "500"))
q_cut <- as.numeric(get_arg("--qvalue", "0.05"))
diff_cut <- as.numeric(get_arg("--diff", "5"))
outdir <- get_arg("--outdir", "DMR_merged")

if (is.null(input_dir)) {
  stop("Usage: Rscript 03_DMR_merge_q05_meth5.R --input-dir methylKit_DMR")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all DMR files
files <- Sys.glob(file.path(input_dir, "**", pattern), dirmark = FALSE)

if (length(files) == 0) {
  stop("No files found matching pattern: ", pattern)
}

for (f in files) {
  contrast <- basename(dirname(f))
  file_base <- tools::file_path_sans_ext(basename(f))
  
  direction <- ifelse(grepl("hyper", file_base, ignore.case = TRUE), "hyper", "hypo")

  dmr_data <- read_tsv(f, show_col_types = FALSE)

  # Filter by thresholds
  dmr_filtered <- dmr_data %>%
    filter(qvalue <= q_cut, abs(meth.diff) >= diff_cut)

  if (nrow(dmr_filtered) == 0) {
    message("Skipping ", f, ": no DMRs passing thresholds")
    next
  }

  # Convert to GRanges
  gr <- GRanges(
    seqnames = dmr_filtered$chr,
    ranges = IRanges(start = dmr_filtered$start, end = dmr_filtered$end),
    meth.diff = dmr_filtered$meth.diff,
    qvalue = dmr_filtered$qvalue,
    pvalue = dmr_filtered$pvalue
  )

  # Merge overlapping or nearby DMRs
  merged_gr <- reduce(gr, min.gapwidth = max_gap, with.revmap = TRUE)

  # Summarize metadata for merged regions
  merged_df <- data.frame(
    chr = as.character(seqnames(merged_gr)),
    start = start(merged_gr),
    end = end(merged_gr),
    width = width(merged_gr),
    n_dmrs = lengths(merged_gr$revmap),
    stringsAsFactors = FALSE
  )

  # Calculate mean methylation difference for merged regions
  merged_df$mean_meth_diff <- sapply(merged_gr$revmap, function(idx) {
    mean(gr$meth.diff[idx])
  })

  merged_df$min_qvalue <- sapply(merged_gr$revmap, function(idx) {
    min(gr$qvalue[idx])
  })

  # Write output
  output_file <- file.path(
    outdir,
    paste0(contrast, "_Merged_DMRs_cov4_", direction, ".txt")
  )

  write_tsv(merged_df, output_file)
  
  message("Merged: ", f)
  message("  Before: ", nrow(dmr_filtered), " DMRs")
  message("  After: ", nrow(merged_df), " merged regions")
  message("  Output: ", output_file)
}

message("[complete] All DMR merging finished.")
