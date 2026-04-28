#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_DMR_feature_methylation_summary.R
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

# 06_DMR_feature_methylation_summary.R
#
# Summarize methylation changes by genomic feature category.
#
# Usage:
#   Rscript 06_DMR_feature_methylation_summary.R \
#     --input-dir DMR_singleFeature \
#     --outdir DMR_summaries

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

input_dir <- get_arg("--input-dir")
outdir <- get_arg("--outdir", input_dir)

if (is.null(input_dir)) {
  stop("Usage: Rscript 06_DMR_feature_methylation_summary.R --input-dir DMR_singleFeature")
}

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all single-feature DMR files
files <- list.files(
  input_dir,
  pattern = "_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No single-feature DMR files found in: ", input_dir)
}

# Combine all data
dt <- rbindlist(lapply(files, function(f) {
  x <- fread(f)
  x[, source_file := basename(f)]
  x[, abs_meth_diff := abs(meth_diff)]
  x
}), fill = TRUE)

# Summary by feature category
summ <- dt[, .(
  n_dmrs = .N,
  mean_meth_diff = mean(abs_meth_diff, na.rm = TRUE),
  median_meth_diff = median(abs_meth_diff, na.rm = TRUE),
  mean_width = mean(end - start, na.rm = TRUE)
), by = .(comparison, class, feature_category)]

# Write summary table
fwrite(
  summ,
  file.path(outdir, "DMR_feature_methylation_summary.tsv"),
  sep = "\t"
)

# Write Excel workbook with multiple sheets
wb <- createWorkbook()

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summ)

addWorksheet(wb, "AllDMRs")
writeData(wb, "AllDMRs", dt)

saveWorkbook(
  wb,
  file.path(outdir, "DMR_feature_summary.xlsx"),
  overwrite = TRUE
)

message("[done] Feature methylation summaries written to: ", outdir)
