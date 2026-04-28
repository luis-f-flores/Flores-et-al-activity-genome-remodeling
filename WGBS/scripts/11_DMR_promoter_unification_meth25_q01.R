#!/usr/bin/env Rscript

# ==============================================================================
# Script: 11_DMR_promoter_unification_meth25_q01.R
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

# 11_DMR_promoter_unification_meth25_q01.R
#
# Extract promoter DMRs with stringent thresholds (q <= 0.01, |delta meth| >= 25%).
#
# Usage:
#   Rscript 11_DMR_promoter_unification_meth25_q01.R \
#     --input-dir DMR_singleFeature \
#     --outdir DMG_promoter_meth25_q01

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

input_dir <- get_arg("--input-dir")
outdir <- get_arg("--outdir", "DMG_promoter_meth25_q01")
q_cut <- as.numeric(get_arg("--qvalue", "0.01"))
meth_cut <- as.numeric(get_arg("--meth-diff", "25"))

if (is.null(input_dir)) {
  stop("Usage: Rscript 11_DMR_promoter_unification_meth25_q01.R --input-dir DMR_singleFeature")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all single-feature DMR files
dmr_files <- list.files(
  input_dir,
  pattern = "_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

if (length(dmr_files) == 0) {
  stop("No DMR files found in: ", input_dir)
}

for (f in dmr_files) {
  dmr_data <- fread(f)

  # Filter for promoter DMRs with stringent thresholds
  promoter_dmr <- dmr_data[
    feature_category == "Promoter" &
    qvalue <= q_cut &
    abs(meth_diff) >= meth_cut
  ]

  if (nrow(promoter_dmr) == 0) {
    message("No promoter DMRs passing thresholds in: ", basename(f))
    next
  }

  # One row per gene (if gene has multiple DMRs, take the one with largest |meth_diff|)
  promoter_dmr[, abs_meth_diff := abs(meth_diff)]
  unified <- promoter_dmr[order(-abs_meth_diff), .SD[1], by = gene_name]

  output_file <- file.path(
    outdir,
    sub("_singleFeature_q05_meth5\\.tsv$", "_promoter_DMG_meth25_q01.tsv", basename(f))
  )

  fwrite(unified, output_file, sep = "\t")

  message("Extracted: ", basename(f))
  message("  Promoter DMGs: ", nrow(unified))
}

message("[done] Promoter DMG extraction complete")
