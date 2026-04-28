#!/usr/bin/env Rscript

# ==============================================================================
# Script: 16_Extract_log2fc_FPKM_DMG.R
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

# 16_Extract_log2fc_FPKM_DMG.R
#
# Extract expression data (log2FC FPKM) for DMGs.
#
# Usage:
#   Rscript 16_Extract_log2fc_FPKM_DMG.R \
#     --dmg-dir DMG_promoter_meth25_q01 \
#     --fpkm-file combined_log2FC_FPKM_GlobalDEG.csv \
#     --outdir DMG_with_expression

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

dmg_dir <- get_arg("--dmg-dir")
fpkm_file <- get_arg("--fpkm-file")
outdir <- get_arg("--outdir", "DMG_with_expression")

if (is.null(dmg_dir) || is.null(fpkm_file)) {
  stop("Required arguments: --dmg-dir and --fpkm-file")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load expression data
fpkm_data <- fread(fpkm_file)

if (!"gene_name" %in% names(fpkm_data)) {
  if (ncol(fpkm_data) > 0) {
    setnames(fpkm_data, 1, "gene_name")
  }
}

# Find all DMG files
dmg_files <- list.files(
  dmg_dir,
  pattern = "_promoter_DMG_meth25_q01\\.tsv$",
  full.names = TRUE
)

if (length(dmg_files) == 0) {
  stop("No DMG files found in: ", dmg_dir)
}

for (f in dmg_files) {
  dmg_data <- fread(f)

  # Extract timepoint from filename
  timepoint <- sub(".*_(\\d+hr)_.*", "\\1", basename(f))
  fpkm_col <- paste0("log2FC_", timepoint)

  if (!fpkm_col %in% names(fpkm_data)) {
    warning("FPKM column not found: ", fpkm_col)
    next
  }

  # Merge DMG data with expression data
  merged <- merge(
    dmg_data,
    fpkm_data[, c("gene_name", fpkm_col), with = FALSE],
    by = "gene_name",
    all.x = TRUE
  )

  setnames(merged, fpkm_col, "log2FC_fpkm")

  # Write output
  output_file <- file.path(
    outdir,
    sub(
      "_promoter_DMG_meth25_q01\\.tsv$",
      "_with_FPKM.tsv",
      basename(f)
    )
  )

  fwrite(merged, output_file, sep = "\t")

  message("Merged: ", basename(f))
  message("  Genes with expression data: ", sum(!is.na(merged$log2FC_fpkm)))
}

message("[done] Expression data extraction complete")
