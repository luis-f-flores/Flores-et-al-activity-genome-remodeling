#!/usr/bin/env Rscript

# ==============================================================================
# Script: 12_DMG_overlap_summary.R
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

# 12_DMG_overlap_summary.R
#
# Summarize overlap of differentially methylated genes (DMGs) across timepoints.
#
# Usage:
#   Rscript 12_DMG_overlap_summary.R \
#     --input-dir DMG_promoter_meth25_q01 \
#     --outdir DMG_overlap_summaries

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
outdir <- get_arg("--outdir", "DMG_overlap_summaries")

if (is.null(input_dir)) {
  stop("Usage: Rscript 12_DMG_overlap_summary.R --input-dir DMG_promoter_meth25_q01")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all promoter DMG files
dmg_files <- list.files(
  input_dir,
  pattern = "_promoter_DMG_meth25_q01\\.tsv$",
  full.names = TRUE
)

if (length(dmg_files) == 0) {
  stop("No DMG files found in: ", input_dir)
}

# Load all DMG data
dmg_list <- list()

for (f in dmg_files) {
  dt <- fread(f)
  dt[, source_file := basename(f)]
  dmg_list[[basename(f)]] <- dt
}

all_dmg <- rbindlist(dmg_list, fill = TRUE)

# Summarize overlap across comparisons
gene_counts <- all_dmg[, .(
  n_comparisons = .N,
  comparisons = paste(unique(comparison), collapse = ","),
  classes = paste(unique(class), collapse = ","),
  mean_meth_diff = mean(abs(meth_diff))
), by = gene_name]

# Write summary
fwrite(
  gene_counts[order(-n_comparisons)],
  file.path(outdir, "DMG_gene_overlap_summary.tsv"),
  sep = "\t"
)

# Shared vs unique genes
shared_genes <- gene_counts[n_comparisons > 1, gene_name]
unique_genes <- gene_counts[n_comparisons == 1, gene_name]

fwrite(
  data.table(
    category = c("shared", "unique", "total"),
    n_genes = c(length(shared_genes), length(unique_genes), nrow(gene_counts))
  ),
  file.path(outdir, "DMG_sharing_summary.tsv"),
  sep = "\t"
)

message("[done] DMG overlap summaries written to: ", outdir)
message("  Total unique DMGs: ", nrow(gene_counts))
message("  Shared across >1 comparison: ", length(shared_genes))
