#!/usr/bin/env Rscript

# ==============================================================================
# Script: 14_DMR_overlap_upset_plot.R
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

# 14_DMR_overlap_upset_plot.R
#
# Create UpSet plot showing overlap patterns across DMR sets.
#
# Usage:
#   Rscript 14_DMR_overlap_upset_plot.R \
#     --input-dir DMG_promoter_meth25_q01 \
#     --outdir DMG_upset_plots

suppressPackageStartupMessages({
  library(UpSetR)
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
outdir <- get_arg("--outdir", "DMG_upset_plots")

if (is.null(input_dir)) {
  stop("Usage: Rscript 14_DMR_overlap_upset_plot.R --input-dir DMG_promoter_meth25_q01")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all DMG files
dmg_files <- list.files(
  input_dir,
  pattern = "_promoter_DMG_meth25_q01\\.tsv$",
  full.names = TRUE
)

if (length(dmg_files) < 2) {
  stop("Need at least 2 DMG files for UpSet plot")
}

# Load gene sets
gene_sets <- list()

for (f in dmg_files) {
  dt <- fread(f)
  label <- sub("_promoter_DMG_meth25_q01\\.tsv$", "", basename(f))
  gene_sets[[label]] <- unique(dt$gene_name[!is.na(dt$gene_name)])
}

# Create binary matrix for UpSet
all_genes <- unique(unlist(gene_sets))
upset_data <- data.frame(gene = all_genes)

for (set_name in names(gene_sets)) {
  upset_data[[set_name]] <- as.integer(upset_data$gene %in% gene_sets[[set_name]])
}

# Generate UpSet plot
pdf(file.path(outdir, "DMG_upset_plot.pdf"), width = 12, height = 8)
upset(
  upset_data,
  sets = names(gene_sets),
  keep.order = FALSE,
  order.by = "freq"
)
dev.off()

message("[done] UpSet plot saved to: ", file.path(outdir, "DMG_upset_plot.pdf"))
