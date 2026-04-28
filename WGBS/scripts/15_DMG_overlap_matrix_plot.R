#!/usr/bin/env Rscript

# ==============================================================================
# Script: 15_DMG_overlap_matrix_plot.R
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

# 15_DMG_overlap_matrix_plot.R
#
# Create heatmap of DMG overlap matrix.
#
# Usage:
#   Rscript 15_DMG_overlap_matrix_plot.R \
#     --overlap-matrix DMR_pairwise_overlap_matrix.tsv \
#     --outdir DMG_overlap_plots

suppressPackageStartupMessages({
  library(pheatmap)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

overlap_matrix_file <- get_arg("--overlap-matrix")
outdir <- get_arg("--outdir", "DMG_overlap_plots")

if (is.null(overlap_matrix_file)) {
  stop("Usage: Rscript 15_DMG_overlap_matrix_plot.R --overlap-matrix matrix.tsv")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load overlap matrix
overlap_data <- fread(overlap_matrix_file)

# Convert to matrix
mat <- as.matrix(overlap_data[, -1])
rownames(mat) <- overlap_data[[1]]

# Create heatmap
pdf(file.path(outdir, "DMG_overlap_heatmap.pdf"), width = 10, height = 10)
pheatmap(
  mat,
  display_numbers = TRUE,
  number_color = "black",
  fontsize_number = 8,
  color = colorRampPalette(c("white", "blue"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "DMG Pairwise Overlap"
)
dev.off()

message("[done] Overlap heatmap saved to: ", file.path(outdir, "DMG_overlap_heatmap.pdf"))
