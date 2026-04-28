#!/usr/bin/env Rscript

# ==============================================================================
# Script: 13_DMR_overlap_any.R
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

# 13_DMR_overlap_any.R
#
# Compute pairwise overlaps between DMR sets using GenomicRanges.
#
# Usage:
#   Rscript 13_DMR_overlap_any.R \
#     --input-dir DMR_singleFeature \
#     --outdir DMR_pairwise_overlaps

suppressPackageStartupMessages({
  library(GenomicRanges)
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
outdir <- get_arg("--outdir", "DMR_pairwise_overlaps")
min_overlap <- as.integer(get_arg("--min-overlap", "1"))

if (is.null(input_dir)) {
  stop("Usage: Rscript 13_DMR_overlap_any.R --input-dir DMR_singleFeature")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all DMR files
dmr_files <- list.files(
  input_dir,
  pattern = "_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

if (length(dmr_files) < 2) {
  stop("Need at least 2 DMR files for pairwise comparison")
}

# Load all DMR sets as GRanges
dmr_sets <- list()

for (f in dmr_files) {
  dt <- fread(f)
  label <- sub("_singleFeature_q05_meth5\\.tsv$", "", basename(f))

  dmr_sets[[label]] <- GRanges(
    seqnames = dt$chr,
    ranges = IRanges(start = dt$start, end = dt$end)
  )
}

# Compute pairwise overlaps
overlap_matrix <- matrix(
  0,
  nrow = length(dmr_sets),
  ncol = length(dmr_sets),
  dimnames = list(names(dmr_sets), names(dmr_sets))
)

for (i in seq_along(dmr_sets)) {
  for (j in seq_along(dmr_sets)) {
    overlaps <- findOverlaps(
      dmr_sets[[i]],
      dmr_sets[[j]],
      minoverlap = min_overlap
    )

    overlap_matrix[i, j] <- length(unique(queryHits(overlaps)))
  }
}

# Write overlap matrix
fwrite(
  as.data.table(overlap_matrix, keep.rownames = "DMR_set"),
  file.path(outdir, "DMR_pairwise_overlap_matrix.tsv"),
  sep = "\t"
)

message("[done] Pairwise DMR overlaps written to: ", outdir)
