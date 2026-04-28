#!/usr/bin/env Rscript

# ==============================================================================
# Script: 07_DMR_CGI_methylation_summary.R
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

# 07_DMR_CGI_methylation_summary.R
#
# Summarize DMR overlap with CpG islands.
#
# Usage:
#   Rscript 07_DMR_CGI_methylation_summary.R \
#     --dmr-dir DMR_singleFeature \
#     --cgi-bed mm10_cpg_islands.bed \
#     --outdir DMR_CGI_summaries

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

dmr_dir <- get_arg("--dmr-dir")
cgi_bed <- get_arg("--cgi-bed")
outdir <- get_arg("--outdir", "DMR_CGI_summaries")

if (is.null(dmr_dir) || is.null(cgi_bed)) {
  stop("Required arguments: --dmr-dir and --cgi-bed")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Load CpG islands
cgi <- fread(cgi_bed)
setnames(cgi, c("chr", "start", "end"), c("chr", "start", "end")[1:min(3, ncol(cgi))])

cgi_gr <- GRanges(
  seqnames = cgi$chr,
  ranges = IRanges(start = cgi$start, end = cgi$end)
)

# Find all DMR files
dmr_files <- list.files(
  dmr_dir,
  pattern = "_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

if (length(dmr_files) == 0) {
  stop("No DMR files found in: ", dmr_dir)
}

results <- list()

for (f in dmr_files) {
  dmr_data <- fread(f)
  comparison <- unique(dmr_data$comparison)[1]
  class <- unique(dmr_data$class)[1]

  dmr_gr <- GRanges(
    seqnames = dmr_data$chr,
    ranges = IRanges(start = dmr_data$start, end = dmr_data$end)
  )

  overlaps <- findOverlaps(dmr_gr, cgi_gr)
  n_overlapping <- length(unique(queryHits(overlaps)))

  results[[basename(f)]] <- data.table(
    comparison = comparison,
    class = class,
    total_dmrs = length(dmr_gr),
    dmrs_overlapping_cgi = n_overlapping,
    percent_overlapping = 100 * n_overlapping / length(dmr_gr)
  )

  message("Processed: ", comparison, " (", class, ")")
}

summary_table <- rbindlist(results)
fwrite(
  summary_table,
  file.path(outdir, "DMR_CGI_overlap_summary.tsv"),
  sep = "\t"
)

message("[done] CGI overlap summaries written to: ", outdir)
