#!/usr/bin/env Rscript


# ==============================================================================
# Script: 20_generate_normalized_bedgraph_from_methylKit.R
# Package: WGBS / methylation analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for WGBS DMC/DMR/DMG analysis, annotation, overlap analysis, methylation-expression integration, signal tracks, and OLS modeling.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.

# ==============================================================================

# 20_generate_normalized_bedgraph_from_methylKit.R
#
# Generate normalized WGBS methylation bedGraph files from Bismark CpG reports.
#
# Usage:
#   Rscript 20_generate_normalized_bedgraph_from_methylKit.R \
#     --reports-dir /path/to/CpG_reports \
#     --sample-sheet sample_sheet_wgbs.csv \
#     --outdir methylKit_bedgraph_bigwig \
#     --mincov 4 \
#     --hi-perc 99.9
#
# sample_sheet_wgbs.csv must contain:
#   sample_id,condition,file,treatment

suppressPackageStartupMessages({
  
library(methylKit)
  
library(readr)

})


# -----------------------------------------------------------------------------
# Parse command-line arguments
# -----------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

reports_dir <- get_arg("--reports-dir", ".")
sample_sheet <- get_arg("--sample-sheet")
outdir <- get_arg("--outdir", "methylKit_bedgraph_bigwig")
mincov <- as.integer(get_arg("--mincov", "4"))
hi_perc <- as.numeric(get_arg("--hi-perc", "99.9"))

if (is.null(sample_sheet)) {
  stop("Usage: Rscript 20_generate_normalized_bedgraph_from_methylKit.R --sample-sheet sample_sheet_wgbs.csv")
}
if (!file.exists(sample_sheet)) stop("sample sheet not found: ", sample_sheet)


# -----------------------------------------------------------------------------
# Create output directories / initialize outputs
# -----------------------------------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- read_csv(sample_sheet, show_col_types = FALSE)
required <- c("sample_id", "file")
missing <- setdiff(required, names(samples))
if (length(missing) > 0) stop("sample sheet missing columns: ", paste(missing, collapse = ", "))

files <- file.path(reports_dir, samples$file)
if (!all(file.exists(files))) {
  stop("Missing CpG report files: ", paste(files[!file.exists(files)], collapse = ", "))
}

treatment <- if ("treatment" %in% names(samples)) samples$treatment else rep(0, nrow(samples))

myobj <- methRead(
  as.list(files),
  sample.id = as.list(samples$sample_id),
  assembly = "mm10",
  treatment = treatment,
  pipeline = "bismarkCytosineReport",
  mincov = mincov
)

myobj <- filterByCoverage(myobj, hi.perc = hi_perc)
myobj <- normalizeCoverage(myobj)


# -----------------------------------------------------------------------------
# Main analysis loop
# -----------------------------------------------------------------------------

for (i in seq_along(myobj)) {
  out_file <- file.path(outdir, paste0(samples$sample_id[i], "_normalized_cov", mincov, ".bedGraph"))
  bedgraph(myobj[[i]], file.name = out_file, col.name = "perc.meth")
  message("Wrote: ", out_file)
}

writeLines(
  c(
    paste0("sample_sheet=", normalizePath(sample_sheet)),
    paste0("reports_dir=", normalizePath(reports_dir)),
    paste0("mincov=", mincov),
    paste0("hi_perc=", hi_perc),
    "normalization=methylKit::normalizeCoverage",
    "bedGraph_value=percent_methylation"
  ),
  file.path(outdir, "bedgraph_generation_parameters.txt")
)
