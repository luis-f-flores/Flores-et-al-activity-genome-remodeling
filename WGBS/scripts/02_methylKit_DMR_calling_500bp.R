#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_methylKit_DMR_calling_500bp.R
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

# 02_methylKit_DMR_calling_500bp.R
#
# Identify differentially methylated regions (DMRs) using methylKit with
# sliding windows.
#
# Usage:
#   Rscript 02_methylKit_DMR_calling_500bp.R \
#     --reports-dir /path/to/bismark_CpG_reports \
#     --sample-sheet samples_wgbs.csv \
#     --contrast-sheet contrasts_wgbs.csv \
#     --outdir methylKit_DMR \
#     --win 500 \
#     --step 250 \
#     --cov-bases 3 \
#     --mincov 4 \
#     --qvalue 0.05 \
#     --diff 5

suppressPackageStartupMessages({
  library(methylKit)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

reports_dir <- get_arg("--reports-dir", ".")
sample_sheet <- get_arg("--sample-sheet")
contrast_sheet <- get_arg("--contrast-sheet")
outdir <- get_arg("--outdir", "methylKit_DMR")
win <- as.integer(get_arg("--win", "500"))
step <- as.integer(get_arg("--step", "250"))
cov_bases <- as.integer(get_arg("--cov-bases", "3"))
mincov <- as.integer(get_arg("--mincov", "4"))
q_cut <- as.numeric(get_arg("--qvalue", "0.05"))
diff_cut <- as.numeric(get_arg("--diff", "5"))

if (is.null(sample_sheet) || is.null(contrast_sheet)) {
  stop("Required arguments: --sample-sheet and --contrast-sheet")
}

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- read_csv(sample_sheet, show_col_types = FALSE)
contrasts <- read_csv(contrast_sheet, show_col_types = FALSE)

# -----------------------------------------------------------------------------
# Main analysis loop: process each contrast
# -----------------------------------------------------------------------------

for (r in seq_len(nrow(contrasts))) {
  comp <- contrasts$contrast[r]
  ctrl <- contrasts$control_condition[r]
  trt <- contrasts$treatment_condition[r]

  # Filter samples for this contrast
  ss <- samples %>%
    filter(condition %in% c(ctrl, trt)) %>%
    arrange(condition != ctrl)

  files <- file.path(reports_dir, ss$file)
  treatment <- ifelse(ss$condition == trt, 1, 0)

  contrast_dir <- file.path(outdir, paste0("DMR_", comp, "_win", win))
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

  # Read methylation data with windowing
  myobj <- methRead(
    as.list(files),
    sample.id = as.list(ss$sample_id),
    assembly = "mm10",
    treatment = treatment,
    pipeline = "bismarkCytosineReport",
    mincov = mincov
  )

  # Filter by coverage
  myobj <- filterByCoverage(myobj, lo.count = NULL, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)

  # Tile genome into windows
  tiles <- tileMethylCounts(myobj, win.size = win, step.size = step, cov.bases = cov_bases)

  # Normalize
  tiles <- normalizeCoverage(tiles)

  # Unite samples
  meth <- unite(tiles, destrand = FALSE, min.per.group = 2L)

  # Calculate differential methylation for regions
  dmrs <- calculateDiffMeth(meth, test = "Chisq", adjust = "SLIM", mc.cores = 1)

  all_data <- getData(dmrs)

  # -----------------------------------------------------------------------------
  # Write outputs
  # -----------------------------------------------------------------------------

  # Write all DMRs
  write.table(
    all_data,
    file.path(contrast_dir, paste0("all_DMRs_", comp, "_win", win, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Write filtered DMRs
  hyper_dmr <- getMethylDiff(dmrs, difference = diff_cut, qvalue = q_cut, type = "hyper")
  hypo_dmr <- getMethylDiff(dmrs, difference = diff_cut, qvalue = q_cut, type = "hypo")

  if (nrow(getData(hyper_dmr)) > 0) {
    write.table(
      getData(hyper_dmr),
      file.path(contrast_dir, paste0("DMRs_", comp, "_win", win, "_hyper_q", q_cut, "_diff", diff_cut, ".txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  if (nrow(getData(hypo_dmr)) > 0) {
    write.table(
      getData(hypo_dmr),
      file.path(contrast_dir, paste0("DMRs_", comp, "_win", win, "_hypo_q", q_cut, "_diff", diff_cut, ".txt")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  message("[done] Processed DMR contrast: ", comp)
  message("  Hypermethylated DMRs: ", nrow(getData(hyper_dmr)))
  message("  Hypomethylated DMRs: ", nrow(getData(hypo_dmr)))
}

message("[complete] All DMR analyses finished.")
