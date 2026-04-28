#!/usr/bin/env Rscript

# ==============================================================================
# Script: 00_methylKit_DMC_calling.R
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

# 00_methylKit_DMC_calling.R
#
# Identify differentially methylated cytosines (DMCs) using methylKit.
#
# Usage:
#   Rscript 00_methylKit_DMC_calling.R \
#     --reports-dir /path/to/bismark_CpG_reports \
#     --sample-sheet samples_wgbs.csv \
#     --contrast-sheet contrasts_wgbs.csv \
#     --outdir DMC_outputs \
#     --mincov 4 \
#     --hi-perc 99.9 \
#     --diff 5 \
#     --pvalue 0.001 \
#     --qvalue 0.1 \
#     --cores 4

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
outdir <- get_arg("--outdir", "DMC_outputs")
mincov <- as.integer(get_arg("--mincov", "4"))
hi_perc <- as.numeric(get_arg("--hi-perc", "99.9"))
diff_cut <- as.numeric(get_arg("--diff", "5"))
p_cut <- as.numeric(get_arg("--pvalue", "0.001"))
q_cut <- as.numeric(get_arg("--qvalue", "0.1"))
cores <- as.integer(get_arg("--cores", "1"))

if (is.null(sample_sheet) || is.null(contrast_sheet)) {
  stop("Required arguments: --sample-sheet and --contrast-sheet")
}

# -----------------------------------------------------------------------------
# Create output directories / initialize outputs
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

  contrast_dir <- file.path(outdir, paste0("DMC_", comp, "_analysis"))
  dir.create(contrast_dir, showWarnings = FALSE, recursive = TRUE)

  # Read methylation data
  myobj <- methRead(
    as.list(files),
    sample.id = as.list(ss$sample_id),
    assembly = "mm10",
    treatment = treatment,
    pipeline = "bismarkCytosineReport",
    mincov = mincov
  )

  # Filter and normalize
  myobj <- filterByCoverage(myobj, hi.perc = hi_perc)
  myobj <- normalizeCoverage(myobj)

  # Unite samples
  meth <- unite(myobj, destrand = FALSE, min.per.group = 2L)

  # Calculate differential methylation
  dmcs <- calculateDiffMeth(
    meth,
    test = "Chisq",
    adjust = "SLIM",
    effect = "wmean",
    mc.cores = cores,
    save.db = FALSE
  )

  all_data <- getData(dmcs)

  # -----------------------------------------------------------------------------
  # Write outputs
  # -----------------------------------------------------------------------------

  # Write all DMCs
  write.table(
    all_data,
    file.path(contrast_dir, paste0("all_DMCs_", comp, ".txt")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Write q-value filtered DMCs
  q_data <- getData(getMethylDiff(dmcs, difference = diff_cut, qvalue = q_cut))
  write.table(
    q_data,
    file.path(
      contrast_dir,
      paste0("significant_DMCs_", comp, "_q", q_cut, "_diff", diff_cut, ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Write manual p-value filtered DMCs
  manual_all <- subset(all_data, abs(meth.diff) >= diff_cut & pvalue <= p_cut)
  manual_hyper <- subset(all_data, meth.diff >= diff_cut & pvalue <= p_cut)
  manual_hypo <- subset(all_data, meth.diff <= -diff_cut & pvalue <= p_cut)

  write.table(
    manual_all,
    file.path(contrast_dir, "manual_filtered_DMCs_p0001_diff5_all.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    manual_hyper,
    file.path(contrast_dir, "manual_filtered_DMCs_p0001_diff5_hyper.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  write.table(
    manual_hypo,
    file.path(contrast_dir, "manual_filtered_DMCs_p0001_diff5_hypo.txt"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Write summary counts
  write.csv(
    data.frame(
      contrast = comp,
      threshold = c(
        paste0("p<=", p_cut, " & |delta|>=", diff_cut),
        paste0("q<=", q_cut, " & |delta|>=", diff_cut)
      ),
      hypomethylated = c(
        nrow(manual_hypo),
        sum(q_data$meth.diff <= -diff_cut)
      ),
      hypermethylated = c(
        nrow(manual_hyper),
        sum(q_data$meth.diff >= diff_cut)
      )
    ),
    file.path(contrast_dir, "DMC_direction_summary_counts.csv"),
    row.names = FALSE
  )

  message("[done] Processed contrast: ", comp)
}

message("[complete] All DMC analyses finished.")
