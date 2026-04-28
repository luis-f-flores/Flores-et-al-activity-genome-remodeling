#!/usr/bin/env Rscript

# ==============================================================================
# Script: 08_Global_CpG_methylation.R
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

# 08_Global_CpG_methylation.R
#
# Calculate global methylation levels from Bismark CpG reports.
#
# Usage:
#   Rscript 08_Global_CpG_methylation.R \
#     --reports-dir /path/to/CpG_reports \
#     --sample-sheet samples_wgbs.csv \
#     --outdir global_methylation

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

reports_dir <- get_arg("--reports-dir")
sample_sheet <- get_arg("--sample-sheet")
outdir <- get_arg("--outdir", "global_methylation")
mincov <- as.integer(get_arg("--mincov", "4"))

if (is.null(reports_dir) || is.null(sample_sheet)) {
  stop("Required arguments: --reports-dir and --sample-sheet")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

samples <- read_csv(sample_sheet, show_col_types = FALSE)

global_meth <- list()

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[i]
  report_file <- file.path(reports_dir, samples$file[i])

  if (!file.exists(report_file)) {
    warning("Report not found: ", report_file)
    next
  }

  # Read Bismark CpG report
  cpg <- fread(
    report_file,
    col.names = c("chr", "pos", "strand", "meth_count", "unmeth_count", "context", "trinuc")
  )

  # Filter by coverage
  cpg <- cpg[meth_count + unmeth_count >= mincov]

  # Calculate methylation percentage
  total_c <- sum(cpg$meth_count)
  total_cov <- sum(cpg$meth_count + cpg$unmeth_count)
  global_perc <- 100 * total_c / total_cov

  global_meth[[sample_id]] <- data.table(
    sample_id = sample_id,
    condition = if ("condition" %in% names(samples)) samples$condition[i] else NA_character_,
    total_cpgs = nrow(cpg),
    total_methylated_calls = total_c,
    total_coverage = total_cov,
    global_methylation_percent = global_perc
  )

  message("Processed: ", sample_id, " (", round(global_perc, 2), "%)")
}

# Combine results
global_summary <- rbindlist(global_meth)

fwrite(
  global_summary,
  file.path(outdir, "global_CpG_methylation.tsv"),
  sep = "\t"
)

# Plot if condition info available
if ("condition" %in% names(global_summary) && !all(is.na(global_summary$condition))) {
  p <- ggplot(global_summary, aes(x = condition, y = global_methylation_percent, fill = condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    labs(
      title = "Global CpG Methylation by Condition",
      x = "Condition",
      y = "Global Methylation (%)"
    ) +
    theme_bw() +
    theme(legend.position = "none")

  ggsave(
    file.path(outdir, "global_methylation_boxplot.pdf"),
    p,
    width = 6,
    height = 5
  )
}

message("[done] Global methylation summary written to: ", outdir)
