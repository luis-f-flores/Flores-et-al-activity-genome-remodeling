#!/usr/bin/env Rscript

# ==============================================================================
# Script: 17_DMG_methylation_expression_correlation.R
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

# 17_DMG_methylation_expression_correlation.R
#
# Compute correlation between promoter methylation changes and gene expression changes.
#
# Usage:
#   Rscript 17_DMG_methylation_expression_correlation.R \
#     --input-dir DMG_with_expression \
#     --outdir DMG_correlation_results

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

input_dir <- get_arg("--input-dir")
outdir <- get_arg("--outdir", "DMG_correlation_results")

if (is.null(input_dir)) {
  stop("Usage: Rscript 17_DMG_methylation_expression_correlation.R --input-dir DMG_with_expression")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Find all DMG+FPKM files
dmg_files <- list.files(
  input_dir,
  pattern = "_with_FPKM\\.tsv$",
  full.names = TRUE
)

if (length(dmg_files) == 0) {
  stop("No DMG+FPKM files found in: ", input_dir)
}

correlation_results <- list()

for (f in dmg_files) {
  dt <- fread(f)

  # Filter to genes with both methylation and expression data
  complete_data <- dt[!is.na(meth_diff) & !is.na(log2FC_fpkm)]

  if (nrow(complete_data) < 3) {
    message("Skipping: ", basename(f), " - insufficient data")
    next
  }

  # Compute correlation
  cor_test <- cor.test(
    complete_data$meth_diff,
    complete_data$log2FC_fpkm,
    method = "pearson"
  )

  label <- sub("_with_FPKM\\.tsv$", "", basename(f))

  correlation_results[[label]] <- data.table(
    comparison = label,
    n_genes = nrow(complete_data),
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    conf_low = cor_test$conf.int[1],
    conf_high = cor_test$conf.int[2]
  )

  # Create scatter plot
  p <- ggplot(complete_data, aes(x = meth_diff, y = log2FC_fpkm)) +
    geom_point(alpha = 0.5, size = 2) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    labs(
      title = paste0("Methylation-Expression Correlation: ", label),
      subtitle = sprintf(
        "r = %.3f, p = %.2e, n = %d",
        cor_test$estimate,
        cor_test$p.value,
        nrow(complete_data)
      ),
      x = "Promoter Methylation Difference (%)",
      y = "Gene Expression Change (log2FC FPKM)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )

  ggsave(
    file.path(outdir, paste0(label, "_correlation_plot.pdf")),
    p,
    width = 8,
    height = 6
  )

  message("Analyzed: ", label)
  message("  Correlation: r = ", round(cor_test$estimate, 3), ", p = ", format(cor_test$p.value, scientific = TRUE))
}

# Write summary table
if (length(correlation_results) > 0) {
  summary_table <- rbindlist(correlation_results)

  fwrite(
    summary_table,
    file.path(outdir, "methylation_expression_correlation_summary.tsv"),
    sep = "\t"
  )

  message("[done] Correlation analyses complete")
  message("  Summary written to: ", file.path(outdir, "methylation_expression_correlation_summary.tsv"))
}
