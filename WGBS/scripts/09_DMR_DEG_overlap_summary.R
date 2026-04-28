#!/usr/bin/env Rscript

# ==============================================================================
# Script: 09_DMR_DEG_overlap_summary.R
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

# 09_DMR_DEG_overlap_summary.R
#
# Summarize overlap between DMRs and differentially expressed genes (DEGs).
#
# Usage:
#   Rscript 09_DMR_DEG_overlap_summary.R \
#     --dmr-dir DMR_singleFeature \
#     --deg-dir DESeq2_results \
#     --outdir DMR_DEG_overlap

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

dmr_dir <- get_arg("--dmr-dir")
deg_dir <- get_arg("--deg-dir")
outdir <- get_arg("--outdir", "DMR_DEG_overlap")
padj_cut <- as.numeric(get_arg("--padj", "0.01"))
fc_cut <- as.numeric(get_arg("--fc", "1.3"))

if (is.null(dmr_dir) || is.null(deg_dir)) {
  stop("Required arguments: --dmr-dir and --deg-dir")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

log2fc_cut <- log2(fc_cut)

# Find DMR files
dmr_files <- list.files(
  dmr_dir,
  pattern = "_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

# Find DEG files
deg_files <- list.files(
  deg_dir,
  pattern = "\\.tsv$|\\.csv$",
  full.names = TRUE
)

if (length(dmr_files) == 0 || length(deg_files) == 0) {
  stop("No DMR or DEG files found")
}

overlap_results <- list()

for (dmr_file in dmr_files) {
  dmr_data <- fread(dmr_file)
  dmr_comparison <- unique(dmr_data$comparison)[1]
  dmr_class <- unique(dmr_data$class)[1]

  dmr_genes <- unique(dmr_data$gene_name[!is.na(dmr_data$gene_name)])
  dmr_genes <- unlist(strsplit(dmr_genes, ","))
  dmr_genes <- unique(dmr_genes[dmr_genes != ""])

  for (deg_file in deg_files) {
    deg_data <- fread(deg_file)

    if (!"padj" %in% names(deg_data) || !"log2FoldChange" %in% names(deg_data)) {
      next
    }

    deg_sig <- deg_data %>%
      filter(
        !is.na(padj),
        padj <= padj_cut,
        !is.na(log2FoldChange),
        abs(log2FoldChange) >= log2fc_cut
      )

    deg_genes <- unique(deg_sig$gene_name)

    overlap_genes <- intersect(dmr_genes, deg_genes)

    overlap_results[[paste(dmr_comparison, dmr_class, basename(deg_file))]] <- data.table(
      dmr_comparison = dmr_comparison,
      dmr_class = dmr_class,
      deg_file = basename(deg_file),
      n_dmr_genes = length(dmr_genes),
      n_deg_genes = length(deg_genes),
      n_overlap = length(overlap_genes),
      percent_dmr_overlap = 100 * length(overlap_genes) / max(1, length(dmr_genes)),
      overlap_genes = paste(overlap_genes, collapse = ",")
    )
  }
}

summary_table <- rbindlist(overlap_results)

fwrite(
  summary_table,
  file.path(outdir, "DMR_DEG_overlap_summary.tsv"),
  sep = "\t"
)

message("[done] DMR-DEG overlap summary written to: ", outdir)
