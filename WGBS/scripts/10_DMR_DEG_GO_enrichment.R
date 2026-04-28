#!/usr/bin/env Rscript

# ==============================================================================
# Script: 10_DMR_DEG_GO_enrichment.R
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

# 10_DMR_DEG_GO_enrichment.R
#
# Perform GO enrichment analysis on genes with both DMRs and differential expression.
#
# Usage:
#   Rscript 10_DMR_DEG_GO_enrichment.R \
#     --overlap-file DMR_DEG_overlap_summary.tsv \
#     --outdir GO_enrichment_DMR_DEG

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
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

overlap_file <- get_arg("--overlap-file")
outdir <- get_arg("--outdir", "GO_enrichment_DMR_DEG")
ont <- toupper(get_arg("--ont", "BP"))
p_cut <- as.numeric(get_arg("--pvalue", "0.05"))

if (is.null(overlap_file)) {
  stop("Usage: Rscript 10_DMR_DEG_GO_enrichment.R --overlap-file DMR_DEG_overlap_summary.tsv")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

overlap_data <- fread(overlap_file)

for (i in seq_len(nrow(overlap_data))) {
  row <- overlap_data[i, ]

  if (row$n_overlap < 5) {
    message("Skipping: ", row$dmr_comparison, " (", row$dmr_class, ") - too few genes (", row$n_overlap, ")")
    next
  }

  genes <- unlist(strsplit(row$overlap_genes, ","))
  genes <- genes[genes != ""]

  if (length(genes) < 5) next

  # Convert to Entrez IDs
  gene_ids <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Mm.eg.db
  )

  if (nrow(gene_ids) < 5) {
    message("Skipping: too few mapped genes")
    next
  }

  # GO enrichment
  go_result <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Mm.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cut,
    readable = TRUE
  )

  if (is.null(go_result) || nrow(go_result@result) == 0) {
    message("No significant GO terms for: ", row$dmr_comparison, " (", row$dmr_class, ")")
    next
  }

  output_prefix <- paste0(row$dmr_comparison, "_", row$dmr_class, "_", sub("\\.tsv$|\\.csv$", "", row$deg_file))

  fwrite(
    as.data.table(go_result@result),
    file.path(outdir, paste0(output_prefix, "_GO_enrichment.tsv")),
    sep = "\t"
  )

  p <- dotplot(go_result, showCategory = 15)
  ggsave(
    file.path(outdir, paste0(output_prefix, "_GO_dotplot.pdf")),
    p,
    width = 8,
    height = 6
  )

  message("GO enrichment: ", output_prefix)
}

message("[done] GO enrichment analyses complete")
