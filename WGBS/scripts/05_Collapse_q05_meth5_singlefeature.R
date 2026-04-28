#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_Collapse_q05_meth5_singlefeature.R
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

# 05_Collapse_q05_meth5_singlefeature.R
#
# Collapse annotated DMRs to a single genomic feature per region based on
# priority ordering.
#
# Usage:
#   Rscript 05_Collapse_q05_meth5_singlefeature.R \
#     --input-dir DMR_merged \
#     --outdir DMR_singleFeature

suppressPackageStartupMessages({
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
outdir <- get_arg("--outdir", file.path(input_dir, "SingleFeature_q05_meth5"))

if (is.null(input_dir)) {
  stop("Usage: Rscript 05_Collapse_q05_meth5_singlefeature.R --input-dir DMR_merged")
}

# -----------------------------------------------------------------------------
# Create output directories
# -----------------------------------------------------------------------------

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Feature priority (first match wins)
priority <- c("Promoter", "Exon", "Intron", "Intergenic", "5'UTR", "3'UTR")

# Map feature names to categories
map_feature <- function(x) {
  x <- tolower(x)
  ifelse(grepl("promoter|tss", x), "Promoter",
  ifelse(grepl("exon", x), "Exon",
  ifelse(grepl("intron", x), "Intron",
  ifelse(grepl("intergenic", x), "Intergenic",
  ifelse(grepl("fiveutr|5utr", x), "5'UTR",
  ifelse(grepl("threeutr|3utr", x), "3'UTR", NA_character_))))))
}

# Find all merged DMR files
files <- list.files(
  input_dir,
  pattern = "_Merged_DMRs_cov4_(hyper|hypo)\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No DMR files found in: ", input_dir)
}

summary_list <- list()

# -----------------------------------------------------------------------------
# Process each DMR file
# -----------------------------------------------------------------------------

for (raw_file in files) {
  class <- ifelse(grepl("hyper", basename(raw_file)), "hyper", "hypo")
  comparison <- sub("_Merged_DMRs_cov4_(hyper|hypo)\\.txt$", "", basename(raw_file))

  # Find corresponding annotation file
  ann_file <- file.path(
    dirname(raw_file),
    paste0("GTF_Annot_", basename(raw_file))
  )

  if (!file.exists(ann_file)) {
    warning("Annotation file not found: ", ann_file)
    next
  }

  # Load DMRs and annotations
  dmrs <- fread(raw_file)

  if ("meth.diff" %in% names(dmrs)) {
    setnames(dmrs, "meth.diff", "meth_diff")
  }

  dmrs[, dmr_id := paste(chr, start, end, sep = "_")]

  ann <- fread(ann_file)
  ann[, dmr_id := paste(chr, start, end, sep = "_")]
  ann[, feature_category := map_feature(feature)]

  # Collapse to single feature per DMR using priority
  collapsed <- ann[!is.na(feature_category), {
    categories <- unique(feature_category)
    best_category <- priority[min(match(categories, priority))]

    gene_names <- if ("gene_name" %in% names(.SD)) {
      unique(gene_name[feature_category == best_category])
    } else {
      NA_character_
    }

    gene_names <- gene_names[!is.na(gene_names) & gene_names != ""]

    .(
      strand = strand[1],
      feature_category = best_category,
      gene_name = ifelse(
        length(gene_names),
        paste(gene_names, collapse = ","),
        NA_character_
      )
    )
  }, by = .(dmr_id, chr, start, end)]

  # Add intergenic DMRs (those not in annotation)
  missing_dmrs <- setdiff(dmrs$dmr_id, collapsed$dmr_id)

  if (length(missing_dmrs) > 0) {
    intergenic <- dmrs[dmr_id %in% missing_dmrs, .(dmr_id, chr, start, end)]
    intergenic[, `:=`(
      strand = "*",
      feature_category = "Intergenic",
      gene_name = NA_character_
    )]

    collapsed <- rbind(collapsed, intergenic, use.names = TRUE, fill = TRUE)
  }

  # Merge with methylation data
  collapsed <- merge(
    collapsed,
    dmrs[, .(dmr_id, meth_diff, qvalue)],
    by = "dmr_id",
    all.x = TRUE
  )

  collapsed[, `:=`(comparison = comparison, class = class)]

  # Write collapsed output
  output_file <- file.path(
    outdir,
    paste0(comparison, "_Merged_DMRs_cov4_", class, "_singleFeature_q05_meth5.tsv")
  )

  fwrite(collapsed[, !"dmr_id"], output_file, sep = "\t")

  # Create summary
  summary <- collapsed[, .N, by = feature_category][, `:=`(
    comparison = comparison,
    class = class
  )]

  summary_list[[paste(comparison, class)]] <- summary

  fwrite(
    summary,
    file.path(outdir, paste0("Summary_SingleFeature_", comparison, "_", class, "_q05_meth5.tsv")),
    sep = "\t"
  )

  message("Processed: ", comparison, " (", class, ")")
  message("  Total DMRs: ", nrow(collapsed))
}

# Write combined summary
if (length(summary_list) > 0) {
  fwrite(
    rbindlist(summary_list),
    file.path(outdir, "SingleFeature_counts_all_q05_meth5.tsv"),
    sep = "\t"
  )
}

message("[complete] All DMR feature collapsing finished.")
