#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_GTF_annotate_DMRs.R
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

# 04_GTF_annotate_DMRs.R
#
# Annotate DMRs with genomic features from GTF-derived GRanges objects.
#
# Usage:
#   Rscript 04_GTF_annotate_DMRs.R \
#     --input-dir DMR_merged \
#     --gtf-features gtf_features_plus_minus1kb.rds \
#     --pattern "*_Merged_DMRs_cov4_*.txt"

suppressPackageStartupMessages({
  library(GenomicRanges)
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

input_dir <- get_arg("--input-dir")
gtf_rds <- get_arg("--gtf-features")
pattern <- get_arg("--pattern", "*_Merged_DMRs_cov4_*.txt")

if (is.null(input_dir) || is.null(gtf_rds)) {
  stop("Required arguments: --input-dir and --gtf-features")
}

# Load GTF feature sets
gtf <- readRDS(gtf_rds)

# Find all DMR files
files <- list.files(
  input_dir,
  pattern = pattern,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files found matching pattern: ", file.path(input_dir, pattern))
}

for (f in files) {
  dmrs <- fread(f)

  # Convert DMRs to GRanges
  dmr_gr <- GRanges(
    seqnames = dmrs$chr,
    ranges = IRanges(start = dmrs$start, end = dmrs$end)
  )

  # Annotate with each feature type
  hits_all <- list()

  for (feature_name in names(gtf)) {
    feat <- gtf[[feature_name]]

    # Handle GRangesLists
    if (inherits(feat, "GRangesList")) {
      feat <- unlist(feat, use.names = FALSE)
    }

    # Find overlaps
    overlaps <- findOverlaps(dmr_gr, feat, minoverlap = 1, ignore.strand = TRUE)

    if (length(overlaps) > 0) {
      metadata_cols <- mcols(feat)[subjectHits(overlaps), , drop = FALSE]

      hits_all[[feature_name]] <- data.frame(
        chr = as.character(seqnames(dmr_gr[queryHits(overlaps)])),
        start = start(dmr_gr[queryHits(overlaps)]),
        end = end(dmr_gr[queryHits(overlaps)]),
        strand = as.character(strand(feat)[subjectHits(overlaps)]),
        feature = feature_name,
        tx_id = if ("tx_name" %in% colnames(metadata_cols)) {
          metadata_cols$tx_name
        } else {
          NA_character_
        },
        gene_name = if ("gene_id" %in% colnames(metadata_cols)) {
          as.character(metadata_cols$gene_id)
        } else {
          NA_character_
        },
        stringsAsFactors = FALSE
      )
    }
  }

  # Combine all annotations
  annotations <- bind_rows(hits_all)

  # Add DMRs with no feature overlap as "intergenic"
  dmrs$key <- paste(dmrs$chr, dmrs$start, dmrs$end, sep = "_")

  annotated_keys <- if (nrow(annotations) > 0) {
    unique(paste(annotations$chr, annotations$start, annotations$end, sep = "_"))
  } else {
    character()
  }

  intergenic_dmrs <- dmrs[!dmrs$key %in% annotated_keys, .(chr, start, end)]

  if (nrow(intergenic_dmrs) > 0) {
    intergenic_dmrs[, `:=`(
      strand = "*",
      feature = "intergenic",
      tx_id = NA_character_,
      gene_name = NA_character_
    )]

    annotations <- bind_rows(annotations, as.data.frame(intergenic_dmrs))
  }

  # Write annotated output
  output_file <- file.path(
    dirname(f),
    paste0("GTF_Annot_", basename(f))
  )

  fwrite(as.data.table(annotations), output_file, sep = "\t")

  message("Annotated: ", basename(f))
  message("  Total regions: ", nrow(annotations))
  message("  Output: ", output_file)
}

message("[complete] All DMR annotations finished.")
