#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_DamID_DiffBind_GTF_annotation.R
# Package: DamID / LAG analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for DamID preprocessing, DiffBind/GTF annotation, LAG/DEG overlap analysis, and GO enrichment.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# 01_DamID_DiffBind_GTF_annotation.R
#

# Run DiffBind analysis for DamID peak inputs and annotate differential peaks
# using custom mm10 GTF-derived feature sets.
#
# Usage:
#   Rscript 01_DamID_DiffBind_GTF_annotation.R \
#     --workdir /path/to/diffbind_folder \
#     --sample-sheet Table1_24.csv \
#     --gtf-features gtf_dmrs_features_plus_minus1kb.rds \
#     --output DiffBind_B1_annotated_full.csv \
#     --reference-condition Mannitol \
#     --report-threshold 1
#
# The GTF features RDS should contain GRanges objects named:
#   promoters_plus_minus1kb, exons, introns, fiveUTRs, threeUTRs, genes

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(DiffBind)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

workdir <- get_arg("--workdir", ".")
sample_sheet <- get_arg("--sample-sheet", "Table1.csv")
gtf_rds <- get_arg("--gtf-features", "gtf_dmrs_features_plus_minus1kb.rds")
output_csv <- get_arg("--output", "DiffBind_B1_annotated_full.csv")
reference_condition <- get_arg("--reference-condition", "Mannitol")
report_threshold <- as.numeric(get_arg("--report-threshold", "1"))

if (!dir.exists(workdir)) stop("workdir does not exist: ", workdir)
setwd(workdir)

if (!file.exists(sample_sheet)) stop("sample sheet not found: ", sample_sheet)
if (!file.exists(gtf_rds)) stop("GTF feature RDS not found: ", gtf_rds)

gtf_features <- readRDS(gtf_rds)

required_features <- c(
  "promoters_plus_minus1kb",
  "exons",
  "introns",
  "fiveUTRs",
  "threeUTRs",
  "genes"
)

missing_features <- setdiff(required_features, names(gtf_features))
if (length(missing_features) > 0) {
  stop("Missing required GTF feature sets: ", paste(missing_features, collapse = ", "))
}

features <- list(
  promoter = trim(gtf_features[["promoters_plus_minus1kb"]]),
  exon     = trim(gtf_features[["exons"]]),
  intron   = trim(gtf_features[["introns"]]),
  fiveUTR  = trim(gtf_features[["fiveUTRs"]]),
  threeUTR = trim(gtf_features[["threeUTRs"]]),
  gene     = trim(gtf_features[["genes"]])
)

message("[DiffBind] Loading sample sheet: ", sample_sheet)
db <- dba(sampleSheet = sample_sheet)

db <- dba.count(
  db,
  bSubControl = TRUE,
  minOverlap = 1,
  fragmentSize = 0,
  summits = TRUE
)

db <- dba.normalize(db)
db <- dba.contrast(db, reorderMeta = list(Condition = reference_condition))
db <- dba.analyze(db)

db_all <- trim(dba.report(db, th = report_threshold))

seqlevelsStyle(db_all) <- "UCSC"
for (nm in names(features)) {
  seqlevelsStyle(features[[nm]]) <- "UCSC"
}

peak_gr <- GRanges(seqnames = seqnames(db_all), ranges = ranges(db_all))

annotate_with_feature <- function(feature_gr, feature_name) {
  hits <- findOverlaps(peak_gr, feature_gr, minoverlap = 1)
  if (length(hits) == 0) return(NULL)

  gene_col <- if ("gene_id" %in% colnames(mcols(feature_gr))) {
    mcols(feature_gr)$gene_id[subjectHits(hits)]
  } else if ("gene_name" %in% colnames(mcols(feature_gr))) {
    mcols(feature_gr)$gene_name[subjectHits(hits)]
  } else if ("transcript_id" %in% colnames(mcols(feature_gr))) {
    mcols(feature_gr)$transcript_id[subjectHits(hits)]
  } else {
    rep(NA_character_, length(subjectHits(hits)))
  }

  data.frame(
    peak_id = queryHits(hits),
    feature = feature_name,
    gene_id = gene_col,
    stringsAsFactors = FALSE
  )
}

# Priority is determined by the feature order below. The first annotation per peak is kept.
feature_order <- c("promoter", "exon", "intron", "fiveUTR", "threeUTR", "gene")

all_hits <- bind_rows(lapply(feature_order, function(nm) {
  annotate_with_feature(features[[nm]], nm)
})) %>%
  distinct(peak_id, .keep_all = TRUE)

db_df <- as.data.frame(db_all)
db_df$peak_id <- seq_len(nrow(db_df))

final_annot <- left_join(db_df, all_hits, by = "peak_id")
write.csv(final_annot, output_csv, row.names = FALSE)

message("[done] wrote annotated DiffBind peak table: ", output_csv)
message("[done] rows: ", nrow(final_annot))
