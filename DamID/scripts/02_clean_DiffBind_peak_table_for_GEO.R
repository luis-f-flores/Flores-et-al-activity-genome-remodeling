#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_clean_DiffBind_peak_table_for_GEO.R
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

# 02_clean_DiffBind_peak_table_for_GEO.R
#
# Clean an annotated DiffBind peak table into the compact processed peak CSV
# used for GEO deposition.
#
# Usage:
#   Rscript 02_clean_DiffBind_peak_table_for_GEO.R \
#     --input DiffBind_B1_annotated_full.csv \
#     --output DamID_24hr_KCl_vs_mannitol_DiffBind_peaks.csv

suppressPackageStartupMessages({
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

input_file <- get_arg("--input")
output_file <- get_arg("--output")

if (is.null(input_file) || is.null(output_file)) {
  stop("Usage: Rscript 02_clean_DiffBind_peak_table_for_GEO.R --input <annotated.csv> --output <clean.csv>")
}
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

df <- read_csv(input_file, show_col_types = FALSE)

if (!"gene_id.value" %in% names(df)) df$gene_id.value <- NA_character_
if (!"gene_id" %in% names(df)) df$gene_id <- NA_character_
if (!"feature" %in% names(df)) df$feature <- NA_character_
if (!"peak_id" %in% names(df)) df$peak_id <- seq_len(nrow(df))

required <- c("seqnames", "start", "end", "width", "strand", "Conc", "Fold", "p.value", "FDR", "peak_id")
missing_required <- setdiff(required, names(df))
if (length(missing_required) > 0) {
  stop("Missing required columns: ", paste(missing_required, collapse = ", "))
}

if (!"Conc_KCl" %in% names(df)) df$Conc_KCl <- NA_real_
if (!"Conc_Mannitol" %in% names(df)) df$Conc_Mannitol <- NA_real_

df_clean <- df %>%
  mutate(
    gene_id.value = as.character(.data$gene_id.value),
    gene_id = as.character(.data$gene_id),
    gene_name = case_when(
      !is.na(.data$gene_id.value) & .data$gene_id.value != "NA" & .data$gene_id.value != "" ~ .data$gene_id.value,
      !is.na(.data$gene_id) & .data$gene_id != "NA" & .data$gene_id != "" ~ .data$gene_id,
      TRUE ~ NA_character_
    ),
    feature = if_else(
      is.na(.data$feature) | .data$feature == "NA" | .data$feature == "",
      "intergenic_or_unannotated",
      as.character(.data$feature)
    ),
    direction = case_when(
      .data$Fold > 0 ~ "KCl_gain",
      .data$Fold < 0 ~ "KCl_loss",
      TRUE ~ "no_change"
    )
  ) %>%
  select(
    seqnames, start, end, width, strand,
    Conc, Conc_KCl, Conc_Mannitol,
    Fold, p.value, FDR,
    peak_id, direction, feature, gene_name
  )

write_csv(df_clean, output_file)

message("[done] wrote: ", output_file)
message("[info] rows: ", nrow(df_clean))
