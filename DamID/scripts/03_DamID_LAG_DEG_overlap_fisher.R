#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_DamID_LAG_DEG_overlap_fisher.R
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

# 03_DamID_LAG_DEG_overlap_fisher.R
#
# Define LAG gene sets from annotated DiffBind peaks, define DEG sets from a
# DESeq2 table, and compute LAG/DEG Fisher exact tests.
#
# Usage:
#   Rscript 03_DamID_LAG_DEG_overlap_fisher.R \
#     --lag-table DamID_24hr_KCl_vs_mannitol_DiffBind_peaks.csv \
#     --deg-table Hour_24_KClvsHour_24_Mann_deg.tsv \
#     --outdir DamID_DEG_overlap_24hr \
#     --lag-fdr 0.05 \
#     --deg-padj 0.01 \
#     --deg-fc 1.3
#
# Optional:
#   --universe expressed_gene_universe.txt

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

lag_table <- get_arg("--lag-table")
deg_table <- get_arg("--deg-table")
outdir <- get_arg("--outdir", "DamID_DEG_overlap_output")
universe_file <- get_arg("--universe", NULL)

lag_gene_col <- get_arg("--lag-gene-col", "gene_name")
lag_direction_col <- get_arg("--lag-direction-col", "direction")
lag_feature_col <- get_arg("--lag-feature-col", "feature")
lag_fdr_col <- get_arg("--lag-fdr-col", "FDR")
lag_fdr_cut <- as.numeric(get_arg("--lag-fdr", "0.05"))
promoter_only <- get_arg("--promoter-only", "TRUE")

deg_gene_col <- get_arg("--deg-gene-col", "gene_name")
deg_log2fc_col <- get_arg("--deg-log2fc-col", "log2FoldChange")
deg_padj_col <- get_arg("--deg-padj-col", "padj")
deg_padj_cut <- as.numeric(get_arg("--deg-padj", "0.01"))
deg_fc_cut <- as.numeric(get_arg("--deg-fc", "1.3"))
deg_log2fc_cut <- log2(deg_fc_cut)

if (is.null(lag_table) || is.null(deg_table)) {
  stop("Usage: Rscript 03_DamID_LAG_DEG_overlap_fisher.R --lag-table <peaks.csv> --deg-table <DESeq2.tsv/csv>")
}
if (!file.exists(lag_table)) stop("lag-table not found: ", lag_table)
if (!file.exists(deg_table)) stop("deg-table not found: ", deg_table)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

writeLines(
  c(
    paste0("lag_fdr_cut=", lag_fdr_cut),
    paste0("deg_padj_cut=", deg_padj_cut),
    paste0("deg_fc_cut=", deg_fc_cut),
    paste0("deg_log2fc_cut=log2(", deg_fc_cut, ")=", deg_log2fc_cut),
    paste0("promoter_only=", promoter_only)
  ),
  file.path(outdir, "analysis_parameters.txt")
)

read_any <- function(path) {
  if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    read_csv(path, show_col_types = FALSE)
  } else {
    read_tsv(path, show_col_types = FALSE)
  }
}

lag <- read_any(lag_table)
deg <- read_any(deg_table)

required_lag <- c(lag_gene_col, lag_direction_col, lag_fdr_col)
required_deg <- c(deg_gene_col, deg_log2fc_col, deg_padj_col)
missing_lag <- setdiff(required_lag, names(lag))
missing_deg <- setdiff(required_deg, names(deg))
if (length(missing_lag) > 0) stop("Missing columns in lag-table: ", paste(missing_lag, collapse = ", "))
if (length(missing_deg) > 0) stop("Missing columns in deg-table: ", paste(missing_deg, collapse = ", "))

lag2 <- lag %>%
  transmute(
    gene = as.character(.data[[lag_gene_col]]),
    direction = as.character(.data[[lag_direction_col]]),
    feature = if (lag_feature_col %in% names(lag)) as.character(.data[[lag_feature_col]]) else NA_character_,
    FDR = suppressWarnings(as.numeric(.data[[lag_fdr_col]]))
  ) %>%
  filter(!is.na(gene), gene != "", gene != "NA", !is.na(FDR), FDR < lag_fdr_cut)

if (toupper(promoter_only) %in% c("TRUE", "T", "1", "YES", "Y")) {
  lag2 <- lag2 %>% filter(.data$feature == "promoter")
}

deg2 <- deg %>%
  transmute(
    gene = as.character(.data[[deg_gene_col]]),
    log2FoldChange = suppressWarnings(as.numeric(.data[[deg_log2fc_col]])),
    padj = suppressWarnings(as.numeric(.data[[deg_padj_col]]))
  ) %>%
  filter(!is.na(gene), gene != "", gene != "NA")

if (!is.null(universe_file)) {
  universe <- unique(readLines(universe_file))
  universe <- universe[!is.na(universe) & universe != ""]
} else {
  universe <- unique(deg2$gene)
}

b1_loss_genes <- lag2 %>%
  filter(.data$direction %in% c("KCl_loss", "B1_loss", "loss")) %>%
  pull(gene) %>%
  unique() %>%
  intersect(universe)

b1_gain_genes <- lag2 %>%
  filter(.data$direction %in% c("KCl_gain", "B1_gain", "gain")) %>%
  pull(gene) %>%
  unique() %>%
  intersect(universe)

deg_up_genes <- deg2 %>%
  filter(!is.na(.data$padj), .data$padj < deg_padj_cut,
         !is.na(.data$log2FoldChange), .data$log2FoldChange >= deg_log2fc_cut) %>%
  pull(gene) %>%
  unique() %>%
  intersect(universe)

deg_down_genes <- deg2 %>%
  filter(!is.na(.data$padj), .data$padj < deg_padj_cut,
         !is.na(.data$log2FoldChange), .data$log2FoldChange <= -deg_log2fc_cut) %>%
  pull(gene) %>%
  unique() %>%
  intersect(universe)

run_fisher <- function(set_a, set_b, label_a, label_b) {
  a <- length(intersect(set_a, set_b))
  b <- length(setdiff(set_a, set_b))
  c <- length(setdiff(set_b, set_a))
  d <- length(setdiff(universe, union(set_a, set_b)))

  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))

  tibble(
    comparison = paste(label_a, label_b, sep = "_vs_"),
    overlap = a,
    set_a_only = b,
    set_b_only = c,
    neither = d,
    n_set_a = length(set_a),
    n_set_b = length(set_b),
    n_universe = length(universe),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value,
    conf_low = ft$conf.int[1],
    conf_high = ft$conf.int[2]
  )
}

results <- bind_rows(
  run_fisher(b1_loss_genes, deg_up_genes, "b1_loss", "deg_up"),
  run_fisher(b1_loss_genes, deg_down_genes, "b1_loss", "deg_down"),
  run_fisher(b1_gain_genes, deg_up_genes, "b1_gain", "deg_up"),
  run_fisher(b1_gain_genes, deg_down_genes, "b1_gain", "deg_down")
) %>%
  mutate(FDR_BH = p.adjust(.data$p_value, method = "BH"))

write_csv(results, file.path(outdir, "DamID_LAG_DEG_overlap_fisher_results.csv"))
writeLines(sort(b1_loss_genes), file.path(outdir, "b1_loss_genes.txt"))
writeLines(sort(b1_gain_genes), file.path(outdir, "b1_gain_genes.txt"))
writeLines(sort(deg_up_genes), file.path(outdir, "deg_up_genes.txt"))
writeLines(sort(deg_down_genes), file.path(outdir, "deg_down_genes.txt"))
writeLines(sort(universe), file.path(outdir, "gene_universe.txt"))

message("[done] wrote results to: ", normalizePath(outdir))
print(results)
