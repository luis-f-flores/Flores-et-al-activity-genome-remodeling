#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_IAP_location_summary.R
# Package: Transposable element analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for TE catalog construction, TE differential-expression annotation, TE-gene overlap, IAP summaries, methylation metaprofiles, and heatmaps.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# Summarize genomic locations of significantly upregulated IAP-family TEs.
# Usage:
#   Rscript 04_IAP_location_summary.R --de TE_DE_results_annotated.tsv --location-summary TE_family_location_summary.tsv --out IAP_up_location_3part.xlsx --padj 0.05 --log2fc-min 0

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(writexl)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value for ", flag)
  args[[idx + 1]]
}

de_tsv <- get_arg("--de")
loc_tsv <- get_arg("--location-summary")
out_xlsx <- get_arg("--out", "IAP_up_location_3part.xlsx")
padj_cut <- as.numeric(get_arg("--padj", "0.05"))
log2fc_min <- as.numeric(get_arg("--log2fc-min", "0"))

if (is.null(de_tsv) || is.null(loc_tsv)) {
  stop("Usage: Rscript 04_IAP_location_summary.R --de <DE.tsv> --location-summary <summary.tsv> --out <out.xlsx>")
}

de <- read_tsv(de_tsv, show_col_types = FALSE) %>% filter(!is.na(padj), !is.na(log2FoldChange))
loc <- read_tsv(loc_tsv, show_col_types = FALSE) %>%
  mutate(
    n_total      = as.numeric(n_total),
    n_promoter   = as.numeric(n_promoter),
    n_exon       = as.numeric(n_exon),
    n_intron     = as.numeric(n_intron),
    n_intergenic = as.numeric(n_intergenic)
  )

iap_up <- de %>%
  filter(grepl("^IAP", family), padj < padj_cut, log2FoldChange > log2fc_min) %>%
  select(family, log2FoldChange, padj) %>%
  arrange(desc(log2FoldChange))

if (nrow(iap_up) == 0) stop("No upregulated IAP families found using the requested thresholds.")

dat <- iap_up %>% inner_join(loc, by = "family")
if (nrow(dat) == 0) stop("Upregulated IAP families were not found in the location summary table.")

sum_promoter <- sum(dat$n_promoter, na.rm = TRUE)
sum_exon <- sum(dat$n_exon, na.rm = TRUE)
sum_intron <- sum(dat$n_intron, na.rm = TRUE)
sum_intergenic <- sum(dat$n_intergenic, na.rm = TRUE)
sum_total <- sum(dat$n_total, na.rm = TRUE)

three_part <- tibble(
  Category = c("Promoter", "Intron+Exon", "Intergenic"),
  Count = c(sum_promoter, sum_exon + sum_intron, sum_intergenic)
) %>% mutate(Pct_of_total = round(Count / sum_total, 4))

families_used <- dat %>%
  select(family, log2FoldChange, padj, n_total, n_promoter, n_exon, n_intron, n_intergenic) %>%
  arrange(desc(log2FoldChange))

readme <- tibble(
  Item = c("DE table", "Location summary", "Family filter", "Significance", "Direction", "Output categories"),
  Value = c(de_tsv, loc_tsv, "^IAP.*", paste0("padj < ", padj_cut), paste0("log2FC > ", log2fc_min), "Promoter | Intron+Exon | Intergenic")
)

write_xlsx(list("00_README" = readme, "01_3part_counts" = three_part, "02_IAP_families" = families_used), out_xlsx)
message("Wrote: ", out_xlsx)
print(three_part)
