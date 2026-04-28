#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_GO_enrichment_B1_bound_DEGs.R
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

# 04_GO_enrichment_B1_bound_DEGs.R
#
# Perform GO enrichment on B1-bound / LAG-associated DEG gene lists.
#
# The input can be either:
#   1) one gene per line, with gene symbols or Entrez IDs, or
#   2) a two-column table with gene symbol in column 1 and Entrez ID in column 2.
#
# Usage:
#   Rscript 04_GO_enrichment_B1_bound_DEGs.R \
#     --input-dir gene_lists \
#     --pattern "*.txt" \
#     --outdir GO_Plots \
#     --id-type SYMBOL
#
# Optional:
#   --universe expressed_gene_universe.txt
#   --ont BP
#   --p-cutoff 0.05
#   --q-cutoff 0.05
#   --show-category 12

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(enrichplot)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

input_dir <- get_arg("--input-dir", ".")
pattern <- get_arg("--pattern", "*.txt")
outdir <- get_arg("--outdir", "GO_Plots")
id_type <- toupper(get_arg("--id-type", "SYMBOL"))  # SYMBOL or ENTREZID
universe_file <- get_arg("--universe", NULL)
ont <- toupper(get_arg("--ont", "BP"))
p_cutoff <- as.numeric(get_arg("--p-cutoff", "0.05"))
q_cutoff <- as.numeric(get_arg("--q-cutoff", "0.05"))
show_category <- as.integer(get_arg("--show-category", "12"))

if (!dir.exists(input_dir)) stop("input-dir not found: ", input_dir)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

files <- Sys.glob(file.path(input_dir, pattern))
if (length(files) == 0) stop("No input files matched: ", file.path(input_dir, pattern))

read_gene_file <- function(path, id_type) {
  x <- read.table(path, header = FALSE, stringsAsFactors = FALSE, sep = "", quote = "", comment.char = "")
  if (ncol(x) >= 2 && id_type == "ENTREZID") {
    genes <- as.character(x[[2]])
  } else {
    genes <- as.character(x[[1]])
  }
  unique(genes[!is.na(genes) & genes != "" & genes != "NA"])
}

convert_to_entrez <- function(genes, id_type) {
  if (id_type == "ENTREZID") return(unique(as.character(genes)))

  mapped <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = unique(genes),
    columns = c("ENTREZID"),
    keytype = id_type
  )

  unique(mapped$ENTREZID[!is.na(mapped$ENTREZID)])
}

universe_entrez <- NULL
if (!is.null(universe_file)) {
  universe_genes <- readLines(universe_file)
  universe_genes <- unique(universe_genes[!is.na(universe_genes) & universe_genes != ""])
  universe_entrez <- convert_to_entrez(universe_genes, id_type)
}

simplify_if_possible <- function(go_obj) {
  if (!is.null(go_obj) && nrow(as.data.frame(go_obj)) > 0) {
    simplify(go_obj, cutoff = 0.7, by = "p.adjust", select_fun = min)
  } else {
    go_obj
  }
}

save_dotplot <- function(enrichment_result, filename_base, title) {
  df <- as.data.frame(enrichment_result)
  if (nrow(df) == 0) {
    message("[no terms] ", title)
    return(invisible(NULL))
  }

  plot_obj <- dotplot(enrichment_result, showCategory = show_category, title = title) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8))

  pdf(file = file.path(outdir, paste0(filename_base, ".pdf")), width = 7, height = 5)
  print(plot_obj)
  dev.off()

  png(file = file.path(outdir, paste0(filename_base, ".png")), width = 2000, height = 1600, res = 300)
  print(plot_obj)
  dev.off()
}

for (file in files) {
  genes_raw <- read_gene_file(file, id_type)
  genes_entrez <- convert_to_entrez(genes_raw, id_type)

  base_filename <- tools::file_path_sans_ext(basename(file))

  if (length(genes_entrez) == 0) {
    message("[skip] no valid Entrez IDs for: ", basename(file))
    next
  }

  ego <- enrichGO(
    gene = genes_entrez,
    universe = universe_entrez,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    readable = TRUE
  )

  ego <- simplify_if_possible(ego)

  out_table <- file.path(outdir, paste0("GO_", ont, "_", base_filename, "_enrichment.tsv"))
  write.table(as.data.frame(ego), file = out_table, sep = "\t", quote = FALSE, row.names = FALSE)

  save_dotplot(
    ego,
    paste0("GO_", ont, "_", base_filename, "_dotplot"),
    paste0("GO ", ont, ": ", base_filename)
  )

  message("[done] ", basename(file), " -> ", out_table)
}

message("[done] GO enrichment complete. Output directory: ", normalizePath(outdir))
