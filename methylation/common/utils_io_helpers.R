#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

read_deg_table <- function(path) {
  df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  df$gene_symbol <- df$gene_name

  fpkm_cols <- grep("_fpkm$", names(df), value = TRUE)
  if (length(fpkm_cols) == 0) {
    warning("No _fpkm columns found in ", path, " -- check real column names.")
    df$max_fpkm <- NA_real_
  } else {
    df$max_fpkm <- apply(df[, fpkm_cols, drop = FALSE], 1, max, na.rm = TRUE)
  }

  df
}

filter_deg_significant <- function(deg_df, deg_qval = DEG_QVAL, deg_fc = DEG_FC) {
  deg_df[deg_df$padj <= deg_qval & abs(deg_df$log2FoldChange) >= log2(deg_fc), ]
}

expressed_genes_fpkm <- function(deg_df, cutoff) {
  unique(deg_df$gene_id[deg_df$max_fpkm > cutoff])
}

read_promoter_dmg_table <- function(path_list, dmg_only = TRUE,
                                     dmg_delta = DMG_DELTA, dmg_qval = DMG_QVAL) {
  hyper <- utils::read.delim(path_list$hyper, stringsAsFactors = FALSE)
  hypo  <- utils::read.delim(path_list$hypo, stringsAsFactors = FALSE)
  hyper$direction <- "hyper"
  hypo$direction  <- "hypo"
  df <- rbind(hyper, hypo)
  df <- dplyr::rename(df, meth.diff = meth_diff, gene_symbol = gene_name)
  df$gene_id <- df$gene_symbol

  if (dmg_only) {
    df <- df[abs(df$meth.diff) >= dmg_delta & df$qvalue <= dmg_qval, ]
  }
  df
}
