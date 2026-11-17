#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
})


q_cutoff    <- 0.01
meth_cutoff <- 25

base_dir <- getwd()

promoter_files <- c(
  "1M_vs_1K_hyper_Promoter_DMRs_ForVolcano.txt",
  "1M_vs_1K_hypo_Promoter_DMRs_ForVolcano.txt",
  "6M_vs_6K_hyper_Promoter_DMRs_ForVolcano.txt",
  "6M_vs_6K_hypo_Promoter_DMRs_ForVolcano.txt",
  "24M_vs_24K_hyper_Promoter_DMRs_ForVolcano.txt",
  "24M_vs_24K_hypo_Promoter_DMRs_ForVolcano.txt"
)

for (fn in promoter_files) {
  in_path  <- file.path(base_dir, fn)
  if (!file.exists(in_path)) {
    warning("Skipping missing file: ", in_path)
    next
  }

  message("Processing: ", fn)
  df <- read.delim(in_path, check.names = FALSE)


  required_cols <- c("gene_name", "meth_diff", "qvalue")
  miss <- setdiff(required_cols, colnames(df))
  if (length(miss) > 0) {
    warning("  -> Missing columns in ", fn, ": ",
            paste(miss, collapse = ", "), "  (skipping)")
    next
  }

  n_in <- nrow(df)

  df_filt <- df %>%
    filter(!is.na(gene_name), gene_name != "") %>%
    filter(!is.na(qvalue), qvalue <= q_cutoff) %>%
    filter(!is.na(meth_diff), abs(meth_diff) >= meth_cutoff)

  n_out <- nrow(df_filt)

  out_path <- file.path(
    base_dir,
    sub("\\.txt$", "_meth25_q01.txt", fn)
  )

  write.table(
    df_filt,
    file      = out_path,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )

  message("  kept ", n_out, " / ", n_in,
          " rows  ->  wrote: ", basename(out_path))
}

message("Done filtering promoter DMR files.")
