#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


source("config.R")

suppressPackageStartupMessages({
  library(dplyr)
})

all_diff <- do.call(rbind, lapply(TIMEPOINTS, function(tp) {
  path <- file.path(OUTPUT_TABLES_DIR, sprintf("08_genebody_methylkit_direct_%s_all.csv", tp))
  if (!file.exists(path)) stop("Run gene_body_methylation.R for ", tp, " first -- ", path, " not found.")
  df <- read.csv(path, stringsAsFactors = FALSE)
  df$timepoint <- tp
  df
}))

all_diff$significant <- abs(all_diff$meth.diff) >= DMR_DELTA & all_diff$qvalue <= DMR_QVAL
all_diff$direction <- ifelse(all_diff$meth.diff > 0, "Hyper", "Hypo")
all_diff$timepoint <- factor(all_diff$timepoint, levels = TIMEPOINTS)

summary_tbl <- all_diff %>%
  dplyr::group_by(timepoint) %>%
  dplyr::summarise(
    n_genes_tested = dplyr::n(),
    n_significant = sum(significant),
    n_hyper = sum(significant & direction == "Hyper"),
    n_hypo = sum(significant & direction == "Hypo"),
    median_delta_significant = stats::median(meth.diff[significant]),
    min_delta_significant = min(meth.diff[significant]),
    max_delta_significant = max(meth.diff[significant]),
    .groups = "drop"
  )

write.csv(summary_tbl,
          file.path(OUTPUT_TABLES_DIR, "12_genebody_significance_summary.csv"),
          row.names = FALSE)

cat("Gene-body significance summary (q<=", DMR_QVAL, "/|delta|>=", DMR_DELTA, "%):\n")
print(summary_tbl)

cat("Gene-body significance summary written to", OUTPUT_TABLES_DIR, "\n")
