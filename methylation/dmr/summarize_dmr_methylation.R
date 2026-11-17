#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


base_dirs <- c(
  gb_path("wgbs/methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp")
)


sf_subdir <- file.path("Merged_q05_meth5", "SingleFeature_q05_meth5")


map_feature <- function(x) {
  x <- tolower(x)
  out <- ifelse(grepl("intergenic", x), "Intergenic",
         ifelse(grepl("promoter|tss", x), "Promoter/TSS",
         ifelse(grepl("exon", x), "Exon",
         ifelse(grepl("intron", x), "Intron",

         ifelse(grepl("\\b5.?utr|utr5", x), "Promoter/TSS",
         ifelse(grepl("\\b3.?utr|utr3", x), "Exon", NA_character_))))))
  factor(out, levels = c("Promoter/TSS","Exon","Intron","Intergenic"))
}

class_factor <- function(x) factor(x, levels = c("hypo","hyper"))


summarize_one_contrast <- function(base_dir) {
  sf_dir <- file.path(base_dir, sf_subdir)
  if (!dir.exists(sf_dir)) {
    cat("Single-feature dir missing:", sf_dir, "\n")
    return(invisible(NULL))
  }

  files <- list.files(
  sf_dir,
  pattern = "_Merged_DMRs_cov4_(hyper|hypo)_singleFeature_q05_meth5\\.tsv$",
  full.names = TRUE
)

  if (!length(files)) {
    cat(" No single-feature TSVs in:", sf_dir, "\n")
    return(invisible(NULL))
  }

  comparison <- sub("^DMR_", "", sub("_multiwindow$", "", basename(dirname(base_dir))))
  out_dir <- sf_dir

  cat("\nContrast:", comparison, "\n")
  cat("   Input dir:", sf_dir, "\n")
  cat("   Files:", length(files), "\n")

  dt_list <- lapply(files, function(f) {
    x <- fread(f, sep = "\t", header = TRUE, data.table = TRUE, showProgress = FALSE)


    if ("meth.diff" %in% names(x)) setnames(x, "meth.diff", "meth_diff")


    feat_col <- if ("feature_category" %in% names(x)) "feature_category" else
                if ("feature" %in% names(x)) "feature" else NA_character_

    if (is.na(feat_col)) {
      stop("No feature or feature_category column in ", f)
    }


    req <- c("chr","start","end","meth_diff","qvalue", feat_col)
    missing <- setdiff(req, names(x))
    if (length(missing)) {
      stop("Missing required columns in ", f, " : ", paste(missing, collapse = ", "))
    }


    if ("DMR_type" %in% names(x)) {
      class_vec <- tolower(x$DMR_type)
    } else if ("class" %in% names(x)) {
      class_vec <- tolower(x$class)
    } else if (grepl("hyper", basename(f), ignore.case = TRUE)) {
      class_vec <- rep("hyper", nrow(x))
    } else if (grepl("hypo", basename(f), ignore.case = TRUE)) {
      class_vec <- rep("hypo", nrow(x))
    } else {
      stop("No DMR_type/class info in ", f, " and cannot infer from file name.")
    }

    x[, `:=`(
      contrast         = comparison,
      class            = class_vec,
      dmr_id           = paste(chr, start, end, sep = "_"),
      abs_meth_diff    = abs(meth_diff),
      feature_category = map_feature(get(feat_col))
    )][]
  })

  long_dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)


  long_dt <- long_dt[!is.na(feature_category)]
  long_dt[, class := class_factor(class)]


  summary_dt <- long_dt[
    , .(
        N             = .N,
        median_abs    = as.numeric(median(abs_meth_diff, na.rm = TRUE)),
        IQR25_abs     = as.numeric(quantile(abs_meth_diff, 0.25, na.rm = TRUE)),
        IQR75_abs     = as.numeric(quantile(abs_meth_diff, 0.75, na.rm = TRUE)),
        median_signed = as.numeric(median(meth_diff, na.rm = TRUE))
      ),
    by = .(comparison, feature_category, class)
  ][order(feature_category, class)]


  long_out    <- file.path(out_dir, paste0(comparison, "_feature_methylation_long.tsv"))
  summary_out <- file.path(out_dir, paste0(comparison, "_feature_methylation_summary.tsv"))

  fwrite(long_dt[, .(comparison, class, feature_category,
                     chr, start, end, qvalue, meth_diff, abs_meth_diff, dmr_id)],
         long_out, sep = "\t", quote = FALSE, na = "NA")

  fwrite(summary_dt, summary_out, sep = "\t", quote = FALSE, na = "NA")


  xlsx_out <- file.path(out_dir, paste0(comparison, "_feature_methylation.xlsx"))
  wb <- createWorkbook()
  addWorksheet(wb, "long")
  writeData(wb, "long",
            as.data.frame(long_dt[, .(
              comparison, class, feature_category,
              chr, start, end, qvalue, meth_diff, abs_meth_diff, dmr_id
            )]))
  addWorksheet(wb, "summary")
  writeData(wb, "summary", as.data.frame(summary_dt))
  saveWorkbook(wb, xlsx_out, overwrite = TRUE)

  cat("", comparison, " →\n  ", long_out, "\n  ",
      summary_out, "\n  ", xlsx_out, "\n", sep = "")
  invisible(TRUE)
}


for (bd in base_dirs) summarize_one_contrast(bd)

cat("DONE\n")
