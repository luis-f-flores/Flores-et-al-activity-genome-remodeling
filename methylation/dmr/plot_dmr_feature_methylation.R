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
  library(ggplot2)
})


base_dirs <- c(
  gb_path("wgbs/methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp")
)

sf_subdir   <- file.path("Merged_q05_meth5", "SingleFeature_q05_meth5")
feat_pat    <- "_Merged_DMRs_cov4_(hyper|hypo)_singleFeature_q05_meth5\\.tsv$"
categories  <- c("Promoter","Exon","Intron","Intergenic","CpG island")


get_script_dir <- function(){
  args <- commandArgs(trailingOnly = FALSE)
  farg <- "--file="; p <- sub(farg,"",args[grepl(farg,args)])
  if (length(p)) return(normalizePath(dirname(p)))
  if (!is.null(sys.frames()[[1]]$ofile)) return(normalizePath(dirname(sys.frames()[[1]]$ofile)))
  tryCatch(normalizePath(getwd()), error=function(e) getwd())
}

script_dir <- get_script_dir()
out_root <- file.path(script_dir, "q05_meth5_SingleFeature_CGI_medianIQR_outputs_signed")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_root, "\n")


contrast_name <- function(bdir) sub("^DMR_", "", sub("_multiwindow$", "", basename(dirname(bdir))))
tp_from_contrast <- function(cn){
  if (grepl("^1M_", cn))  return("1h")
  if (grepl("^6M_", cn))  return("6h")
  if (grepl("^24M_", cn)) return("24h")
  return(cn)
}


map_feature <- function(x) {
  x <- tolower(x)
  out <- ifelse(grepl("intergenic", x), "Intergenic",
         ifelse(grepl("promoter|tss", x), "Promoter",
         ifelse(grepl("exon", x), "Exon",
         ifelse(grepl("intron", x), "Intron", NA_character_))))
  factor(out, levels = c("Promoter","Exon","Intron","Intergenic"))
}

summarize_signed <- function(dt){

  dt[, .(
    N              = .N,
    median_signed  = as.numeric(median(meth_diff, na.rm = TRUE)),
    IQR25_signed   = as.numeric(quantile(meth_diff, 0.25, na.rm = TRUE)),
    IQR75_signed   = as.numeric(quantile(meth_diff, 0.75, na.rm = TRUE))
  ), by = .(contrast, feature_category, class, timepoint)]
}


read_single_feature_for_contrast <- function(sf_dir, comp, tp){
  files <- list.files(sf_dir, pattern = feat_pat, full.names = TRUE)
  if (!length(files)) {
    cat(" No SingleFeature q05 meth5 DMR files in:", sf_dir, "\n")
    return(NULL)
  }

  dt_list <- lapply(files, function(f) {
    x <- fread(f, sep = "\t")
    if ("meth.diff" %in% names(x)) setnames(x, "meth.diff", "meth_diff")

    needed <- c("chr", "start", "end", "meth_diff")
    if (!all(needed %in% names(x))) {
      stop("Missing required columns in ", f,
           " (need at least: ", paste(needed, collapse=", "), ")")
    }


    class <- if (grepl("_hyper_", basename(f))) "hyper" else "hypo"


    if ("feature_category" %in% names(x)) {
      x[, feature_category := as.character(feature_category)]
    } else if ("feature" %in% names(x)) {
      x[, feature_category := as.character(map_feature(feature))]
    } else {
      stop("No feature or feature_category column in ", f)
    }

    x[, `:=`(
      contrast      = comp,
      timepoint     = tp,
      class         = class
    )]


    x[, .(contrast,
          timepoint,
          class,
          feature_category,
          chr,
          start,
          end,
          meth_diff)]
  })

  dt <- rbindlist(dt_list, use.names = TRUE, fill = TRUE)

  dt[feature_category %in% c("Promoter","Exon","Intron","Intergenic")]
}


read_cgi_for_contrast <- function(sf_dir, comp, tp){
  cgi_long <- file.path(sf_dir, paste0(comp, "_cgi_methylation_long.tsv"))
  if (!file.exists(cgi_long)) {
    cat(" Missing CGI long file:", cgi_long, "\n")
    return(NULL)
  }
  C <- fread(cgi_long, sep = "\t")
  if ("meth.diff" %in% names(C)) setnames(C, "meth.diff", "meth_diff")

  needed <- c("class", "meth_diff")
  missing <- setdiff(needed, names(C))
  if (length(missing)) {
    stop("Missing columns in CGI long file ", cgi_long,
         " (need: ", paste(needed, collapse=", "), ")")
  }

  if (!"contrast" %in% names(C)) C[, contrast := comp]
  if (!"timepoint" %in% names(C)) C[, timepoint := tp]

  C[, feature_category := "CpG island"]
  C[, .(contrast,
        timepoint,
        class,
        feature_category,
        chr,
        start,
        end,
        meth_diff)]
}


all_sum <- list()

for (bd in base_dirs){
  comp <- contrast_name(bd)
  tp   <- tp_from_contrast(comp)

  sf_dir <- file.path(bd, sf_subdir)
  if (!dir.exists(sf_dir)) {
    cat("Missing SingleFeature dir:", sf_dir, "\n")
    next
  }

  cat("\nContrast:", comp, " (", tp, ")\n   SingleFeature dir:", sf_dir, "\n")

  F <- read_single_feature_for_contrast(sf_dir, comp, tp)
  if (is.null(F) || !nrow(F)) {
    cat(" No single-feature rows for:", comp, "\n")
    next
  }

  C <- read_cgi_for_contrast(sf_dir, comp, tp)
  if (is.null(C) || !nrow(C)) {
    cat(" No CGI rows for:", comp, "\n")
    next
  }


  F[, `:=`(
    feature_category = factor(feature_category,
                              levels = c("Promoter","Exon","Intron","Intergenic")),
    class = factor(class, levels = c("hypo","hyper")),
    timepoint = factor(tp, levels = c("1h","6h","24h"))
  )]

  C[, `:=`(
    feature_category = factor(feature_category, levels = categories),
    class = factor(class, levels = c("hypo","hyper")),
    timepoint = factor(tp, levels = c("1h","6h","24h"))
  )]

  S_feat <- summarize_signed(F)
  S_cgi  <- summarize_signed(C)

  all_sum[[comp]] <- rbind(S_feat, S_cgi, use.names = TRUE, fill = TRUE)
}

if (!length(all_sum)) stop("No summaries computed — ensure SingleFeature q05 meth5 and CGI long TSVs exist for each contrast.")
sum_dt <- rbindlist(all_sum, use.names = TRUE, fill = TRUE)


sum_dt[, subgroup := factor(
  paste0(ifelse(class == "hypo", "Hypo", "Hyper"), "-", timepoint),
  levels = c("Hypo-1h","Hypo-6h","Hypo-24h",
             "Hyper-1h","Hyper-6h","Hyper-24h"))
]


sum_dt[, class          := factor(class, levels = c("hypo","hyper"))]
sum_dt[, timepoint      := factor(as.character(timepoint),
                                  levels = c("1h","6h","24h"))]
sum_dt[, feature_category := factor(as.character(feature_category),
                                    levels = categories)]


plot_tbl <- sum_dt[, .(feature_category, class, timepoint, subgroup,
                       median_signed, IQR25_signed, IQR75_signed)]

tsv_out  <- file.path(out_root, "Combined_q05_meth5_SingleFeature_CGI_medianIQR_signed.tsv")
xlsx_out <- file.path(out_root, "Combined_q05_meth5_SingleFeature_CGI_medianIQR_signed.xlsx")

fwrite(plot_tbl, tsv_out, sep = "\t", quote = FALSE, na = "NA")

wb <- createWorkbook()
addWorksheet(wb, "median_IQR_signed")
writeData(wb, "median_IQR_signed", as.data.frame(plot_tbl))
saveWorkbook(wb, xlsx_out, overwrite = TRUE)


p <- ggplot(
  sum_dt,
  aes(x     = feature_category,
      y     = median_signed,
      ymin  = IQR25_signed,
      ymax  = IQR75_signed,
      color = class,
      shape = timepoint,
      group = subgroup)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_pointrange(position = position_dodge(width = 0.9), size = 0.6) +
  scale_color_manual(
    name   = "Class",
    values = c(hypo = "#F8766D",
               hyper = "#00BFC4"),
    breaks = c("hyper","hypo"),
    labels = c("Hyper (Δ > 0)", "Hypo (Δ < 0)")
  ) +
  scale_shape_manual(
    name   = "Timepoint",
    values = c(`1h` = 16, `6h` = 15, `24h` = 17),
    breaks = c("1h","6h","24h")
  ) +
  labs(
    title = "Δmethylation q05 (meth_diff ≥ 5): median Δ with IQR by category\nSingleFeature DMRs + CpG islands",
    x     = "Category",
    y     = "Median Δmethylation (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x       = element_text(angle = 25, hjust = 1),
    panel.grid.minor  = element_blank()
  )

pdf_out <- file.path(out_root, "Combined_q05_meth5_SingleFeature_CGI_medianIQR_signed.pdf")
png_out <- file.path(out_root, "Combined_q05_meth5_SingleFeature_CGI_medianIQR_signed.png")

ggsave(pdf_out, plot = p, width = 11, height = 5)
ggsave(png_out, plot = p, width = 11, height = 5, dpi = 300)

cat("Wrote:\n  ", tsv_out, "\n  ", xlsx_out, "\n  ",
    pdf_out, "\n  ", png_out, "\n", sep = "")


