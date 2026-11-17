#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(annotatr)
  library(GenomicRanges)
  library(data.table)
  library(openxlsx)
})


base_dirs <- c(
  gb_path("wgbs/methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp")
)


merged_subdir <- "Merged_q05_meth5"

merged_pat <- "^[0-9]+M_vs_[0-9]+K_Merged_DMRs_cov4_(hyper|hypo)\\.txt$"

out_subdir    <- file.path(merged_subdir, "SingleFeature_q05_meth5")


cgi_annot <- build_annotations(
  genome = "mm10",
  annotations = c("mm10_cpg_islands")
)


contrast_name <- function(bdir) sub("^DMR_", "", sub("_multiwindow$", "", basename(dirname(bdir))))
tp_from_contrast <- function(cn) {
  if (grepl("^1M_", cn))  return("1h")
  if (grepl("^6M_", cn))  return("6h")
  if (grepl("^24M_", cn)) return("24h")
  return(cn)
}

read_q05 <- function(path) {
  dt <- fread(path)
  if ("meth.diff" %in% names(dt)) setnames(dt, "meth.diff", "meth_diff")
  req <- c("chr","start","end","qvalue","meth_diff")
  if (!all(req %in% names(dt))) stop("Missing columns in ", path, " (need: ", paste(req, collapse=","), ")")
  dt[]
}

annotate_cgi <- function(dmrs_dt) {
  gr <- GRanges(seqnames = dmrs_dt$chr, ranges = IRanges(dmrs_dt$start, dmrs_dt$end))
  ann <- annotate_regions(gr, annotations = cgi_annot, ignore.strand = TRUE)
  as.data.table(ann)
}

summarize_dt <- function(dt) {

  dt[, .(
    N             = .N,
    median_abs    = as.numeric(median(abs_meth_diff, na.rm = TRUE)),
    IQR25_abs     = as.numeric(quantile(abs_meth_diff, 0.25, na.rm = TRUE)),
    IQR75_abs     = as.numeric(quantile(abs_meth_diff, 0.75, na.rm = TRUE)),
    median_signed = as.numeric(median(meth_diff, na.rm = TRUE))
  ), by = .(contrast, class)][order(class)]
}


run_contrast <- function(base_dir) {
  comp <- contrast_name(base_dir)
  tp   <- tp_from_contrast(comp)


  out_dir <- file.path(base_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  merged_dir <- file.path(base_dir, merged_subdir)
  files <- list.files(merged_dir, pattern = merged_pat, full.names = TRUE)
  if (!length(files)) {
    cat(" No q05/meth5 merged files in:", merged_dir, "\n")
    return(invisible(NULL))
  }

  long_list <- lapply(files, function(f) {
    class <- if (grepl("_hyper\\.txt$", f)) "hyper" else "hypo"
    dmrs  <- read_q05(f)
    if (!nrow(dmrs)) return(NULL)

    ann_dt <- annotate_cgi(dmrs)
    if (!nrow(ann_dt)) return(NULL)


    setnames(ann_dt, old = c("seqnames","start","end"), new = c("chr","start","end"))


    dmrs[, dmr_id := paste(chr, start, end, sep = "_")]
    ann_dt[, dmr_id := paste(chr, start, end, sep = "_")]
    ann_dt <- merge(
      ann_dt,
      unique(dmrs[, .(dmr_id, qvalue, meth_diff)]),
      by = "dmr_id", all.x = TRUE, sort = FALSE
    )


    ann_dt <- ann_dt[!duplicated(dmr_id)]

    ann_dt[, `:=`(
      contrast      = comp,
      timepoint     = tp,
      class         = class,
      cgi_category  = "CpG island",
      abs_meth_diff = abs(meth_diff),
      source        = basename(f)
    )]


    ann_dt[, .(contrast, timepoint, class, cgi_category,
               chr, start, end, qvalue, meth_diff, abs_meth_diff, dmr_id, source)]
  })

  long_dt <- rbindlist(long_list, use.names = TRUE, fill = TRUE)
  if (!nrow(long_dt)) {
    cat(" No CpG island overlaps for:", comp, "\n")
    return(invisible(NULL))
  }


  summary_dt <- summarize_dt(long_dt)


  long_out    <- file.path(out_dir, paste0(comp, "_cgi_methylation_long.tsv"))
  summary_out <- file.path(out_dir, paste0(comp, "_cgi_methylation_summary.tsv"))
  xlsx_out    <- file.path(out_dir, paste0(comp, "_cgi_methylation.xlsx"))

  fwrite(long_dt,    long_out,    sep = "\t", quote = FALSE, na = "NA")
  fwrite(summary_dt, summary_out, sep = "\t", quote = FALSE, na = "NA")

  wb <- createWorkbook()
  addWorksheet(wb, "long");    writeData(wb, "long",    as.data.frame(long_dt))
  addWorksheet(wb, "summary"); writeData(wb, "summary", as.data.frame(summary_dt))
  saveWorkbook(wb, xlsx_out, overwrite = TRUE)

  cat("", comp, " →\n  ", long_out, "\n  ", summary_out, "\n  ", xlsx_out, "\n", sep = "")
  invisible(TRUE)
}


for (bd in base_dirs) run_contrast(bd)
cat("DONE — CpG ISLAND tables & summaries written into SingleFeature_q05_meth5/.\n")
