#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



library(methylKit)


setwd(gb_path("wgbs/Reports"))


output_dir <- "../methylKit_bedgraph_bigwig"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


file.list <- list(
  "A0_CKDL240037499-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A1K_CKDL240037501-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A1M_CKDL240037500-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A24K_CKDL240037505-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A24M_CKDL240037504-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A6K_CKDL240037503-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "A6M_CKDL240037502-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B0_CKDL240037506-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1K_CKDL240037508-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1M_CKDL240037507-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24K_CKDL240037512-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24M_CKDL240037511-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "B6K_CKDL240037510-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "B6M_CKDL240037509-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C0_CKDL240037513-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1K_CKDL240037515-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1M_CKDL240037514-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24K_CKDL240037519-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24M_CKDL240037518-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C6K_CKDL240037517-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C6M_CKDL240037516-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"
)


sample.id <- list(
  "A0", "A1K", "A1M", "A24K", "A24M", "A6K", "A6M",
  "B0", "B1K", "B1M", "B24K", "B24M", "B6K", "B6M",
  "C0", "C1K", "C1M", "C24K", "C24M", "C6K", "C6M"
)


treatment <- rep(0, length(sample.id))


myobj <- methRead(
  file.list,
  sample.id = sample.id,
  assembly = "mm10",
  treatment = treatment,
  pipeline = "bismarkCytosineReport",
  mincov = 4
)


myobj <- filterByCoverage(myobj, hi.perc = 99.9)


myobj <- normalizeCoverage(myobj)


for (i in 1:length(myobj)) {
  out_file <- file.path(output_dir, paste0(sample.id[[i]], "_normalized_cov4.bedGraph"))

  bedgraph(myobj[[i]],
           file.name = out_file,
           col.name = "perc.meth")

  cat("Saved bedGraph:", out_file, "\n")
}

cat("BedGraph generation complete for all samples (normalized, cov ≥ 4)\n")
