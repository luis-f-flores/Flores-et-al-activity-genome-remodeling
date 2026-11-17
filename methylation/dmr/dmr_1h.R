#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



library(methylKit)


setwd(gb_path("wgbs/Reports"))


output_base <- "../methylKit_DMR/DMR_1M_vs_1K_multiwindow"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)


file.list <- list(
  "A1M_CKDL240037500-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1M_CKDL240037507-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1M_CKDL240037514-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A1K_CKDL240037501-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1K_CKDL240037508-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1K_CKDL240037515-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"
)
sample.id <- list("A1M", "B1M", "C1M", "A1K", "B1K", "C1K")
treatment <- c(0, 0, 0, 1, 1, 1)


myobj <- methRead(file.list, sample.id = sample.id, assembly = "mm10", treatment = treatment,
                  pipeline = "bismarkCytosineReport", mincov = 4)
myobj <- filterByCoverage(myobj, hi.perc = 99.9)
myobj <- normalizeCoverage(myobj)


window_configs <- list(
  list(win = 500, step = 250, cov = 3, label = "500bp")


)


diff_cutoff <- 5
q_cutoff <- 0.05


for (cfg in window_configs) {
  cat("Running DMR analysis for window:", cfg$label, "\n")

  cfg_dir <- file.path(output_base, paste0("DMR_win_", cfg$label))
  dir.create(cfg_dir, showWarnings = FALSE)

  tiles <- tileMethylCounts(myobj, win.size = cfg$win, step.size = cfg$step, cov.bases = cfg$cov)
  meth_tiled <- unite(tiles, destrand = TRUE, min.per.group = 2L)
  myDiff <- calculateDiffMeth(meth_tiled)
  diffData <- getData(myDiff)


  write.table(diffData,
              file.path(cfg_dir, paste0("DMRs_all_unfiltered_win_", cfg$label, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)


  diff_label <- paste0("diff", diff_cutoff, "_q005")
  q_dmrs   <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = q_cutoff)
  q_hyper  <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = q_cutoff, type = "hyper")
  q_hypo   <- getMethylDiff(myDiff, difference = diff_cutoff, qvalue = q_cutoff, type = "hypo")

  write.table(q_dmrs,  file.path(cfg_dir, paste0("DMRs_", diff_label, ".txt")),       sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(q_hyper, file.path(cfg_dir, paste0("DMRs_", diff_label, "_hyper.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(q_hypo,  file.path(cfg_dir, paste0("DMRs_", diff_label, "_hypo.txt")),  sep = "\t", quote = FALSE, row.names = FALSE)

  cat("", cfg$label, "-", diff_label,
      ": total =", nrow(q_dmrs),
      "| hyper =", nrow(q_hyper),
      "| hypo =", nrow(q_hypo), "\n")

  pdf(file.path(cfg_dir, "DMR_distribution.pdf"))
  diffMethPerChr(myDiff, plot = TRUE)
  dev.off()

  cat("Finished window:", cfg$label, "\n\n")
}

cat("1M vs 1K | DMR scan complete (500bp, q ≤ 0.05, diff ≥ 5%)\n")
