#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



library(methylKit)


setwd(gb_path("wgbs/Reports"))


output_base <- "../methylKit_DMR/DMR_24M_vs_24K_multiwindow"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)


file.list <- list(
  "A24M_CKDL240037504-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24M_CKDL240037511-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24M_CKDL240037518-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A24K_CKDL240037505-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24K_CKDL240037512-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24K_CKDL240037519-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"
)
sample.id <- list("A24M", "B24M", "C24M", "A24K", "B24K", "C24K")
treatment <- c(0, 0, 0, 1, 1, 1)


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


window_configs <- list(
  list(win = 500, step = 250, cov = 3, label = "500bp")


)


thresholds <- list(
  list(diff = 5)


)


for (cfg in window_configs) {
  cat("Running DMR analysis for window:", cfg$label, "\n")

  cfg_dir <- file.path(output_base, paste0("DMR_win_", cfg$label))
  dir.create(cfg_dir, showWarnings = FALSE)

  tiles <- tileMethylCounts(myobj,
                            win.size = cfg$win,
                            step.size = cfg$step,
                            cov.bases = cfg$cov)

  meth_tiled <- unite(tiles, destrand = TRUE, min.per.group = 2L)
  myDiff <- calculateDiffMeth(meth_tiled)
  diffData <- getData(myDiff)


  write.table(diffData,
              file.path(cfg_dir, paste0("DMRs_all_unfiltered_win_", cfg$label, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)


  for (th in thresholds) {
    diff_label <- paste0("diff", th$diff, "_q005")

    q_dmrs   <- getMethylDiff(myDiff, difference = th$diff, qvalue = 0.05)
    q_hyper  <- getMethylDiff(myDiff, difference = th$diff, qvalue = 0.05, type = "hyper")
    q_hypo   <- getMethylDiff(myDiff, difference = th$diff, qvalue = 0.05, type = "hypo")

    write.table(q_dmrs,  file.path(cfg_dir, paste0("DMRs_", diff_label, ".txt")),       sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(q_hyper, file.path(cfg_dir, paste0("DMRs_", diff_label, "_hyper.txt")), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(q_hypo,  file.path(cfg_dir, paste0("DMRs_", diff_label, "_hypo.txt")),  sep = "\t", quote = FALSE, row.names = FALSE)

    cat("", cfg$label, "-", diff_label, ": q ≤ 0.05 =", nrow(q_dmrs), "\n")
  }

  pdf(file.path(cfg_dir, "DMR_distribution.pdf"))
  diffMethPerChr(myDiff, plot = TRUE)
  dev.off()

  cat("Finished window:", cfg$label, "\n\n")
}

cat("24M vs 24K DMR scan complete with 5% diff and q ≤ 0.05!\n")
