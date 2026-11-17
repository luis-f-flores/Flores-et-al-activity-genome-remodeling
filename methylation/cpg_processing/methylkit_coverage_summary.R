#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}

library(methylKit)


setwd(gb_path("wgbs/Reports"))


output_dir <- gb_path("wgbs/methylKit_analysis")
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


sample.id <- list("A0", "A1K", "A1M", "A24K", "A24M", "A6K", "A6M",
                  "B0", "B1K", "B1M", "B24K", "B24M", "B6K", "B6M",
                  "C0", "C1K", "C1M", "C24K", "C24M", "C6K", "C6M")


treatment <- c(0, 1, 1, 1, 1, 1, 1,
               0, 1, 1, 1, 1, 1, 1,
               0, 1, 1, 1, 1, 1, 1)


myobj_list <- list()
failed_samples <- c()

cat("Loading each sample individually with mincov = 4...\n")
for (i in seq_along(file.list)) {
  tryCatch({
    this_obj <- methRead(
      file.list[[i]],
      sample.id = sample.id[[i]],
      assembly = "mm10",
      treatment = treatment[[i]],
      pipeline = "bismarkCytosineReport",
      mincov = 4
    )
    myobj_list[[sample.id[[i]]]] <- this_obj
    cat(paste0("Loaded: ", sample.id[[i]], "\n"))
  }, error = function(e) {
    cat(paste0("Failed: ", sample.id[[i]], " — ", conditionMessage(e), "\n"))
    failed_samples <<- c(failed_samples, sample.id[[i]])
  })
}
myobj <- new("methylRawList", myobj_list)

if (length(failed_samples) > 0) {
  writeLines(failed_samples, file.path(output_dir, "failed_samples.txt"))
  cat("Failed samples written to failed_samples.txt\n")
}


myobj <- filterByCoverage(myobj, hi.perc = 99.9)


cat("Normalizing coverage...\n")
myobj <- normalizeCoverage(myobj)


coverage_summary <- data.frame()
for (i in seq_along(myobj)) {
  cov_vals <- myobj[[i]]$coverage
  coverage_summary <- rbind(coverage_summary, data.frame(
    Sample = myobj[[i]]@sample.id,
    Min_Coverage = min(cov_vals),
    Max_Coverage = max(cov_vals),
    Mean_Coverage = mean(cov_vals),
    Median_Coverage = median(cov_vals),
    Percent_5X = mean(cov_vals >= 5) * 100,
    Percent_10X = mean(cov_vals >= 10) * 100,
    Percent_20X = mean(cov_vals >= 20) * 100
  ))
}
write.csv(coverage_summary, file.path(output_dir, "coverage_summary.csv"), row.names=FALSE)


cat("Coverage summary written.\n")
