#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


library(methylKit)


setwd(gb_path("wgbs/Reports"))


dmc_output_dir <- "../DMC_24M_vs_24hrKCl_analysis_cov4"
dir.create(dmc_output_dir, showWarnings = FALSE, recursive = TRUE)


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


cat("Reading methylation data...\n")
myobj <- methRead(file.list, sample.id = sample.id, assembly = "mm10",
                  treatment = treatment, pipeline = "bismarkCytosineReport",
                  mincov = 4)


myobj <- filterByCoverage(myobj, hi.perc = 99.9)


cat("Normalizing coverage...\n")
myobj <- normalizeCoverage(myobj)


cat("Merging methylation data...\n")
meth <- unite(myobj, destrand = FALSE, min.per.group = 2L)


cat("Performing differential methylation analysis...\n")
dmcs <- calculateDiffMeth(meth, test = "Chisq", adjust = "SLIM",
                          effect = "wmean", mc.cores = 4, save.db = FALSE)


write.table(getData(dmcs),
            file.path(dmc_output_dir, "all_DMCs_24M_vs_24hrKCl.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


write.table(getData(getMethylDiff(dmcs, difference = 5, qvalue = 0.1)),
            file.path(dmc_output_dir, "significant_DMCs_24M_vs_24hrKCl_q0.1_diff5.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)


cat("Applying p-value and methylation-difference filters (p ≤ 0.001, Δmeth ≥ 5%)...\n")
diffData <- getData(dmcs)
diffData_all   <- subset(diffData, abs(meth.diff) >= 5 & pvalue <= 0.001)
diffData_hyper <- subset(diffData, meth.diff >= 5 & pvalue <= 0.001)
diffData_hypo  <- subset(diffData, meth.diff <= -5 & pvalue <= 0.001)

write.table(diffData_all,
            file.path(dmc_output_dir, "manual_filtered_DMCs_p0001_diff5_all.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(diffData_hyper,
            file.path(dmc_output_dir, "manual_filtered_DMCs_p0001_diff5_hyper.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(diffData_hypo,
            file.path(dmc_output_dir, "manual_filtered_DMCs_p0001_diff5_hypo.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Filtering complete: ",
    "\n   Total:", nrow(diffData_all),
    "\n   Hypermethylated:", nrow(diffData_hyper),
    "\n   Hypomethylated:", nrow(diffData_hypo), "\n")


cat("Generating annotated DMC summary CSV...\n")
all_data <- diffData
all_data$direction <- ifelse(all_data$meth.diff > 0, "hypermethylated",
                             ifelse(all_data$meth.diff < 0, "hypomethylated", "no_change"))


meth_percent <- percMethylation(meth)
meth_percent <- as.data.frame(meth_percent)
meth_percent_t <- t(meth_percent)
control_means <- meth_percent_t[treatment == 0, , drop = FALSE]
treated_means <- meth_percent_t[treatment == 1, , drop = FALSE]

try({
  all_data$mean_control <- colMeans(control_means, na.rm = TRUE)
  all_data$mean_treated <- colMeans(treated_means, na.rm = TRUE)
  all_data$median_control <- apply(control_means, 2, median, na.rm = TRUE)
  all_data$median_treated <- apply(treated_means, 2, median, na.rm = TRUE)
}, silent = TRUE)

write.csv(all_data,
          file.path(dmc_output_dir, "all_CpGs_DMC_summary_annotated.csv"),
          row.names = FALSE)


cat("Summarizing DMC direction counts by filter...\n")

manual_counts <- data.frame(
  threshold = "p≤0.001 & Δmeth≥5",
  hypomethylated = nrow(diffData_hypo),
  hypermethylated = nrow(diffData_hyper)
)

dmc_q <- getData(getMethylDiff(dmcs, difference = 5, qvalue = 0.1))
q_hypo  <- sum(dmc_q$meth.diff <= -5)
q_hyper <- sum(dmc_q$meth.diff >= 5)

slim_counts <- data.frame(
  threshold = "q<0.1 & Δmeth≥5",
  hypomethylated = q_hypo,
  hypermethylated = q_hyper
)

direction_summary <- rbind(manual_counts, slim_counts)
write.csv(direction_summary,
          file.path(dmc_output_dir, "DMC_direction_summary_counts.csv"),
          row.names = FALSE)

cat("DMC analysis complete for 24M vs 24K. Summary CSVs written.\n")
