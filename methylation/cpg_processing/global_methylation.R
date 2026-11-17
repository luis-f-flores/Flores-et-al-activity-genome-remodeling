#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}

library(methylKit)
library(dplyr)


setwd(gb_path("wgbs/Reports"))


output_dir <- "../global_cpg_methylation_split_safe"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


file.list <- list(
  "A0_CKDL240037499-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B0_CKDL240037506-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "C0_CKDL240037513-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A1K_CKDL240037501-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1K_CKDL240037508-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1K_CKDL240037515-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A1M_CKDL240037500-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B1M_CKDL240037507-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "C1M_CKDL240037514-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A6K_CKDL240037503-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B6K_CKDL240037510-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C6K_CKDL240037517-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A6M_CKDL240037502-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B6M_CKDL240037509-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C6M_CKDL240037516-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A24K_CKDL240037505-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24K_CKDL240037512-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24K_CKDL240037519-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "A24M_CKDL240037504-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt",
  "B24M_CKDL240037511-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt",
  "C24M_CKDL240037518-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"
)


sample.id <- list(
  "A0", "B0", "C0", "A1K", "B1K", "C1K", "A1M", "B1M", "C1M",
  "A6K", "B6K", "C6K", "A6M", "B6M", "C6M",
  "A24K", "B24K", "C24K", "A24M", "B24M", "C24M"
)


conditions <- rep(c("0hr", "1K", "1M", "6K", "6M", "24K", "24M"), each = 3)


sample_info <- data.frame(
  File = unlist(file.list),
  Sample = unlist(sample.id),
  Condition = conditions,
  stringsAsFactors = FALSE
)


all_samples <- list()
all_conditions <- list()

for (cond in unique(sample_info$Condition)) {
  cat("\nProcessing:", cond, "\n")

  group <- sample_info[sample_info$Condition == cond, ]

  myobj <- methRead(
    location = as.list(group$File),
    sample.id = as.list(group$Sample),
    assembly = "mm10",
    treatment = rep(1, 3),
    pipeline = "bismarkCytosineReport",
    mincov = 4
  )

  myobj <- filterByCoverage(myobj, hi.perc = 99.9)
  myobj <- normalizeCoverage(myobj)

  cat("Running unite()...\n")
  meth_united <- unite(myobj, destrand = FALSE, min.per.group = 2L)
  meth_matrix <- percMethylation(meth_united)

  sample_df <- data.frame(
    Sample = group$Sample,
    Condition = cond,
    Methylation_Percentage = colMeans(meth_matrix, na.rm = TRUE),
    Num_CpGs = colSums(!is.na(meth_matrix))
  )

  all_samples[[cond]] <- sample_df

  wmeth <- sum(sample_df$Num_CpGs * sample_df$Methylation_Percentage) / sum(sample_df$Num_CpGs)
  all_conditions[[cond]] <- data.frame(Condition = cond, Weighted_Methylation = wmeth)
}


final_sample_df <- do.call(rbind, all_samples)
final_condition_df <- do.call(rbind, all_conditions)

write.csv(final_sample_df, file.path(output_dir, "Sample_Methylation.csv"), row.names = FALSE)
write.csv(final_condition_df, file.path(output_dir, "Condition_Weighted_Methylation.csv"), row.names = FALSE)

cat("\nAll methylation outputs written to:", output_dir, "\n")
