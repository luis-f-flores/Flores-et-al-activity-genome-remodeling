#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



library(methylKit)
library(ggplot2)
library(ggrepel)


setwd(gb_path("wgbs/Reports"))


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


cat("Reading samples with mincov = 4...\n")
myobj <- methRead(file.list,
                  sample.id = sample.id,
                  assembly = "mm10",
                  treatment = treatment,
                  pipeline = "bismarkCytosineReport",
                  mincov = 4)


cat("Filtering high-coverage outliers...\n")
myobj <- filterByCoverage(myobj, hi.perc = 99.9)


cat("Normalizing coverage...\n")
myobj <- normalizeCoverage(myobj)


cat("Uniting samples...\n")
meth <- unite(myobj, destrand = FALSE)


qc_dir <- "methylKit_QC"
dir.create(qc_dir, showWarnings = FALSE)


pca_obj <- PCASamples(meth, obj.return = TRUE)

if (!is.null(pca_obj)) {
  cat("PCA complete — generating the ggplot2 PDF...\n")


  pca_var <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)
  labels_df <- data.frame(pca_obj$x[,1:2], Sample = meth@sample.ids)
  write.csv(labels_df, file.path(qc_dir, "PCA_samples_coordinates.csv"), row.names = FALSE)
  write.csv(data.frame(component = c("PC1", "PC2"), percent = pca_var), file.path(qc_dir, "PCA_variance.csv"), row.names = FALSE)


  p <- ggplot(labels_df, aes(PC1, PC2, label = Sample)) +
    geom_point(color = "darkblue", size = 3) +
    geom_text_repel(max.overlaps = 50, size = 3, box.padding = 0.3, point.padding = 0.2) +
    labs(
      title = "PCA of CpG Methylation — All Samples",
      x = paste0("PC1 (", pca_var[1], "%)"),
      y = paste0("PC2 (", pca_var[2], "%)")
    ) +
    theme_minimal(base_size = 14) +
    expand_limits(x = max(labels_df$PC1) + 500)


  cat(sprintf("PC1 range: %.1f to %.1f\n", min(labels_df$PC1), max(labels_df$PC1)))
  cat(sprintf("PC2 range: %.1f to %.1f\n", min(labels_df$PC2), max(labels_df$PC2)))
  cat("Writing PCA plot to file...\n")


  pdf(file.path(qc_dir, "PCA_samples.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
  cat("PCA_samples.pdf written successfully.\n")

} else {
  cat("PCA failed — check sample overlap or data integrity.\n")
}

cat("PCA outputs written to", qc_dir, "\n")
