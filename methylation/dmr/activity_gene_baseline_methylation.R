#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


source("config.R")
source("utils_annotation_helpers.R")

suppressPackageStartupMessages({
  library(methylKit)
  library(ggplot2)
  library(dplyr)
})

gene_models <- load_collapsed_gene_models(gene_body_bed_file, tss_bed_file)
promoters_gr <- build_promoter_windows(gene_models$tss, window = PROMOTER_WINDOW)

baseline_files <- raw_cpg_report_files[["0hr"]]
cat("Reading 3 baseline (0hr) raw Bismark CpG reports...\n")
myobj <- methRead(as.list(baseline_files),
                   sample.id = as.list(names(baseline_files)),
                   assembly = "mm10",
                   treatment = c(0, 0, 0),
                   pipeline = "bismarkCytosineReport",
                   mincov = 4)

myobj <- filterByCoverage(myobj, hi.perc = 99.9)

cat("Aggregating per-CpG counts into", length(promoters_gr), "promoter regions...\n")
region_obj <- regionCounts(myobj, regions = promoters_gr)

cat("Uniting samples...\n")
meth_regions <- unite(region_obj, destrand = FALSE, min.per.group = 1L)

pm <- percMethylation(meth_regions)
mean_pct <- rowMeans(pm, na.rm = TRUE)


region_coords <- getData(meth_regions)[, c("chr", "start", "end")]
region_gr <- GenomicRanges::GRanges(region_coords$chr,
                                     IRanges::IRanges(region_coords$start, region_coords$end))
match_idx <- GenomicRanges::findOverlaps(region_gr, promoters_gr, type = "equal", select = "first")

per_gene_baseline <- data.frame(
  gene_id     = promoters_gr$gene_id[match_idx],
  gene_symbol = promoters_gr$gene_symbol[match_idx],
  mean_promoter_mCG = mean_pct
)
per_gene_baseline <- per_gene_baseline[stats::complete.cases(per_gene_baseline), ]

classify_gene <- function(symbol) {
  if (symbol %in% activity_gene_lists$dPRG) return("dPRG")
  if (symbol %in% activity_gene_lists$SRG)  return("SRG")
  if (symbol %in% activity_gene_lists$rPRG) return("rPRG")
  "Other"
}

per_gene_baseline$gene_class <- vapply(per_gene_baseline$gene_symbol, classify_gene, character(1))
per_gene_baseline$gene_class <- factor(per_gene_baseline$gene_class,
                                        levels = c("rPRG", "dPRG", "SRG", "Other"))

background <- per_gene_baseline$mean_promoter_mCG[per_gene_baseline$gene_class == "Other"]

stat_summary <- lapply(c("dPRG", "SRG", "rPRG"), function(cls) {
  vals <- per_gene_baseline$mean_promoter_mCG[per_gene_baseline$gene_class == cls]
  if (length(vals) < 3) {
    return(data.frame(gene_class = cls, n = length(vals), median_mCG = NA, p_value = NA))
  }
  wt <- stats::wilcox.test(vals, background)
  data.frame(
    gene_class = cls,
    n = length(vals),
    median_mCG = stats::median(vals, na.rm = TRUE),
    p_value = wt$p.value
  )
})
stat_summary <- do.call(rbind, stat_summary)
stat_summary <- rbind(
  stat_summary,
  data.frame(gene_class = "Other", n = length(background),
             median_mCG = stats::median(background, na.rm = TRUE), p_value = NA)
)

write.csv(stat_summary,
          file.path(OUTPUT_TABLES_DIR, "01_baseline_methylation_by_class.csv"),
          row.names = FALSE)

p <- ggplot(per_gene_baseline, aes(x = gene_class, y = mean_promoter_mCG, fill = gene_class)) +

  geom_violin(trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(title = "Baseline (0hr) promoter methylation by activity gene class",
       subtitle = "methylKit regionCounts() on raw A0/B0/C0 CpG reports, -1kb/+100bp promoter",
       x = NULL, y = "Mean promoter mCG (%)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_FIGURES_DIR, "01_baseline_methylation_by_class.pdf"), p, width = 5, height = 4)

cat("Baseline promoter methylation complete. Genes by class:\n")
print(table(per_gene_baseline$gene_class))
cat("Output in", OUTPUT_TABLES_DIR, "and", OUTPUT_FIGURES_DIR, "\n")
