#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
PROJECT_ROOT <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(PROJECT_ROOT, "R", "utils_figure_paths.R"))

run_id <- "2026-09-19_p2609_salmonte"
run_dir <- file.path(PROJECT_ROOT, "runs", run_id)
input_dir <- file.path(PROJECT_ROOT, "inputs", "P2609-01_fastq")
figure_dir <- revision_figure_dir(run_id)
table_dir <- file.path(run_dir, "tables")
dir.create(table_dir, showWarnings = FALSE)

samples <- read.delim(file.path(input_dir, "samples.tsv"), stringsAsFactors = FALSE,
                      check.names = FALSE)
qc <- read.delim(file.path(run_dir, "QC_SUMMARY.tsv"), stringsAsFactors = FALSE)
if (nrow(samples) != 9L || nrow(qc) != 9L || any(qc$status != "complete"))
  stop("Expected nine completed and uniquely labeled quantifications")
if (anyDuplicated(samples$sample_id) || anyDuplicated(qc$sample_id) ||
    !setequal(samples$sample_id, qc$sample_id)) stop("Sample IDs do not match")
samples <- samples[match(qc$sample_id, samples$sample_id), ]
condition <- factor(samples$condition, levels = c("TTX", "Bic", "Bic_Nimo"))
if (!identical(as.integer(table(condition)), c(3L, 3L, 3L)))
  stop("Expected three samples in each condition")
if (!all(qc$target_count == 301L) || !all(qc$library_types == "IU"))
  stop("Salmon target count or inferred library type differs between libraries")

quant_files <- file.path(run_dir, "quant", qc$sample_id, "quant.sf")
quants <- lapply(quant_files, function(path) read.delim(path, check.names = FALSE,
                                                       stringsAsFactors = FALSE))
targets <- quants[[1L]]$Name
if (length(targets) != 301L || anyDuplicated(targets) ||
    !all(vapply(quants, function(x) identical(x$Name, targets), logical(1))))
  stop("The nine Salmon quant.sf target lists are not identical")

counts <- do.call(cbind, lapply(quants, function(x) x$NumReads))
storage.mode(counts) <- "double"
rownames(counts) <- targets
colnames(counts) <- qc$sample_id
if (any(!is.finite(counts)) || any(counts < 0) ||
    any(!is.finite(qc$num_processed)) || any(qc$num_processed <= 0))
  stop("Non-finite or negative input values")

cpm_input <- sweep(counts, 2, qc$num_processed, "/") * 1e6
log_cpm <- log2(cpm_input + 1)
expressed <- rowSums(counts >= 10) >= 3
if (sum(expressed) < 3L) stop("Too few expressed TE targets for PCA")

pca <- prcomp(t(log_cpm[expressed, , drop = FALSE]), center = TRUE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
scores <- data.frame(sample_id = qc$sample_id, condition = condition,
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     processed_pairs = qc$num_processed,
                     TE_mapped_pairs = qc$num_mapped,
                     TE_mapped_percent = qc$percent_mapped,
                     stringsAsFactors = FALSE)
write.table(scores, file.path(table_dir, "PCA_SCORES.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

group_mean <- sapply(levels(condition), function(group) {
  rowMeans(cpm_input[, condition == group, drop = FALSE])
})
family <- data.frame(Name = targets,
                     mean_CPM_TTX = group_mean[, "TTX"],
                     mean_CPM_Bic = group_mean[, "Bic"],
                     mean_CPM_Bic_Nimo = group_mean[, "Bic_Nimo"],
                     log2_ratio_Bic_TTX = log2((group_mean[, "Bic"] + 1) /
                                           (group_mean[, "TTX"] + 1)),
                     log2_ratio_Bic_Nimo_TTX = log2((group_mean[, "Bic_Nimo"] + 1) /
                                                (group_mean[, "TTX"] + 1)),
                     log2_ratio_Bic_Nimo_Bic = log2((group_mean[, "Bic_Nimo"] + 1) /
                                                (group_mean[, "Bic"] + 1)),
                     expressed_for_PCA = expressed,
                     stringsAsFactors = FALSE)
variation <- apply(log_cpm, 1, var)
top <- names(sort(variation[expressed], decreasing = TRUE))[seq_len(min(30L, sum(expressed)))]
family$top30_heatmap <- family$Name %in% top
write.table(family, file.path(table_dir, "FAMILY_SUMMARY.tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)

palette <- c(TTX = "#315F9A", Bic = "#C15B3B", Bic_Nimo = "#639A64")
labels <- c(TTX = "TTX", Bic = "Bic", Bic_Nimo = "Bic + Nimo")
pca_plot <- ggplot(scores, aes(PC1, PC2, color = condition, shape = condition)) +
  geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample_id), color = "black", size = 3,
                           box.padding = 0.25, point.padding = 0.12,
                           max.overlaps = Inf, min.segment.length = 0,
                           seed = 20260919, show.legend = FALSE) +
  scale_color_manual(values = palette, labels = labels) +
  scale_shape_manual(values = c(TTX = 16, Bic = 17, Bic_Nimo = 15), labels = labels) +
  labs(title = "SalmonTE expression across P2609-01 libraries",
       subtitle = sprintf("%d of 301 TE targets; >=10 assigned fragments in >=3 samples", sum(expressed)),
       x = sprintf("PC1 (%.1f%%)", 100 * variance[1]),
       y = sprintf("PC2 (%.1f%%)", 100 * variance[2]),
       color = "Condition", shape = "Condition",
       caption = "log2(TE fragments per million input pairs + 1); centered, not scaled.") +
  scale_x_continuous(expand = expansion(mult = c(0.13, 0.23))) +
  scale_y_continuous(expand = expansion(mult = c(0.13, 0.13))) +
  coord_fixed(clip = "off") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        plot.caption = element_text(size = 8, hjust = 0),
        legend.position = "bottom", plot.margin = margin(10, 18, 10, 10))

pdf(file.path(figure_dir, "P2609_SalmonTE_PCA.pdf"), width = 8, height = 5.8,
    useDingbats = FALSE)
print(pca_plot)
dev.off()
png(file.path(figure_dir, "P2609_SalmonTE_PCA.png"), width = 1920, height = 1392,
    res = 240, type = "cairo")
print(pca_plot)
dev.off()


heat <- log_cpm[top, , drop = FALSE]
heat <- sweep(heat, 1, rowMeans(heat), "-")
limit <- max(0.5, ceiling(max(abs(heat)) * 10) / 10)
colors <- colorRampPalette(c("#315F9A", "white", "#B84A40"))(101)
annotation <- data.frame(Condition = factor(labels[as.character(condition)],
                                            levels = unname(labels)))
rownames(annotation) <- qc$sample_id
heat_args <- list(mat = heat, color = colors, breaks = seq(-limit, limit, length.out = 102),
                  cluster_rows = TRUE, cluster_cols = FALSE, gaps_col = c(3, 6),
                  annotation_col = annotation,
                  annotation_colors = list(Condition = setNames(unname(palette), unname(labels))),
                  border_color = NA, fontsize_row = 8, fontsize_col = 9,
                  main = "30 most variable expressed TE targets\nrow-centered log2(input-pair CPM + 1)")
pdf(file.path(figure_dir, "P2609_SalmonTE_family_heatmap.pdf"), width = 8.3,
    height = 8.7, useDingbats = FALSE)
do.call(pheatmap::pheatmap, heat_args)
dev.off()
png(file.path(figure_dir, "P2609_SalmonTE_family_heatmap.png"),
    width = 1826, height = 1914, res = 220, type = "cairo")
do.call(pheatmap::pheatmap, heat_args)
dev.off()

cat(sprintf("PCA: %d expressed targets, PC1 %.1f%%, PC2 %.1f%%\n",
            sum(expressed), 100 * variance[1], 100 * variance[2]))
cat("Figures:", figure_dir, "\n")
cat("Tables:", table_dir, "\n")
