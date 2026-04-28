#!/usr/bin/env Rscript
# 06_RNAseq_cluster_heatmap.R
#
# Plot a z-scored RNA-seq trajectory heatmap ordered by cluster.

suppressPackageStartupMessages({
    library(pheatmap)
    library(dplyr)
    library(readr)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
    idx <- match(flag, args)

    if (is.na(idx)) {
        return(default)
    }

    if (idx == length(args)) {
        stop("Missing value after ", flag)
    }

    args[[idx + 1]]
}

input_file <- get_arg("--input")
outdir <- get_arg("--outdir", "cluster_heatmap")
prefix <- get_arg("--prefix", "RNAseq_clustered_trajectory")
cluster_col <- get_arg("--cluster-col", "Cluster")
gene_label_col <- get_arg("--gene-label-col", "gene_name")

if (is.null(input_file) || !file.exists(input_file)) {
    stop("Usage: Rscript 06_RNAseq_cluster_heatmap.R --input <Clustered_Genes_Hierarchical_Euclidean_Z_zscores.csv>")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, check.names = FALSE)

if (!cluster_col %in% colnames(df)) {
    stop("Cluster column not found: ", cluster_col)
}

log2fc_cols <- grep("^log2FC_", colnames(df), value = TRUE)

if (length(log2fc_cols) < 2) {
    stop("Need at least two columns beginning with log2FC_")
}

df <- df %>%
    arrange(.data[[cluster_col]])

mat <- as.matrix(df[, log2fc_cols])
mode(mat) <- "numeric"

rownames(mat) <- if (gene_label_col %in% colnames(df)) {
    make.unique(as.character(df[[gene_label_col]]))
} else {
    make.unique(as.character(df$gene_id))
}

annotation_row <- data.frame(
    Cluster = factor(df[[cluster_col]])
)

rownames(annotation_row) <- rownames(mat)

pdf(
    file.path(outdir, paste0(prefix, "_heatmap.pdf")),
    width = 5.5,
    height = 8
)

pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    annotation_row = annotation_row,
    main = "RNA-seq DEG trajectory clusters"
)

dev.off()

png(
    file.path(outdir, paste0(prefix, "_heatmap.png")),
    width = 1800,
    height = 2600,
    res = 300
)

pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    annotation_row = annotation_row,
    main = "RNA-seq DEG trajectory clusters"
)

dev.off()

cluster_means <- df %>%
    group_by(.data[[cluster_col]]) %>%
    summarize(
        across(all_of(log2fc_cols), mean, na.rm = TRUE),
        n_genes = n(),
        .groups = "drop"
    )

write_csv(
    cluster_means,
    file.path(outdir, paste0(prefix, "_cluster_mean_zscores.csv"))
)

cluster_means_long <- tidyr::pivot_longer(
    cluster_means,
    cols = all_of(log2fc_cols),
    names_to = "timepoint",
    values_to = "mean_zscore"
)

cluster_means_long$timepoint <- gsub("^log2FC_", "", cluster_means_long$timepoint)

p <- ggplot(
    cluster_means_long,
    aes(
        x = timepoint,
        y = mean_zscore,
        group = .data[[cluster_col]],
        color = as.factor(.data[[cluster_col]])
    )
) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    labs(
        x = "Timepoint",
        y = "Mean z-scored log2FC",
        color = "Cluster"
    ) +
    theme_minimal(base_size = 12)

ggsave(
    file.path(outdir, paste0(prefix, "_cluster_mean_trajectories.pdf")),
    p,
    width = 6,
    height = 4
)

ggsave(
    file.path(outdir, paste0(prefix, "_cluster_mean_trajectories.png")),
    p,
    width = 6,
    height = 4,
    dpi = 300
)

message("Wrote heatmap and cluster mean trajectories to: ", normalizePath(outdir))
