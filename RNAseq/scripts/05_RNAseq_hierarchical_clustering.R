#!/usr/bin/env Rscript
# 05_RNAseq_hierarchical_clustering.R
#
# Z-score per-gene RNA-seq trajectories and perform Euclidean hierarchical clustering.

suppressPackageStartupMessages({
    library(dplyr)
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
outdir <- get_arg("--outdir", "clustering_results")
k <- as.integer(get_arg("--k", "5"))
method <- get_arg("--method", "ward.D2")
distance <- get_arg("--distance", "euclidean")

if (is.null(input_file) || !file.exists(input_file)) {
    stop("Usage: Rscript 05_RNAseq_hierarchical_clustering.R --input <combined_log2FC_FPKM_GlobalDEG.csv>")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, check.names = FALSE)

if (!"gene_id" %in% colnames(df)) {
    colnames(df)[1] <- "gene_id"
}

log2fc_cols <- grep("^log2FC_", colnames(df), value = TRUE)

if (length(log2fc_cols) < 2) {
    stop("Need at least two columns beginning with log2FC_")
}

ts_matrix_raw <- as.matrix(df[, log2fc_cols])
mode(ts_matrix_raw) <- "numeric"

ts_matrix_z <- t(scale(t(ts_matrix_raw), center = TRUE, scale = TRUE))
ts_matrix_z[is.na(ts_matrix_z)] <- 0

z_df <- df
z_df[, log2fc_cols] <- ts_matrix_z

write.csv(
    z_df,
    file.path(outdir, "ZScored_log2FC_Matrix.csv"),
    row.names = FALSE
)

dist_mat <- dist(ts_matrix_z, method = distance)
hc <- hclust(dist_mat, method = method)
clusters <- cutree(hc, k = k)

clustered <- df
clustered$Cluster <- as.factor(clusters)

write.csv(
    clustered,
    file.path(outdir, "Clustered_Genes_Hierarchical_Euclidean_Z.csv"),
    row.names = FALSE
)

clustered_z <- z_df
clustered_z$Cluster <- as.factor(clusters)

write.csv(
    clustered_z,
    file.path(outdir, "Clustered_Genes_Hierarchical_Euclidean_Z_zscores.csv"),
    row.names = FALSE
)

pdf(file.path(outdir, "Dendrogram_Euclidean_Z.pdf"), width = 8, height = 6)
plot(
    hc,
    main = paste0("RNA-seq trajectory clustering, k = ", k),
    xlab = "",
    sub = "",
    labels = FALSE
)
rect.hclust(hc, k = k, border = "blue")
dev.off()

writeLines(
    c(
        paste0("input_file=", normalizePath(input_file)),
        paste0("distance=", distance),
        paste0("method=", method),
        paste0("k=", k),
        paste0("log2fc_columns=", paste(log2fc_cols, collapse = ",")),
        "matrix_transform=per-gene z-score across timepoints"
    ),
    file.path(outdir, "clustering_parameters.txt")
)

message("Wrote clustering results to: ", normalizePath(outdir))
