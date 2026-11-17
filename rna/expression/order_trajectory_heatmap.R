#!/usr/bin/env Rscript
# Order genes within the five normalized-count trajectory clusters.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript order_trajectory_heatmap.R cluster_assignments.csv output.csv")
clusters <- read.csv(args[1], check.names = FALSE)
stopifnot(!anyDuplicated(clusters$gene_id))
values <- as.matrix(clusters[, c("log2FC_1hr", "log2FC_6hr", "log2FC_24hr")])
standardized <- t(scale(t(values)))
stopifnot(all(is.finite(standardized)))
hierarchy <- hclust(dist(standardized, method = "euclidean"), method = "ward.D2")
stopifnot(all(as.integer(clusters$Cluster) == cutree(hierarchy, k = 5)))
row_order <- order(as.integer(clusters$Cluster), match(seq_len(nrow(clusters)), hierarchy$order))
genes <- data.frame(row_order = seq_along(row_order), source_row = row_order,
                    clusters[row_order, c("gene_id", "gene_name", "Cluster")], check.names = FALSE)
rownames(genes) <- NULL
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.csv(genes, args[2], row.names = FALSE)
