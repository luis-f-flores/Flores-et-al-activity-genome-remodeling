#!/usr/bin/env Rscript
# KCl RNA clustering: row-standardized normalized-count trajectories.
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



input_file <- gb_path("rna/K_clustering_Neuron_RNA_seq_11-13-24/filtered_results/combined_log2FC_FPKM_GlobalDEG.csv")
output_dir <- gb_path("rna/K_clustering_Neuron_RNA_seq_11-13-24/filtered_results/")
dir.create(output_dir, showWarnings = FALSE)

log2fc_cols <- c("log2FC_1hr", "log2FC_6hr", "log2FC_24hr")
df <- read.csv(input_file)


ts_matrix_raw <- as.matrix(df[, log2fc_cols])
ts_matrix_z <- t(scale(t(ts_matrix_raw), center = TRUE, scale = TRUE))


z_df <- df
z_df[, log2fc_cols] <- ts_matrix_z
write.csv(z_df, file.path(output_dir, "ZScored_log2FC_Matrix.csv"), row.names = FALSE)
cat("Z-scored matrix saved\n")


run_clustering <- function(matrix, label) {
  cat(paste0("[", label, "] Euclidean hierarchical clustering...\n"))
  dist_mat <- dist(matrix, method = "euclidean")
  hc <- hclust(dist_mat, method = "ward.D2")
  df_hc <- df
  df_hc$Cluster <- as.factor(cutree(hc, k = 5))
  write.csv(df_hc, file.path(output_dir, paste0("Clustered_Genes_Hierarchical_Euclidean_", label, ".csv")), row.names = FALSE)

  pdf(file.path(output_dir, paste0("Dendrogram_Euclidean_", label, ".pdf")))
  plot(hc, main = paste("Euclidean Hierarchical Dendrogram —", label), xlab = "", sub = "")
  rect.hclust(hc, k = 5, border = "blue")
  dev.off()


  cat(paste0("Clustering complete (Ward.D2) for: ", label, "\n"))
}


run_clustering(ts_matrix_z, "Z")
