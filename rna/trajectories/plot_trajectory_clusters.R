#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}



library(dplyr)
library(ggplot2)
library(factoextra)
library(pheatmap)
library(tidyr)


output_dir <- gb_path("rna/K_clustering_Neuron_RNA_seq_11-13-24/filtered_results/")
input_files <- list(
  Hierarchical_Euclidean_Z = "Clustered_Genes_Hierarchical_Euclidean_Z.csv"
)


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


for (method_name in names(input_files)) {
  input_file <- paste0(output_dir, input_files[[method_name]])

  if (!file.exists(input_file)) {
    cat("Skipping", method_name, "— file not found:", input_file, "\n")
    next
  }

  cat("Processing:", method_name, "\n")

  df <- read.csv(input_file)
  df$Cluster <- as.factor(df$Cluster)


  ts_data <- split(df[, c("log2FC_1hr", "log2FC_6hr", "log2FC_24hr")], seq(nrow(df)))
  ts_matrix <- do.call(rbind, ts_data)

  pdf(paste0(output_dir, method_name, "_Clustering_PCA.pdf"), width = 8, height = 6)
  fviz_cluster(list(data = ts_matrix, cluster = as.numeric(df$Cluster)),
               geom = "point",
               ellipse.type = "convex",
               ellipse.level = 0.95,
               ggtheme = theme_minimal(base_size = 14) +
                        theme(panel.grid.minor = element_blank())) +
    labs(title = paste(method_name, "Clustering (PCA)"))
  dev.off()
  cat("PCA plot saved for", method_name, "\n")


cat("Generating line plot...\n")
long_data <- df %>%
  select(starts_with("log2FC"), Cluster) %>%
  pivot_longer(cols = starts_with("log2FC"), names_to = "Timepoint", values_to = "log2FC") %>%
  mutate(Timepoint = factor(gsub("log2FC_", "", Timepoint), levels = c("1hr", "6hr", "24hr")),
         Cluster = as.factor(Cluster))


cat("Unique clusters found:", unique(long_data$Cluster), "\n")

trends <- long_data %>%
  group_by(Cluster, Timepoint) %>%
  summarize(mean_log2FC = mean(log2FC, na.rm = TRUE),
            se_log2FC = sd(log2FC, na.rm = TRUE) / sqrt(n()),
            .groups = "drop")

cat("Trend matrix dims:", dim(trends), "\n")


if (nrow(trends) > 0) {
  pdf(paste0(output_dir, method_name, "_Cluster_Trends.pdf"), width = 8, height = 6)
  print(
    ggplot(trends, aes(x = Timepoint, y = mean_log2FC, group = Cluster, color = Cluster)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = mean_log2FC - se_log2FC, ymax = mean_log2FC + se_log2FC),
                    width = 0.2, linewidth = 0.8) +
      facet_wrap(~ Cluster, ncol = 2) +
      labs(title = paste("Temporal Expression —", method_name),
           x = "Timepoint", y = "Mean log2FC") +
      theme_minimal(base_size = 14) +
      theme(panel.grid.minor = element_blank(),
            legend.position = "none")
  )
  dev.off()
  cat("Line plot saved for", method_name, "\n")
} else {
  cat("No trend data — skipping plot for", method_name, "\n")
}

  summary_stats <- long_data %>%
    group_by(Cluster, Timepoint) %>%
    summarize(mean_log2FC = mean(log2FC),
              se_log2FC = sd(log2FC) / sqrt(n()),
              n_genes = n(),
              .groups = "drop")

  total_genes <- df %>%
    group_by(Cluster) %>%
    summarize(total_genes = n(), .groups = "drop")

  stats_export <- left_join(summary_stats, total_genes, by = "Cluster")

  write.csv(stats_export,
            file = paste0(output_dir, method_name, "_Cluster_Stats.csv"),
            row.names = FALSE)
  cat("Stats CSV exported for", method_name, "\n")


  df_heat <- df %>%
    select(Cluster, log2FC_1hr, log2FC_6hr, log2FC_24hr)
  heatmap_matrix <- as.matrix(df_heat[, -1])

  heatmap_colors <- colorRampPalette(c("blue", "#f0f0f0", "red"))(250)

  pdf(file = paste0(output_dir, method_name, "_Heatmap.pdf"), width = 5, height = 9)
  print(
    pheatmap(heatmap_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = heatmap_colors,
             main = paste("Heatmap —", method_name),
             show_rownames = FALSE,
             show_colnames = TRUE,
             labels_col = c("1hr", "6hr", "24hr"),
             breaks = seq(-2, 2, length.out = 250))
  )
  dev.off()

  png(file = paste0(output_dir, method_name, "_Heatmap.png"), width = 1000, height = 1800, res = 200)
  print(
    pheatmap(heatmap_matrix,
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             color = heatmap_colors,
             main = paste("Heatmap —", method_name),
             show_rownames = FALSE,
             show_colnames = TRUE,
             labels_col = c("1hr", "6hr", "24hr"),
             breaks = seq(-2, 2, length.out = 250))
  )
  dev.off()

  write.csv(df_heat,
            file = paste0(output_dir, method_name, "_Heatmap_Matrix.csv"),
            row.names = FALSE)

  cat("Heatmap and matrix exported for", method_name, "\n\n")
}
