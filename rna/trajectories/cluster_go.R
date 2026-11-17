#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(writexl)


directory <- gb_path("rna_clusters")
clustered_file <- file.path(directory, "Clustered_Genes_Hierarchical_Euclidean_Z.csv")


perform_enrichment_per_cluster <- function(clustered_file, output_dir = "enrichment_results_Euclidean") {
  dir.create(output_dir, showWarnings = FALSE)

  df <- read.csv(clustered_file)

  if (!("Cluster" %in% colnames(df)) || !("gene_id" %in% colnames(df))) {
    stop("File must contain 'Cluster' and 'gene_id' columns.")
  }


  df$EntrezID <- mapIds(org.Mm.eg.db,
                        keys = df$gene_id,
                        column = "ENTREZID",
                        keytype = "ENSEMBL",
                        multiVals = "first")

  df <- df[!is.na(df$EntrezID), ]

  unique_clusters <- unique(df$Cluster)

  for (cluster in unique_clusters) {
    cat("Processing Cluster:", cluster, "\n")

    gene_ids <- df %>% filter(Cluster == cluster) %>% pull(EntrezID)
    if (length(gene_ids) == 0) next

    cluster_output_dir <- file.path(output_dir, paste0("Cluster_", cluster))
    dir.create(cluster_output_dir, showWarnings = FALSE)


    ego_bp <- enrichGO(gene          = gene_ids,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)


     if (!is.null(ego_bp) && nrow(ego_bp@result) > 1) {
       ego_bp <- simplify(ego_bp,
                          cutoff     = 0.7,
                          by         = "p.adjust",
                          select_fun = min)
     }


    ego_cc <- enrichGO(gene          = gene_ids,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)


     if (!is.null(ego_cc) && nrow(ego_cc@result) > 1) {
       ego_cc <- simplify(ego_cc,
                          cutoff     = 0.7,
                          by         = "p.adjust",
                          select_fun = min)
     }


    save_top_terms <- function(enrichment_result, filename) {
      if (!is.null(enrichment_result) && nrow(enrichment_result@result) > 0) {
        enrichment_result@result$p.adjust <- as.numeric(as.character(enrichment_result@result$p.adjust))
        write_xlsx(enrichment_result@result %>% arrange(p.adjust) %>% head(15),
                   file.path(cluster_output_dir, filename))
        cat("Saved top 15 terms to", filename, "in", cluster_output_dir, "\n")
      } else {
        cat("No significant terms found for", filename, "\n")
      }
    }

    save_top_terms(ego_bp, "GO_Biological_Process_Top15.xlsx")
    save_top_terms(ego_cc, "GO_Cellular_Component_Top15.xlsx")
  }

  cat("Enrichment analysis for all clusters completed.\n")
}


perform_enrichment_per_cluster(clustered_file)
