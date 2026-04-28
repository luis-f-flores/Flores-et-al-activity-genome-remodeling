#!/usr/bin/env Rscript
# 07_GO_for_RNAseq_clusters.R
#
# Perform GO enrichment per RNA-seq trajectory cluster.

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    library(dplyr)
    library(writexl)
    library(enrichplot)
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

clustered_file <- get_arg("--clustered-file")
outdir <- get_arg("--outdir", "enrichment_results")
gene_id_col <- get_arg("--gene-id-col", "gene_id")
cluster_col <- get_arg("--cluster-col", "Cluster")
keytype <- toupper(get_arg("--keytype", "ENSEMBL"))
ont <- toupper(get_arg("--ont", "BP"))
p_cutoff <- as.numeric(get_arg("--p-cutoff", "0.05"))
q_cutoff <- as.numeric(get_arg("--q-cutoff", "0.05"))
show_category <- as.integer(get_arg("--show-category", "15"))

if (is.null(clustered_file) || !file.exists(clustered_file)) {
    stop("Usage: Rscript 07_GO_for_RNAseq_clusters.R --clustered-file <Clustered_Genes_Hierarchical_Euclidean_Z.csv>")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(clustered_file, check.names = FALSE)

if (!gene_id_col %in% colnames(df)) {
    stop("Missing gene ID column: ", gene_id_col)
}

if (!cluster_col %in% colnames(df)) {
    stop("Missing cluster column: ", cluster_col)
}

df$EntrezID <- AnnotationDbi::mapIds(
    org.Mm.eg.db,
    keys = as.character(df[[gene_id_col]]),
    column = "ENTREZID",
    keytype = keytype,
    multiVals = "first"
)

df <- df[!is.na(df$EntrezID), ]

clusters <- sort(unique(df[[cluster_col]]))

for (cluster_id in clusters) {
    message("Processing cluster: ", cluster_id)

    gene_ids <- df %>%
        filter(.data[[cluster_col]] == cluster_id) %>%
        pull(EntrezID) %>%
        unique()

    if (length(gene_ids) == 0) {
        next
    }

    ego <- enrichGO(
        gene = gene_ids,
        OrgDb = org.Mm.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        readable = TRUE
    )

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 1) {
        ego <- simplify(
            ego,
            cutoff = 0.7,
            by = "p.adjust",
            select_fun = min
        )
    }

    cluster_outdir <- file.path(outdir, paste0("Cluster_", cluster_id))
    dir.create(cluster_outdir, showWarnings = FALSE, recursive = TRUE)

    ego_df <- as.data.frame(ego)

    write_xlsx(
        list(GO = ego_df %>% arrange(p.adjust) %>% head(show_category)),
        file.path(cluster_outdir, paste0("GO_", ont, "_Top", show_category, ".xlsx"))
    )

    write.table(
        ego_df,
        file = file.path(cluster_outdir, paste0("GO_", ont, "_all_terms.tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    if (nrow(ego_df) > 0) {
        plot_obj <- dotplot(
            ego,
            showCategory = show_category,
            title = paste0("Cluster ", cluster_id, " GO ", ont)
        ) +
            theme(
                axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
                axis.text.y = element_text(size = 8)
            )

        ggsave(
            file.path(cluster_outdir, paste0("GO_", ont, "_dotplot.pdf")),
            plot_obj,
            width = 7,
            height = 5
        )

        ggsave(
            file.path(cluster_outdir, paste0("GO_", ont, "_dotplot.png")),
            plot_obj,
            width = 7,
            height = 5,
            dpi = 300
        )
    }
}

writeLines(
    c(
        paste0("clustered_file=", normalizePath(clustered_file)),
        paste0("gene_id_col=", gene_id_col),
        paste0("cluster_col=", cluster_col),
        paste0("keytype=", keytype),
        paste0("ontology=", ont),
        paste0("p_cutoff=", p_cutoff),
        paste0("q_cutoff=", q_cutoff),
        "GO_simplify_cutoff=0.7"
    ),
    file.path(outdir, "GO_parameters.txt")
)

message("Done. Output: ", normalizePath(outdir))
