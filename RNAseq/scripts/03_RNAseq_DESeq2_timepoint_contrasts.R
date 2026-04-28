#!/usr/bin/env Rscript
# 03_RNAseq_DESeq2_timepoint_contrasts.R
#
# Per-timepoint DESeq2 contrasts for KCl versus mannitol.

suppressPackageStartupMessages({
    library(DESeq2)
    library(readr)
    library(dplyr)
    library(tibble)
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

counts_file <- get_arg("--counts")
metadata_file <- get_arg("--metadata")
contrast_file <- get_arg("--contrast-sheet")
outdir <- get_arg("--outdir", "DESeq2_results")
fc_threshold <- as.numeric(get_arg("--fc", "1.3"))
padj_threshold <- as.numeric(get_arg("--padj", "0.01"))

if (is.null(counts_file) || is.null(metadata_file) || is.null(contrast_file)) {
    stop("Required arguments: --counts, --metadata, --contrast-sheet")
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_table_auto <- function(path) {
    if (grepl("\\.csv$", path, ignore.case = TRUE)) {
        read_csv(path, show_col_types = FALSE)
    } else {
        read_tsv(path, show_col_types = FALSE)
    }
}

counts_df <- read_table_auto(counts_file)
metadata <- read_csv(metadata_file, show_col_types = FALSE)
contrasts <- read_csv(contrast_file, show_col_types = FALSE)

required_metadata <- c("sample_id", "timepoint", "treatment")
missing_metadata <- setdiff(required_metadata, names(metadata))

if (length(missing_metadata) > 0) {
    stop("Metadata missing required columns: ", paste(missing_metadata, collapse = ", "))
}

gene_col <- names(counts_df)[1]

counts <- as.data.frame(counts_df)
rownames(counts) <- counts[[gene_col]]
counts[[gene_col]] <- NULL
counts <- round(as.matrix(counts))

sample_keep <- intersect(metadata$sample_id, colnames(counts))

metadata <- metadata %>%
    filter(sample_id %in% sample_keep)

counts <- counts[, metadata$sample_id, drop = FALSE]

for (i in seq_len(nrow(contrasts))) {
    timepoint <- as.character(contrasts$timepoint[i])
    control <- as.character(contrasts$control_treatment[i])
    test <- as.character(contrasts$test_treatment[i])
    prefix <- as.character(contrasts$output_prefix[i])

    meta_sub <- metadata %>%
        filter(timepoint == !!timepoint, treatment %in% c(control, test)) %>%
        mutate(treatment = factor(treatment, levels = c(control, test)))

    if (nrow(meta_sub) < 4) {
        warning("Skipping ", prefix, ": too few samples")
        next
    }

    count_sub <- counts[, meta_sub$sample_id, drop = FALSE]

    keep <- rowSums(count_sub >= 10) >= 2
    count_sub <- count_sub[keep, , drop = FALSE]

    dds <- DESeqDataSetFromMatrix(
        countData = count_sub,
        colData = as.data.frame(meta_sub),
        design = ~ treatment
    )

    dds <- DESeq(dds)

    res <- results(dds, contrast = c("treatment", test, control))

    res_df <- as.data.frame(res) %>%
        rownames_to_column("gene_id") %>%
        arrange(padj)

    write_tsv(
        res_df,
        file.path(outdir, paste0(prefix, "_DESeq2_all_results.tsv"))
    )

    sig <- res_df %>%
        filter(
            !is.na(padj),
            padj <= padj_threshold,
            !is.na(log2FoldChange),
            abs(log2FoldChange) >= log2(fc_threshold)
        )

    write_tsv(
        sig,
        file.path(outdir, paste0(prefix, "_DESeq2_DEGs_padj", padj_threshold, "_FC", fc_threshold, ".tsv"))
    )

    message("[done] ", prefix, ": all = ", nrow(res_df), "; significant = ", nrow(sig))
}

writeLines(
    c(
        paste0("counts_file=", normalizePath(counts_file)),
        paste0("metadata_file=", normalizePath(metadata_file)),
        paste0("contrast_file=", normalizePath(contrast_file)),
        paste0("fc_threshold=", fc_threshold),
        paste0("log2fc_threshold=log2(", fc_threshold, ")=", log2(fc_threshold)),
        paste0("padj_threshold=", padj_threshold),
        "model=DESeq2 per-timepoint design ~ treatment"
    ),
    file.path(outdir, "DESeq2_parameters.txt")
)
