#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = "")
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Run with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
if (!nzchar(GB_WORK_ROOT)) {
  GB_WORK_ROOT <- normalizePath(file.path(dirname(script_path), "..", "..", "..", "work"), mustWork = FALSE)
}

source_run <- Sys.getenv(
  "TELESCOPE_DE_SOURCE_RUN",
  gb_path("rna/telescope_de")
)
output_run <- Sys.getenv(
  "TELESCOPE_SUBFAMILY_OUTPUT_RUN",
  gb_path("rna/telescope_target_expression")
)
cli<-commandArgs(trailingOnly=TRUE)
if(length(cli)!=5L)stop("Usage: telescope_target_expression.R telescope_DE_directory TE_meta.tsv salmon_quant.sf salmon_DE.csv output_directory")
source_run<-cli[1];output_run<-cli[5]
table_dir <- file.path(output_run, "tables")
qc_dir <- file.path(output_run, "qc")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

te_meta_path <- Sys.getenv(
  "TELESCOPE_TE_META",
  gb_path("rna/TE_catalog_mm10/TE_meta.tsv")
)
salmon_quant_path <- Sys.getenv(
  "TELESCOPE_SALMON_QUANT",
  gb_path("rna/SalmonTE_out_6h/A6M/quant.sf")
)
salmon_de_path <- Sys.getenv(
  "TELESCOPE_SALMON_DE",
  gb_path("rna/SalmonTE_DE_6h/TE_DE_results.csv")
)

te_meta_path<-cli[2];salmon_quant_path<-cli[3];salmon_de_path<-cli[4]

sha256 <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE)
  if (length(line) != 1L) stop("Could not hash: ", path)
  sub(" .*", "", line)
}

read_count_matrix <- function(path) {
  x <- read.delim(gzfile(path), stringsAsFactors = FALSE, check.names = FALSE)
  if (names(x)[1] != "locus_id" || anyDuplicated(x$locus_id)) stop("Unexpected count matrix")
  ids <- x$locus_id
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "numeric"
  if (any(!is.finite(m)) || any(m < 0) || any(abs(m - round(m)) > 1e-8)) {
    stop("Counts are not finite non-negative integers")
  }
  rownames(m) <- ids
  storage.mode(m) <- "integer"
  m
}

counts <- read_count_matrix(file.path(source_run, "counts_final_count.tsv.gz"))
conf_counts <- read_count_matrix(file.path(source_run, "counts_final_conf.tsv.gz"))
if (!identical(dimnames(counts), dimnames(conf_counts))) stop("Count matrices are not aligned")

meta <- read.delim(te_meta_path, stringsAsFactors = FALSE, check.names = FALSE)
required_meta <- c("copyID", "family", "repClass", "repFamily")
if (!all(required_meta %in% names(meta)) || anyDuplicated(meta$copyID)) stop("Unexpected TE_meta.tsv")
meta <- meta[match(rownames(counts), meta$copyID), , drop = FALSE]
if (anyNA(meta$copyID) || !identical(meta$copyID, rownames(counts))) stop("Unmapped Telescope loci")

subfamilies <- unique(meta$family)
one_value <- function(values, label) {
  value <- unique(values)
  if (length(value) != 1L) stop("Multiple ", label, " annotations within subfamily")
  value
}
joined_values <- function(values) paste(sort(unique(values)), collapse = "|")
subfamily_meta <- data.frame(
  repeat_subfamily = subfamilies,
  repeat_class = vapply(subfamilies, function(x) one_value(meta$repClass[meta$family == x], "class"), character(1)),
  repeat_family = vapply(subfamilies, function(x) joined_values(meta$repFamily[meta$family == x]), character(1)),
  repeat_family_ambiguous = vapply(subfamilies, function(x) length(unique(meta$repFamily[meta$family == x])) > 1L, logical(1)),
  represented_loci = as.integer(table(factor(meta$family, levels = subfamilies))),
  stringsAsFactors = FALSE
)

salmon <- read.delim(salmon_quant_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!"Name" %in% names(salmon) || anyDuplicated(tolower(salmon$Name))) stop("Unexpected SalmonTE targets")
salmon_match <- match(tolower(subfamily_meta$repeat_subfamily), tolower(salmon$Name))
subfamily_meta$salmonte_target <- salmon$Name[salmon_match]
subfamily_meta$salmonte_case_normalized_match <- !is.na(salmon_match)
if (!all(subfamily_meta$salmonte_case_normalized_match)) stop("A Telescope subfamily lacks a SalmonTE target")

aggregate_by_subfamily <- function(matrix) {
  out <- rowsum(matrix, group = meta$family, reorder = FALSE)
  out[subfamilies, , drop = FALSE]
}
sub_counts <- aggregate_by_subfamily(counts)
sub_conf <- aggregate_by_subfamily(conf_counts)
if (!identical(dimnames(sub_counts), dimnames(sub_conf))) stop("Collapsed matrices are not aligned")

samples <- read.delim(file.path(source_run, "samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
samples <- samples[match(colnames(sub_counts), samples$sample), , drop = FALSE]
if (anyNA(samples$sample) || !identical(samples$sample, colnames(sub_counts))) stop("Sample table mismatch")
samples$condition <- factor(samples$condition, levels = c("Mannitol", "KCl"))
rownames(samples) <- samples$sample

normalization <- read.delim(file.path(source_run, "normalization_inputs.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
normalization <- normalization[match(colnames(sub_counts), normalization$sample), , drop = FALSE]
size_factors <- normalization$mapped_fragment_size_factor
names(size_factors) <- normalization$sample
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) stop("Invalid size factors")

eligible <- rowSums(sub_counts >= 10L) >= 3L
sub_counts_fit <- sub_counts[eligible, , drop = FALSE]
sub_conf_fit <- sub_conf[eligible, , drop = FALSE]
if (nrow(sub_counts_fit) < 2L) stop("Too few subfamilies after filtering")

fit_model <- function(matrix) {
  dds <- DESeqDataSetFromMatrix(
    countData = matrix,
    colData = samples,
    design = ~ condition
  )
  sizeFactors(dds) <- size_factors[colnames(matrix)]
  dds <- DESeq(dds, quiet = TRUE)
  coefficient <- "condition_KCl_vs_Mannitol"
  if (!coefficient %in% resultsNames(dds)) stop("Expected coefficient missing")
  raw <- results(dds, contrast = c("condition", "KCl", "Mannitol"))
  shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm", quiet = TRUE)
  table <- data.frame(
    repeat_subfamily = rownames(matrix),
    baseMean = raw$baseMean,
    log2FoldChange_MLE = raw$log2FoldChange,
    lfcSE_MLE = raw$lfcSE,
    stat = raw$stat,
    pvalue = raw$pvalue,
    padj_BH_all_tested = p.adjust(raw$pvalue, method = "BH"),
    padj_DESeq2_independent_filter = raw$padj,
    log2FoldChange_apeglm = shrunk$log2FoldChange,
    lfcSE_apeglm = shrunk$lfcSE,
    stringsAsFactors = FALSE
  )
  list(dds = dds, table = table)
}

primary <- fit_model(sub_counts_fit)
confidence <- fit_model(sub_conf_fit)
confidence_table <- confidence$table[, c("repeat_subfamily", "baseMean", "log2FoldChange_apeglm", "lfcSE_apeglm")]
names(confidence_table)[-1] <- paste0(names(confidence_table)[-1], "_final_conf")

result <- merge(primary$table, confidence_table, by = "repeat_subfamily", all.x = TRUE, sort = FALSE)
result <- merge(result, subfamily_meta, by = "repeat_subfamily", all.x = TRUE, sort = FALSE)
result <- result[match(rownames(sub_counts_fit), result$repeat_subfamily), , drop = FALSE]
result$pooled_confident_fraction <- rowSums(sub_conf_fit) / pmax(rowSums(sub_counts_fit), 1)
result$samples_final_count_ge_10 <- rowSums(sub_counts_fit >= 10L)
result$samples_final_conf_ge_10 <- rowSums(sub_conf_fit >= 10L)

salmon_de <- read.csv(salmon_de_path, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
salmon_index <- match(tolower(result$repeat_subfamily), tolower(rownames(salmon_de)))
result$salmonte_6h_log2FoldChange <- salmon_de$log2FoldChange[salmon_index]
result$salmonte_6h_padj <- salmon_de$padj[salmon_index]

write_gz <- function(data, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(data, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
count_table <- data.frame(repeat_subfamily = rownames(sub_counts), sub_counts, check.names = FALSE)
conf_table <- data.frame(repeat_subfamily = rownames(sub_conf), sub_conf, check.names = FALSE)
write_gz(count_table, file.path(table_dir, "subfamily_counts_final_count.tsv.gz"))
write_gz(conf_table, file.path(table_dir, "subfamily_counts_final_conf.tsv.gz"))
write.table(subfamily_meta, file.path(table_dir, "subfamily_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(result, file.path(table_dir, "subfamily_DE_KCl_vs_Mannitol.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

finite_assignment_comparison <- is.finite(result$log2FoldChange_apeglm) & is.finite(result$log2FoldChange_apeglm_final_conf)
finite_salmon <- is.finite(result$log2FoldChange_apeglm) & is.finite(result$salmonte_6h_log2FoldChange)
summary <- list(
  status = "complete",
  unit = "TE subfamily formed by summing Telescope-assigned locus counts with the same TE_meta.tsv family label",
  primary_count = "final_count",
  comparison_count = "final_conf",
  normalization = "total genome-mapped fragments / geometric mean; inherited from validated locus model",
  design = "DESeq2 ~ condition; KCl versus Mannitol; three independent libraries per condition",
  filter = "summed subfamily count >=10 in at least 3 of 6 libraries",
  effect_size = "apeglm-shrunken log2 fold change",
  p_values = "DESeq2 Wald test; BH across all tested subfamilies",
  catalog_loci = nrow(counts),
  catalog_subfamilies = nrow(sub_counts),
  tested_subfamilies = nrow(result),
  salmonte_targets = nrow(salmon),
  salmonte_case_normalized_matches = sum(subfamily_meta$salmonte_case_normalized_match),
  salmonte_targets_absent_from_telescope_catalog = nrow(salmon) - sum(subfamily_meta$salmonte_case_normalized_match),
  significant_up_BH_0_05 = sum(result$padj_BH_all_tested < 0.05 & result$log2FoldChange_apeglm > 0, na.rm = TRUE),
  significant_down_BH_0_05 = sum(result$padj_BH_all_tested < 0.05 & result$log2FoldChange_apeglm < 0, na.rm = TRUE),
  final_count_vs_final_conf_spearman = cor(result$log2FoldChange_apeglm[finite_assignment_comparison], result$log2FoldChange_apeglm_final_conf[finite_assignment_comparison], method = "spearman"),
  telescope_vs_salmonte_6h_spearman_descriptive = cor(result$log2FoldChange_apeglm[finite_salmon], result$salmonte_6h_log2FoldChange[finite_salmon], method = "spearman"),
  cross_method_comparison = "Telescope and SalmonTE retain their respective assignment/reference models. Missing SalmonTE-only targets are not zero-filled.",
  source_hashes = list(
    te_meta = sha256(te_meta_path),
    counts_final_count = sha256(file.path(source_run, "counts_final_count.tsv.gz")),
    counts_final_conf = sha256(file.path(source_run, "counts_final_conf.tsv.gz")),
    salmonte_quant_target_source = sha256(salmon_quant_path),
    salmonte_6h_de = sha256(salmon_de_path)
  )
)
write_json(summary, file.path(qc_dir, "SUMMARY.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))

cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")
