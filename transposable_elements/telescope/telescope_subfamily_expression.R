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
  "TELESCOPE_RMSK_SUBFAMILY_OUTPUT_RUN",
  gb_path("rna/telescope_repeatmasker_subfamily")
)
master_path <- Sys.getenv(
  "TE_ANNOTATION_MASTER",
  gb_path("rna/TE_catalog_mm10/TE_annotation_master_mm10.tsv.gz")
)

cli<-commandArgs(trailingOnly=TRUE)
if(length(cli)!=3L)stop("Usage: telescope_subfamily_expression.R telescope_DE_directory annotation_master output_directory")
source_run<-cli[1];master_path<-cli[2];output_run<-cli[3]
table_dir <- file.path(output_run, "tables")
qc_dir <- file.path(output_run, "qc")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

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

counts_path <- file.path(source_run, "counts_final_count.tsv.gz")
conf_counts_path <- file.path(source_run, "counts_final_conf.tsv.gz")
samples_path <- file.path(source_run, "samples.tsv")
normalization_path <- file.path(source_run, "normalization_inputs.tsv")

counts <- read_count_matrix(counts_path)
conf_counts <- read_count_matrix(conf_counts_path)
if (!identical(dimnames(counts), dimnames(conf_counts))) stop("Count matrices are not aligned")

master <- read.delim(gzfile(master_path), stringsAsFactors = FALSE, check.names = FALSE)
required_master <- c(
  "locus_id", "raw_repeat_name", "raw_repeat_class", "raw_repeat_family",
  "salmonte_target", "is_iap", "iap_component"
)
if (!all(required_master %in% names(master)) || anyDuplicated(master$locus_id)) {
  stop("Unexpected TE annotation master")
}
meta <- master[match(rownames(counts), master$locus_id), required_master, drop = FALSE]
rm(master)
gc(verbose = FALSE)
if (anyNA(meta$locus_id) || !identical(meta$locus_id, rownames(counts))) {
  stop("Telescope loci do not map exactly to the TE annotation master")
}

raw_names <- unique(meta$raw_repeat_name)
one_value <- function(values, label) {
  value <- unique(values)
  if (length(value) != 1L) stop("Multiple ", label, " annotations within raw RepeatMasker name")
  value
}
joined_values <- function(values) paste(sort(unique(values)), collapse = "|")
target_raw_name_count <- table(meta$salmonte_target, meta$raw_repeat_name)
target_raw_name_count <- rowSums(target_raw_name_count > 0)

raw_meta <- data.frame(
  repeatmasker_subfamily = raw_names,
  raw_repeat_class = vapply(raw_names, function(x) one_value(meta$raw_repeat_class[meta$raw_repeat_name == x], "class"), character(1)),
  raw_repeat_family = vapply(raw_names, function(x) one_value(meta$raw_repeat_family[meta$raw_repeat_name == x], "family"), character(1)),
  salmonte_target = vapply(raw_names, function(x) one_value(meta$salmonte_target[meta$raw_repeat_name == x], "SalmonTE target"), character(1)),
  represented_loci = as.integer(table(factor(meta$raw_repeat_name, levels = raw_names))),
  is_iap = vapply(raw_names, function(x) any(meta$is_iap[meta$raw_repeat_name == x]), logical(1)),
  iap_component = vapply(raw_names, function(x) joined_values(meta$iap_component[meta$raw_repeat_name == x]), character(1)),
  stringsAsFactors = FALSE
)
raw_meta$raw_names_per_salmonte_target <- as.integer(target_raw_name_count[raw_meta$salmonte_target])

aggregate_by_raw_name <- function(matrix) {
  out <- rowsum(matrix, group = meta$raw_repeat_name, reorder = FALSE)
  out[raw_names, , drop = FALSE]
}
raw_counts <- aggregate_by_raw_name(counts)
raw_conf <- aggregate_by_raw_name(conf_counts)
if (!identical(dimnames(raw_counts), dimnames(raw_conf))) stop("Aggregated matrices are not aligned")

samples <- read.delim(samples_path, stringsAsFactors = FALSE, check.names = FALSE)
samples <- samples[match(colnames(raw_counts), samples$sample), , drop = FALSE]
if (anyNA(samples$sample) || !identical(samples$sample, colnames(raw_counts))) stop("Sample table mismatch")
samples$condition <- factor(samples$condition, levels = c("Mannitol", "KCl"))
rownames(samples) <- samples$sample

normalization <- read.delim(normalization_path, stringsAsFactors = FALSE, check.names = FALSE)
normalization <- normalization[match(colnames(raw_counts), normalization$sample), , drop = FALSE]
size_factors <- normalization$mapped_fragment_size_factor
names(size_factors) <- normalization$sample
if (any(!is.finite(size_factors)) || any(size_factors <= 0)) stop("Invalid size factors")

eligible <- rowSums(raw_counts >= 10L) >= 3L
raw_counts_fit <- raw_counts[eligible, , drop = FALSE]
raw_conf_fit <- raw_conf[eligible, , drop = FALSE]
if (nrow(raw_counts_fit) < 2L) stop("Too few RepeatMasker subfamilies after filtering")

fit_model <- function(matrix) {
  dds <- DESeqDataSetFromMatrix(countData = matrix, colData = samples, design = ~ condition)
  sizeFactors(dds) <- size_factors[colnames(matrix)]
  dds <- DESeq(dds, quiet = TRUE)
  coefficient <- "condition_KCl_vs_Mannitol"
  if (!coefficient %in% resultsNames(dds)) stop("Expected coefficient missing")
  raw <- results(dds, contrast = c("condition", "KCl", "Mannitol"))
  shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm", quiet = TRUE)
  table <- data.frame(
    repeatmasker_subfamily = rownames(matrix),
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

primary <- fit_model(raw_counts_fit)
confidence <- fit_model(raw_conf_fit)
confidence_table <- confidence$table[, c(
  "repeatmasker_subfamily", "baseMean", "log2FoldChange_apeglm", "lfcSE_apeglm"
)]
names(confidence_table)[-1] <- paste0(names(confidence_table)[-1], "_final_conf")

result <- merge(primary$table, confidence_table, by = "repeatmasker_subfamily", all.x = TRUE, sort = FALSE)
result <- merge(result, raw_meta, by = "repeatmasker_subfamily", all.x = TRUE, sort = FALSE)
result <- result[match(rownames(raw_counts_fit), result$repeatmasker_subfamily), , drop = FALSE]
result$pooled_confident_fraction <- rowSums(raw_conf_fit) / pmax(rowSums(raw_counts_fit), 1)
result$samples_final_count_ge_10 <- rowSums(raw_counts_fit >= 10L)
result$samples_final_conf_ge_10 <- rowSums(raw_conf_fit >= 10L)

write_gz <- function(data, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(data, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
write_gz(data.frame(repeatmasker_subfamily = rownames(raw_counts), raw_counts, check.names = FALSE),
         file.path(table_dir, "repeatmasker_subfamily_counts_final_count.tsv.gz"))
write_gz(data.frame(repeatmasker_subfamily = rownames(raw_conf), raw_conf, check.names = FALSE),
         file.path(table_dir, "repeatmasker_subfamily_counts_final_conf.tsv.gz"))
write.table(raw_meta, file.path(table_dir, "repeatmasker_subfamily_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(result, file.path(table_dir, "repeatmasker_subfamily_DE_KCl_vs_Mannitol.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

finite_assignment_comparison <- is.finite(result$log2FoldChange_apeglm) &
  is.finite(result$log2FoldChange_apeglm_final_conf)
summary <- list(
  status = "complete",
  unit = "original mm10 RepeatMasker raw_repeat_name consensus/subfamily after summing represented Telescope loci",
  upstream_assignment = "Telescope final_count with final_conf high-confidence assignments",
  primary_count = "final_count",
  comparison_count = "final_conf",
  normalization = "total genome-mapped fragments / geometric mean; inherited from validated locus model",
  design = "DESeq2 ~ condition; KCl versus Mannitol; three independent libraries per condition",
  filter = "summed raw RepeatMasker subfamily count >=10 in at least 3 of 6 libraries",
  effect_size = "apeglm-shrunken log2 fold change",
  p_values = "DESeq2 Wald test; BH across all tested raw RepeatMasker subfamilies",
  telescope_catalog_loci = nrow(counts),
  catalog_repeatmasker_subfamilies = nrow(raw_counts),
  tested_repeatmasker_subfamilies = nrow(result),
  significant_up_BH_0_05 = sum(result$padj_BH_all_tested < 0.05 & result$log2FoldChange_apeglm > 0, na.rm = TRUE),
  significant_down_BH_0_05 = sum(result$padj_BH_all_tested < 0.05 & result$log2FoldChange_apeglm < 0, na.rm = TRUE),
  final_count_vs_final_conf_spearman = cor(
    result$log2FoldChange_apeglm[finite_assignment_comparison],
    result$log2FoldChange_apeglm_final_conf[finite_assignment_comparison], method = "spearman"
  ),
  raw_names_collapsed_by_salmonte_aliases = sum(raw_meta$raw_names_per_salmonte_target > 1L),
  analysis_boundary = "downstream regrouping and DE refit only; no Telescope alignment or EM assignment was rerun",
  source_hashes = list(
    te_annotation_master = sha256(master_path),
    counts_final_count = sha256(counts_path),
    counts_final_conf = sha256(conf_counts_path),
    samples = sha256(samples_path),
    normalization_inputs = sha256(normalization_path)
  )
)
write_json(summary, file.path(qc_dir, "SUMMARY.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))
cat(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), "\n")
