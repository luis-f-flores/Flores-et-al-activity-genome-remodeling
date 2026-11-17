#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


script_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the external directory prepared by prepare_telescope_de_inputs.py")
package_root <- normalizePath(args[[1]], mustWork = TRUE)

read_config <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(names(x), c("key", "value"))) stop("analysis_config.tsv must have key and value columns")
  if (anyDuplicated(x$key)) stop("Duplicate keys in analysis_config.tsv")
  setNames(x$value, x$key)
}

cfg <- read_config(file.path(package_root, "analysis_config.tsv"))
required_cfg <- c(
  "analysis_state", "design_mode", "contrast_numerator", "contrast_denominator",
  "prefilter_min_count", "prefilter_min_samples", "alpha",
  "primary_normalization", "primary_count", "lfc_shrinkage"
)
missing_cfg <- setdiff(required_cfg, names(cfg))
if (length(missing_cfg)) stop("Missing configuration keys: ", paste(missing_cfg, collapse = ", "))
if (cfg[["analysis_state"]] != "ready") {
  stop("analysis_state is '", cfg[["analysis_state"]], "'; expected 'ready' for the configured model")
}
if (!cfg[["design_mode"]] %in% c("paired", "unpaired")) {
  stop("design_mode must be 'paired' or 'unpaired'; current value: ", cfg[["design_mode"]])
}
if (cfg[["primary_count"]] != "final_count") stop("Primary count must remain final_count")
if (cfg[["primary_normalization"]] != "mapped_fragments") stop("Primary normalization must remain mapped_fragments")
if (cfg[["lfc_shrinkage"]] != "apeglm") stop("This model uses apeglm shrinkage")

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(jsonlite)
})

min_count <- as.integer(cfg[["prefilter_min_count"]])
min_samples <- as.integer(cfg[["prefilter_min_samples"]])
alpha <- as.numeric(cfg[["alpha"]])
if (!is.finite(min_count) || min_count < 1L) stop("Invalid prefilter_min_count")
if (!is.finite(min_samples) || min_samples < 1L) stop("Invalid prefilter_min_samples")
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("Invalid alpha")

read_count_matrix <- function(path) {
  x <- read.delim(gzfile(path), row.names = 1, check.names = FALSE)
  m <- as.matrix(x)
  storage.mode(m) <- "numeric"
  if (any(!is.finite(m)) || any(m < 0) || any(abs(m - round(m)) > 1e-8)) {
    stop("Count matrix contains invalid values: ", path)
  }
  storage.mode(m) <- "integer"
  m
}

counts <- read_count_matrix(file.path(package_root, "counts_final_count.tsv.gz"))
conf_counts <- read_count_matrix(file.path(package_root, "counts_final_conf.tsv.gz"))
if (!identical(dim(counts), dim(conf_counts)) ||
    !identical(rownames(counts), rownames(conf_counts)) ||
    !identical(colnames(counts), colnames(conf_counts))) {
  stop("final_count and final_conf matrices are not aligned")
}

samples <- read.delim(file.path(package_root, "samples.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
if (anyDuplicated(samples$sample)) stop("Duplicate sample names")
if (!setequal(samples$sample, colnames(counts))) stop("Sample metadata and count columns differ")
samples <- samples[match(colnames(counts), samples$sample), , drop = FALSE]
rownames(samples) <- samples$sample
samples$condition <- factor(samples$condition, levels = c(cfg[["contrast_denominator"]], cfg[["contrast_numerator"]]))
samples$replicate_label <- factor(samples$replicate_label)
if (anyNA(samples$condition)) stop("Configured contrast levels do not match sample conditions")

if (cfg[["design_mode"]] == "paired") {
  block_table <- table(samples$replicate_label, samples$condition)
  if (!all(block_table == 1L)) stop("Paired design requires one sample from each condition in every replicate block")
  design_formula <- ~ replicate_label + condition
} else {
  design_formula <- ~ condition
}
model_matrix <- model.matrix(design_formula, data = samples)
if (qr(model_matrix)$rank != ncol(model_matrix)) stop("Design matrix is not full rank")

norm_input <- read.delim(file.path(package_root, "normalization_inputs.tsv"), stringsAsFactors = FALSE)
norm_input <- norm_input[match(colnames(counts), norm_input$sample), , drop = FALSE]
if (anyNA(norm_input$sample)) stop("Missing normalization rows")
mapped_size_factors <- norm_input$mapped_fragment_size_factor
names(mapped_size_factors) <- norm_input$sample
if (any(!is.finite(mapped_size_factors)) || any(mapped_size_factors <= 0)) stop("Invalid mapped-fragment size factors")

keep <- rowSums(counts >= min_count) >= min_samples
if (sum(keep) < 2L) stop("Technical filter retained fewer than two loci")
counts_fit <- counts[keep, , drop = FALSE]
conf_fit <- conf_counts[keep, , drop = FALSE]

dir.create(file.path(package_root, "results"), showWarnings = FALSE)
dir.create(file.path(package_root, "figures"), showWarnings = FALSE)
dir.create(file.path(package_root, "qc"), showWarnings = FALSE)

fit_deseq <- function(matrix, size_factors, label) {
  dds <- DESeqDataSetFromMatrix(
    countData = matrix,
    colData = samples,
    design = design_formula
  )
  sizeFactors(dds) <- size_factors[colnames(dds)]
  dds <- DESeq(dds, minReplicatesForReplace = Inf, quiet = TRUE)
  coefficient_candidates <- grep("^condition_", resultsNames(dds), value = TRUE)
  if (length(coefficient_candidates) != 1L) {
    stop("Expected one condition coefficient for ", label, "; found: ", paste(coefficient_candidates, collapse = ", "))
  }
  coefficient <- coefficient_candidates[[1]]
  contrast <- c("condition", cfg[["contrast_numerator"]], cfg[["contrast_denominator"]])
  res_all <- results(
    dds, contrast = contrast, alpha = alpha,
    independentFiltering = FALSE, cooksCutoff = TRUE
  )
  res_if <- results(
    dds, contrast = contrast, alpha = alpha,
    independentFiltering = TRUE, cooksCutoff = TRUE
  )
  res_shrunk <- lfcShrink(dds, coef = coefficient, type = "apeglm", quiet = TRUE)
  table <- data.frame(
    locus_id = rownames(dds),
    baseMean = res_all$baseMean,
    log2FoldChange = res_all$log2FoldChange,
    lfcSE = res_all$lfcSE,
    stat = res_all$stat,
    pvalue = res_all$pvalue,
    padj_BH_all_retained = res_all$padj,
    padj_DESeq2_independent_filter = res_if$padj,
    log2FoldChange_shrunk = res_shrunk$log2FoldChange,
    lfcSE_shrunk = res_shrunk$lfcSE,
    cooks_flag = is.na(res_all$pvalue) & res_all$baseMean > 0,
    stringsAsFactors = FALSE
  )
  list(dds = dds, coefficient = coefficient, table = table, shrunk = res_shrunk)
}

primary <- fit_deseq(counts_fit, mapped_size_factors, "final_count/mapped")
conf_result <- fit_deseq(conf_fit, mapped_size_factors, "final_conf/mapped")

support <- read.delim(
  gzfile(file.path(package_root, "locus_support.tsv.gz")),
  row.names = 1, check.names = FALSE
)
support_columns <- c(
  "transcript_length", "sum_final_count", "sum_final_conf",
  "pooled_confident_fraction", "samples_positive_final_count",
  "samples_final_count_ge_10", "samples_final_conf_ge_10",
  "sum_unique_count", "sum_init_best"
)
if (!all(support_columns %in% names(support))) stop("locus_support.tsv.gz lacks required columns")

add_support <- function(result_table) {
  cbind(
    result_table,
    support[result_table$locus_id, support_columns, drop = FALSE],
    row.names = NULL
  )
}

write_gz_table <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

primary_table <- add_support(primary$table)
conf_table <- conf_result$table

write_gz_table(primary_table, file.path(package_root, "results", "DE_KCl_vs_Mannitol_final_count.tsv.gz"))
write_gz_table(conf_table, file.path(package_root, "results", "DE_final_conf.tsv.gz"))

normalization_out <- data.frame(
  sample = colnames(counts_fit),
  condition = samples[colnames(counts_fit), "condition"],
  replicate_label = samples[colnames(counts_fit), "replicate_label"],
  mapped_fragments = norm_input$mapped_fragments,
  mapped_fragment_size_factor = mapped_size_factors[colnames(counts_fit)],
  final_count_total = colSums(counts),
  final_conf_total = colSums(conf_counts),
  stringsAsFactors = FALSE
)
write.table(normalization_out, file.path(package_root, "qc", "normalization_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

filter_summary <- data.frame(
  stage = c("reported_union", "positive_final_count_in_at_least_1", sprintf("final_count_ge_%d_in_at_least_%d", min_count, min_samples)),
  loci = c(nrow(counts), sum(rowSums(counts) > 0), sum(keep))
)
write.table(filter_summary, file.path(package_root, "qc", "filter_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(package_root, "qc", "sessionInfo.txt"))

run_summary <- list(
  status = "complete",
  scope = "Locus-level Telescope differential expression",
  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  samples = colnames(counts_fit),
  design_mode = cfg[["design_mode"]],
  design_formula = paste(deparse(design_formula), collapse = ""),
  contrast = paste(cfg[["contrast_numerator"]], "vs", cfg[["contrast_denominator"]]),
  primary_count = "final_count",
  primary_normalization = "total genome-mapped fragments / geometric mean",
  filter = sprintf("final_count >= %d in at least %d samples", min_count, min_samples),
  loci_tested = nrow(counts_fit),
  alpha_for_report_and_independent_filter = alpha,
  hard_log2FC_cutoff_used = FALSE,
  lfc_shrinkage = cfg[["lfc_shrinkage"]],
  unshrunken_log2FoldChange_retained = TRUE,
  final_conf_comparison = TRUE,
  coefficient = primary$coefficient
)
write_json(run_summary, file.path(package_root, "qc", "ANALYSIS_RUN.json"), pretty = TRUE, auto_unbox = TRUE)

cat(toJSON(run_summary, pretty = TRUE, auto_unbox = TRUE), "\n")
