#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript te_copy_proportions.R ",
    "<TE_DE_results.csv> <TE_annotation_master_mm10.tsv.gz> <output_dir>"
  )
}

de_file <- normalizePath(args[[1]], mustWork = TRUE)
catalog_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, mustWork = TRUE)

alpha <- 0.05
expected_de_targets <- 301L
expected_catalog_loci <- 1163808L
expected_catalog_targets <- 274L

normalize_target <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", trimws(as.character(x))))
}

write_tsv <- function(x, filename) {
  write.table(
    x,
    file = file.path(output_dir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

de <- read.csv(
  de_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA", "NaN")
)

required_de <- c("baseMean", "log2FoldChange", "padj")
if (!all(required_de %in% names(de))) {
  stop("DE table lacks required columns: ", paste(setdiff(required_de, names(de)), collapse = ", "))
}

target_col <- intersect(c("family", "target", "salmonte_target"), names(de))
if (length(target_col) == 0L) {
  de$salmonte_target <- as.character(de[[1L]])
} else {
  target_col <- target_col[[1]]
  de$salmonte_target <- as.character(de[[target_col]])
}

de$target_key <- normalize_target(de$salmonte_target)
if (anyNA(de$target_key) || any(de$target_key == "")) stop("Blank SalmonTE target name in DE table.")
if (anyDuplicated(de$salmonte_target)) stop("Duplicate literal SalmonTE target name in DE table.")
if (anyDuplicated(de$target_key)) stop("Target-name normalization creates a collision in the DE table.")
if (nrow(de) != expected_de_targets) {
  stop("Expected ", expected_de_targets, " SalmonTE targets; found ", nrow(de), ".")
}

de$baseMean <- suppressWarnings(as.numeric(de$baseMean))
de$log2FoldChange <- suppressWarnings(as.numeric(de$log2FoldChange))
de$padj <- suppressWarnings(as.numeric(de$padj))
de$status <- "Not significant"
de$status[!is.na(de$padj) & de$padj < alpha & de$log2FoldChange > 0] <- "Higher in KCl"
de$status[!is.na(de$padj) & de$padj < alpha & de$log2FoldChange < 0] <- "Lower in KCl"

catalog <- read.delim(
  catalog_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

required_catalog <- c(
  "locus_id", "salmonte_target", "raw_repeat_class", "raw_repeat_family"
)
if (!all(required_catalog %in% names(catalog))) {
  stop(
    "Canonical catalog lacks required columns: ",
    paste(setdiff(required_catalog, names(catalog)), collapse = ", ")
  )
}
if (anyDuplicated(catalog$locus_id)) stop("Canonical catalog contains duplicate locus_id values.")
if (nrow(catalog) != expected_catalog_loci) {
  stop("Expected ", expected_catalog_loci, " catalog loci; found ", nrow(catalog), ".")
}

catalog$target_key <- normalize_target(catalog$salmonte_target)
catalog_target_names <- unique(catalog[c("target_key", "salmonte_target")])
if (anyDuplicated(catalog_target_names$target_key)) {
  stop("Target-name normalization creates a collision in the canonical catalog.")
}
if (nrow(catalog_target_names) != expected_catalog_targets) {
  stop(
    "Expected ", expected_catalog_targets,
    " catalog-mapped SalmonTE targets; found ", nrow(catalog_target_names), "."
  )
}

de_index <- match(catalog$target_key, de$target_key)
if (anyNA(de_index)) {
  stop(
    "Canonical catalog contains target(s) absent from the DE table: ",
    paste(unique(catalog$salmonte_target[is.na(de_index)]), collapse = ", ")
  )
}

catalog$status <- de$status[de_index]
catalog$de_target <- de$salmonte_target[de_index]

status_levels <- c("Higher in KCl", "Lower in KCl", "Not significant")
copy_totals <- as.integer(table(factor(catalog$status, levels = status_levels)))
overall <- data.frame(
  status = status_levels,
  copies = copy_totals,
  percent = 100 * copy_totals / sum(copy_totals),
  stringsAsFactors = FALSE
)

target_copy_counts <- as.data.frame(table(catalog$target_key), stringsAsFactors = FALSE)
names(target_copy_counts) <- c("target_key", "catalog_copy_count")
target_copy_counts$catalog_copy_count <- as.integer(target_copy_counts$catalog_copy_count)

target_summary <- merge(
  de[c("salmonte_target", "target_key", "baseMean", "log2FoldChange", "padj", "status")],
  target_copy_counts,
  by = "target_key",
  all.x = TRUE,
  sort = FALSE
)
target_summary$catalog_copy_count[is.na(target_summary$catalog_copy_count)] <- 0L
target_summary <- target_summary[match(de$target_key, target_summary$target_key), ]

unmapped_targets <- target_summary[target_summary$catalog_copy_count == 0L, ]

composition <- aggregate(
  list(copies = catalog$locus_id),
  by = list(status = catalog$status, repeat_family = catalog$raw_repeat_family),
  FUN = length
)
composition$status_total <- ave(composition$copies, composition$status, FUN = sum)
composition$percent_within_status <- 100 * composition$copies / composition$status_total
composition <- composition[
  order(match(composition$status, status_levels), -composition$copies, composition$repeat_family),
]

expected_copies <- c(23957L, 4740L, 1135111L)
if (!identical(overall$copies, expected_copies)) {
  stop(
    "Verified Figure 5C copy totals changed. Expected ",
    paste(expected_copies, collapse = "/"), "; observed ",
    paste(overall$copies, collapse = "/"), "."
  )
}
if (nrow(unmapped_targets) != 27L) {
  stop("Expected 27 SalmonTE targets without a distinct catalog locus; found ", nrow(unmapped_targets), ".")
}

qc <- data.frame(
  metric = c(
    "de_targets", "catalog_loci", "catalog_mapped_targets",
    "de_targets_without_catalog_loci", "higher_expression_targets",
    "lower_expression_targets", "not_significant_targets",
    "higher_catalog_copies", "lower_catalog_copies", "not_significant_catalog_copies",
    "catalog_copy_total", "alpha", "copy_unit"
  ),
  value = c(
    nrow(de), nrow(catalog), nrow(catalog_target_names), nrow(unmapped_targets),
    sum(de$status == "Higher in KCl"), sum(de$status == "Lower in KCl"),
    sum(de$status == "Not significant"), overall$copies[[1]], overall$copies[[2]],
    overall$copies[[3]], sum(overall$copies), alpha,
    "raw mm10 RepeatMasker loci assigned to normalized SalmonTE targets"
  ),
  stringsAsFactors = FALSE
)

write_tsv(overall, "Fig5C_copy_proportions.tsv")
write_tsv(target_summary, "Fig5C_target_copy_counts.tsv")
write_tsv(unmapped_targets, "Fig5C_targets_without_catalog_loci.tsv")
write_tsv(composition, "Fig5C_repeat_family_composition.tsv")
write_tsv(qc, "Fig5C_QC.tsv")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

cat("TE genomic-copy summary\n")
print(overall, row.names = FALSE)
cat("\nOutput:", output_dir, "\n")
