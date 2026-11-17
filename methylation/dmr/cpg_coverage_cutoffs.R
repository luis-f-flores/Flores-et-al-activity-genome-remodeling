#!/usr/bin/env Rscript
# Calculate the per-sample upper coverage threshold used for DMR-union counts.
suppressPackageStartupMessages(library(methylKit))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript cpg_coverage_cutoffs.R samples.tsv output.tsv")
samples <- read.delim(args[1], stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("sample", "report") %in% names(samples)), !anyDuplicated(samples$sample))
sample_dir <- dirname(normalizePath(args[1]))
rows <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  report <- samples$report[i]
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", report)) report <- file.path(sample_dir, report)
  report <- normalizePath(report, mustWork = TRUE)
  x <- methRead(report, sample.id = samples$sample[i], assembly = "mm10",
                pipeline = "bismarkCytosineReport", mincov = 4)
  cutoff <- as.numeric(quantile(getData(x)$coverage, probs = .999))
  rows[[i]] <- data.frame(sample = samples$sample[i], report = report, upper_99.9 = cutoff)
  rm(x)
  gc()
}
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
write.table(do.call(rbind, rows), args[2], sep = "\t", quote = FALSE, row.names = FALSE)
