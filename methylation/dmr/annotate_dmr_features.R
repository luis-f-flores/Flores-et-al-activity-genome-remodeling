#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


library(GenomicRanges)
library(dplyr)
library(stringr)


root_dir <- gb_path("wgbs/methylKit_DMR")
gtf_rds <- gb_path("wgbs/gtf_dmrs_features_plus1kb.rds")
pattern <- "DMR_.+_multiwindow/DMR_win_500bp$"


cat("Loading GTF-based features from:\n", gtf_rds, "\n\n")
gtf_features <- readRDS(gtf_rds)


annotate_dmr_file <- function(file_path) {
  file_name <- basename(file_path)
  working_dir <- dirname(file_path)
  setwd(working_dir)
  cat("Annotating:", file_name, "\n")

  dmrs <- read.delim(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("chr", "start", "end") %in% colnames(dmrs))) {
    cat("Skipping", file_name, "— missing required columns.\n")
    return()
  }

  gr <- GRanges(seqnames = dmrs$chr, ranges = IRanges(start = dmrs$start, end = dmrs$end))
  all_hits <- list()
  summary_log <- data.frame(Feature = character(), Overlap_Count = integer())

  for (feature_name in names(gtf_features)) {
    feature_set <- gtf_features[[feature_name]]
    if (inherits(feature_set, "GRangesList")) {
      feature_set <- unlist(feature_set, use.names = FALSE)
    }
    if (length(feature_set) == 0) next

    hits <- findOverlaps(gr, feature_set, minoverlap = 1)
    hit_count <- length(hits)

    if (hit_count > 0) {
      cat(sprintf("%s overlaps found: %d\n", feature_name, hit_count))
      summary_log <- rbind(summary_log, data.frame(Feature = feature_name, Overlap_Count = hit_count))

      mcols_data <- mcols(feature_set)[subjectHits(hits), ]
      annot_df <- data.frame(
        chr = as.character(seqnames(gr[queryHits(hits)])),
        start = start(gr[queryHits(hits)]),
        end = end(gr[queryHits(hits)]),
        strand = as.character(strand(feature_set)[subjectHits(hits)]),
        feature = feature_name,
        tx_id = if ("tx_name" %in% colnames(mcols_data)) mcols_data$tx_name else NA,
        gene_name = if ("gene_id" %in% colnames(mcols_data)) sapply(mcols_data$gene_id, function(x) paste(x, collapse = ",")) else NA
      )

      all_hits[[feature_name]] <- annot_df
    }
  }


annotated_df <- bind_rows(all_hits)


dmrs$key <- paste(dmrs$chr, dmrs$start, dmrs$end, sep = "_")
annotated_keys <- unique(paste(annotated_df$chr, annotated_df$start, annotated_df$end, sep = "_"))


intergenic_df <- dmrs[!dmrs$key %in% annotated_keys, c("chr", "start", "end")]
intergenic_df$strand <- "*"
intergenic_df$feature <- "intergenic"
intergenic_df$tx_id <- NA
intergenic_df$gene_name <- NA


all_hits[["intergenic"]] <- intergenic_df
summary_log <- rbind(summary_log, data.frame(Feature = "intergenic", Overlap_Count = nrow(intergenic_df)))


if (length(all_hits) > 0) {
  final_df <- bind_rows(all_hits)
  out_annot <- file.path(getwd(), paste0("GTF_Annot_", file_name))
  out_summary <- file.path(getwd(), paste0("Summary_Overlap_", file_name, ".txt"))

  write.table(final_df, out_annot, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(summary_log, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Annotation saved to:", out_annot, "\n")
  cat("Summary saved to:", out_summary, "\n")
} else {
  cat("No overlapping features found for:", file_name, "\n")
}
}


target_dirs <- list.dirs(root_dir, recursive = TRUE, full.names = TRUE)
filtered_dirs <- target_dirs[grepl(pattern, target_dirs)]

cat("Multiwindow DMR directories found:", length(filtered_dirs), "\n\n")

for (dir in filtered_dirs) {
  cat("Checking directory:", dir, "\n")
  dmr_files <- list.files(dir, pattern = "Merged_DMRs_cov.*_(hyper|hypo)\\.txt$", full.names = TRUE)
  cat("Files found:", length(dmr_files), "\n")

  if (length(dmr_files) == 0) {
    cat(" No DMR files found in this directory. Skipping...\n\n")
    next
  }

  for (file in dmr_files) {
    cat("Calling annotation function on:", file, "\n")
    try(annotate_dmr_file(file), silent = FALSE)
  }
}

cat("DONE: All GTF-based DMR annotations (with strand) complete.\n")

