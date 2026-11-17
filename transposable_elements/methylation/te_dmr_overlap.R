#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop("Usage: te_dmr_overlap.R TE_mm10.bed TE_meta.tsv merged_DMR_manifest.tsv output_directory")
te_catalog_bed_file<-args[1];te_catalog_meta_file<-args[2];manifest<-read.delim(args[3],stringsAsFactors=FALSE)
OUTPUT_TABLES_DIR<-args[4];dir.create(OUTPUT_TABLES_DIR,recursive=TRUE,showWarnings=FALSE)
TIMEPOINTS<-c('1hr','6hr','24hr')
dmr_merged_files<-setNames(lapply(c(1,6,24),function(t){x<-manifest[manifest$timepoint==t,];stopifnot(nrow(x)==2L);setNames(as.list(x$path),x$direction)}),TIMEPOINTS)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(dplyr)
})

dmr_all <- do.call(rbind, lapply(TIMEPOINTS, function(tp) {
  hyper <- utils::read.delim(dmr_merged_files[[tp]]$hyper, stringsAsFactors = FALSE)
  hypo  <- utils::read.delim(dmr_merged_files[[tp]]$hypo, stringsAsFactors = FALSE)
  df <- rbind(hyper, hypo)
  df <- dplyr::rename(df, meth.diff = meth_diff)
  df$timepoint <- tp
  df
}))

dmr_gr <- GenomicRanges::GRanges(
  seqnames = dmr_all$chr,
  ranges   = IRanges::IRanges(start = dmr_all$start, end = dmr_all$end),
  timepoint = dmr_all$timepoint,
  meth.diff = dmr_all$meth.diff
)

te_catalog_bed  <- rtracklayer::import(te_catalog_bed_file, format = "BED")
te_catalog_meta <- utils::read.delim(te_catalog_meta_file, stringsAsFactors = FALSE)

te_catalog_bed$TE_family <- te_catalog_meta$family[match(te_catalog_bed$name, te_catalog_meta$copyID)]

hits <- GenomicRanges::findOverlaps(dmr_gr, te_catalog_bed, minoverlap = 1L)

overlap_df <- data.frame(
  dmr_chr      = as.character(GenomicRanges::seqnames(dmr_gr))[S4Vectors::queryHits(hits)],
  dmr_start    = GenomicRanges::start(dmr_gr)[S4Vectors::queryHits(hits)],
  dmr_end      = GenomicRanges::end(dmr_gr)[S4Vectors::queryHits(hits)],
  timepoint    = dmr_gr$timepoint[S4Vectors::queryHits(hits)],
  meth.diff    = dmr_gr$meth.diff[S4Vectors::queryHits(hits)],
  TE_copy_id   = te_catalog_bed$name[S4Vectors::subjectHits(hits)],
  TE_family    = te_catalog_bed$TE_family[S4Vectors::subjectHits(hits)]
)
overlap_df$direction <- ifelse(overlap_df$meth.diff > 0, "hyper", "hypo")

write.csv(overlap_df,
          file.path(OUTPUT_TABLES_DIR, "07_dmr_te_locus_overlap.csv"),
          row.names = FALSE)

family_summary <- overlap_df %>%
  dplyr::group_by(TE_family, timepoint, direction) %>%
  dplyr::summarise(n_copies = dplyr::n_distinct(TE_copy_id), .groups = "drop")

write.csv(family_summary,
          file.path(OUTPUT_TABLES_DIR, "07_te_family_dmr_summary.csv"),
          row.names = FALSE)

family_totals <- te_catalog_meta %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(total_copies = dplyr::n(), total_bp = sum(length, na.rm = TRUE), .groups = "drop") %>%
  dplyr::rename(TE_family = family)

family_summary_norm <- family_summary %>%
  dplyr::left_join(family_totals, by = "TE_family") %>%
  dplyr::mutate(
    pct_copies_overlapping = 100 * n_copies / total_copies,
    overlap_per_mb = n_copies / (total_bp / 1e6)
  )

write.csv(family_summary_norm,
          file.path(OUTPUT_TABLES_DIR, "07_te_family_dmr_normalized.csv"),
          row.names = FALSE)

cat("TE locus and family DMR-overlap tables written to", OUTPUT_TABLES_DIR, "\n")
