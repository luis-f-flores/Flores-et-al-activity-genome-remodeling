#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

load_collapsed_gene_models <- function(gene_body_bed_file, tss_bed_file = NULL) {
  gene_body <- rtracklayer::import(gene_body_bed_file, format = "BED")
  gene_body$gene_symbol <- gene_body$name
  gene_body$gene_id <- gene_body$name

  dup_symbols <- gene_body$gene_symbol[duplicated(gene_body$gene_symbol)]
  if (length(dup_symbols) > 0) {
    message(sprintf(
      "load_collapsed_gene_models(): excluding %d gene symbols with >1 distinct locus in %s (ambiguous join key, not a single gene).",
      length(unique(dup_symbols)), basename(gene_body_bed_file)
    ))
    gene_body <- gene_body[!gene_body$gene_symbol %in% dup_symbols]
  }

  coord_key <- paste(GenomicRanges::seqnames(gene_body), GenomicRanges::start(gene_body),
                      GenomicRanges::end(gene_body))
  dup_coords <- coord_key[duplicated(coord_key)]
  if (length(dup_coords) > 0) {
    message(sprintf(
      "load_collapsed_gene_models(): excluding %d gene symbols sharing a locus with >1 distinct gene_symbol in %s (ambiguous join key on the coordinate axis).",
      sum(coord_key %in% unique(dup_coords)), basename(gene_body_bed_file)
    ))
    gene_body <- gene_body[!coord_key %in% dup_coords]
  }

  if (!is.null(tss_bed_file) && file.exists(tss_bed_file)) {
    tss <- rtracklayer::import(tss_bed_file, format = "BED")
    tss$gene_symbol <- tss$name
    tss <- tss[tss$gene_symbol %in% gene_body$gene_symbol]
    tss$gene_id <- tss$name
  } else {
    tss <- GenomicRanges::resize(gene_body, width = 1, fix = "start")
  }

  list(gene_body = gene_body, tss = tss)
}

build_promoter_windows <- function(tss_gr, window = c(-1000, 0)) {
  GenomicRanges::promoters(tss_gr, upstream = -window[1], downstream = window[2])
}


get_gene_body_windows <- function(gene_models) {
  gene_models$gene_body
}

overlap_regions_to_genes <- function(region_df, gene_windows) {
  region_gr <- GenomicRanges::GRanges(
    seqnames = region_df$chr,
    ranges   = IRanges::IRanges(start = region_df$start, end = region_df$end)
  )
  hits <- GenomicRanges::findOverlaps(region_gr, gene_windows, minoverlap = 1L)

  out <- region_df[S4Vectors::queryHits(hits), , drop = FALSE]
  out$gene_id     <- gene_windows$gene_id[S4Vectors::subjectHits(hits)]
  out$gene_symbol <- gene_windows$gene_symbol[S4Vectors::subjectHits(hits)]
  out
}
