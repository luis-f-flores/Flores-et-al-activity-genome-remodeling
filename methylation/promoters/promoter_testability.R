#!/usr/bin/env Rscript
# Assess promoter overlap with tested methylation windows.
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({library(GenomicRanges); library(data.table)})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==5L)
features <- readRDS(args[1])
prom <- features[["promoters_plus1kb"]]
stopifnot(inherits(prom,"GRanges"), "gene_id" %in% names(mcols(prom)))
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
coords <- unique(data.frame(chr=as.character(seqnames(prom)), start=start(prom),
                           end=end(prom), gene_symbol=as.character(mcols(prom)$gene_id)))
stopifnot(all(coords$start>=1),all(coords$end>=coords$start))
fwrite(coords, file.path(out,"21_original_promoter_intervals.tsv"),sep="\t")

rows <- list()
for(i in seq_len(3L)) {
  tp <- c("1hr","6hr","24hr")[i]
  x <- fread(args[i+2L],select=c("chr","start","end","pvalue","qvalue"))
  x <- x[is.finite(pvalue)&is.finite(qvalue)]
  gr <- GRanges(x$chr,IRanges(x$start,x$end))
  hit <- countOverlaps(prom,gr,ignore.strand=TRUE)>0
  z <- data.table(gene_symbol=as.character(mcols(prom)$gene_id),tested=hit)
  z <- z[,.(wgbs_promoter_tested=any(tested)),by=gene_symbol]
  z[,timepoint:=tp]
  rows[[i]] <- z
  cat(tp,":",nrow(x),"tested windows;",sum(z$wgbs_promoter_tested),"genes with assessed promoters\n")
  rm(x,gr);gc()
}
fwrite(rbindlist(rows),file.path(out,"21_promoter_testability.tsv"),sep="\t")
writeLines(capture.output(sessionInfo()),file.path(out,"21_promoter_testability_session.txt"))
