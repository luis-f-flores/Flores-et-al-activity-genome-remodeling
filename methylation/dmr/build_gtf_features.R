#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(dplyr)


args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("Usage: build_gtf_features.R RefGene.gtf output.rds")
gtf_file <- args[1]
dir.create(dirname(args[2]),recursive=TRUE,showWarnings=FALSE)
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")

tss <- promoters(txdb, upstream = 1000, downstream = 0, columns = c("tx_name", "gene_id"))


exons <- exons(txdb)
introns <- intronsByTranscript(txdb, use.names = TRUE)
five_utrs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
three_utrs <- threeUTRsByTranscript(txdb, use.names = TRUE)


genes <- genes(txdb)


saveRDS(list(
  promoters_plus1kb = tss,
  exons = exons,
  introns = introns,
  fiveUTRs = five_utrs,
  threeUTRs = three_utrs,
  genes = genes
), file = args[2])

message("Saved genomic features: ", args[2])
