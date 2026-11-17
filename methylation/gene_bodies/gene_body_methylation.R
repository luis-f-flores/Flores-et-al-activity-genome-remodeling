#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(methylKit);library(GenomicRanges);library(data.table)})

read_gene_bodies <- function(path) {
  annotation <- read.delim(path, header=FALSE,
    col.names=c('chr','start0','end','locus_id','score','strand','gene_symbol'))
  stopifnot(!anyDuplicated(annotation$locus_id))
  coordinates <- paste(annotation$chr, annotation$start0, annotation$end, sep=':')
  ambiguous <- duplicated(coordinates) | duplicated(coordinates, fromLast=TRUE)
  annotation <- annotation[!ambiguous, ]
  regions <- GRanges(annotation$chr, IRanges(annotation$start0+1L, annotation$end),
                     strand=annotation$strand)
  regions$locus_id <- annotation$locus_id
  regions$gene_symbol <- annotation$gene_symbol
  regions
}

test_gene_bodies <- function(counts, regions, timepoint, cores=2L) {
  united <- unite(counts, destrand=FALSE, min.per.group=2L)
  set.seed(20261007)
  result <- calculateDiffMeth(united, test='Chisq', adjust='SLIM',
                              effect='wmean', mc.cores=cores, save.db=FALSE)
  values <- getData(result)
  hit <- findOverlaps(GRanges(values$chr, IRanges(values$start, values$end)),
                      regions, type='equal', select='first')
  stopifnot(!anyNA(hit))
  values$locus_id <- regions$locus_id[hit]
  values$gene_symbol <- regions$gene_symbol[hit]
  values$timepoint <- timepoint
  values
}

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args)!=4L) stop('Usage: gene_body_methylation.R timepoint samples.tsv gene_bodies.bed output_directory')
  tp <- args[1]
  stopifnot(tp %in% c('1hr','6hr','24hr'))
  samples <- read.delim(args[2], stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(samples$sample))
  ids <- c(paste0(LETTERS[1:3], sub('hr$','',tp), 'M'),
           paste0(LETTERS[1:3], sub('hr$','',tp), 'K'))
  selected <- samples[match(ids,samples$sample), ]
  stopifnot(!anyNA(selected$sample))
  cores <- as.integer(Sys.getenv('GB_CORES', '2'))
  stopifnot(!is.na(cores), cores>=1L)
  data.table::setDTthreads(cores)
  regions <- read_gene_bodies(args[3])
  raw <- methRead(as.list(selected$report), sample.id=as.list(ids), assembly='mm10',
    treatment=c(0,0,0,1,1,1), pipeline='bismarkCytosineReport', mincov=4)
  raw <- filterByCoverage(raw, hi.perc=99.9)
  raw <- normalizeCoverage(raw)
  counts <- regionCounts(raw, regions=regions)
  values <- test_gene_bodies(counts, regions, tp, cores)
  output <- args[4]
  dir.create(output, recursive=TRUE, showWarnings=FALSE)
  write.csv(values, file.path(output, sprintf('08_genebody_methylkit_direct_%s_all.csv',tp)), row.names=FALSE)
  significant <- values[abs(values$meth.diff)>=5 & values$qvalue<=.05, ]
  write.csv(significant, file.path(output, sprintf('08_genebody_significant_%s.csv',tp)), row.names=FALSE)
  cat(tp, ':', nrow(values), 'tested gene bodies;', nrow(significant), 'significant\n')
}

if (sys.nframe()==0L) main()
