#!/usr/bin/env Rscript
# Assign promoter, gene-body and intergenic DMR classes and count CpG-island overlap.
suppressPackageStartupMessages({library(readr);library(dplyr);library(tibble);library(GenomicRanges);library(GenomicFeatures);library(GenomeInfoDb);library(rtracklayer);library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3)stop('Usage: summarize_dmr_annotations.R RefGene.gtf DMR_files.tsv output_dir [CpG_islands.bed]')
gtf_src<-args[1];files<-read.delim(args[2],stringsAsFactors=FALSE);out<-args[3]
stopifnot(all(c('timepoint','direction','path')%in%names(files)),nrow(files)==6L)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
txdb <- suppressWarnings(suppressMessages(makeTxDbFromGFF(gtf_src, format = "gtf")))
gtf_gr <- import(gtf_src)
tx_rows <- gtf_gr[mcols(gtf_gr)$type == "transcript"]
tx2gene <- tibble(
  tx_name = as.character(mcols(tx_rows)$transcript_id),
  gene    = as.character(mcols(tx_rows)$gene_name)
) %>% distinct() %>% filter(!is.na(tx_name), !is.na(gene))

tx <- transcripts(txdb, use.names = TRUE)
mcols(tx)$tx_name <- if (!is.null(mcols(tx)$tx_name)) mcols(tx)$tx_name else names(tx)

prom <- promoters(tx, upstream = 1000, downstream = 0)
mcols(prom)$tx_name <- mcols(tx)$tx_name
prom <- prom[!is.na(mcols(prom)$tx_name)]
prom$gene <- tx2gene$gene[match(prom$tx_name, tx2gene$tx_name)]

ex_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
exons_tx <- unlist(ex_tx); exons_tx$tx_name <- rep(names(ex_tx), elementNROWS(ex_tx))
exons_tx$gene <- tx2gene$gene[match(exons_tx$tx_name, tx2gene$tx_name)]

utr5_tx <- fiveUTRsByTranscript(txdb, use.names = TRUE)
utr5 <- unlist(utr5_tx); utr5$tx_name <- rep(names(utr5_tx), elementNROWS(utr5_tx))
utr5$gene <- tx2gene$gene[match(utr5$tx_name, tx2gene$tx_name)]

utr3_tx <- threeUTRsByTranscript(txdb, use.names = TRUE)
utr3 <- unlist(utr3_tx); utr3$tx_name <- rep(names(utr3_tx), elementNROWS(utr3_tx))
utr3$gene <- tx2gene$gene[match(utr3$tx_name, tx2gene$tx_name)]

intr_tx <- intronsByTranscript(txdb, use.names = TRUE)
intr <- unlist(intr_tx); intr$tx_name <- rep(names(intr_tx), elementNROWS(intr_tx))
intr$gene <- tx2gene$gene[match(intr$tx_name, tx2gene$tx_name)]

gene_features <- c(prom, exons_tx, utr5, utr3, intr)
gene_features <- gene_features[!is.na(mcols(gene_features)$gene)]
gene_features <- keepStandardChromosomes(gene_features, pruning.mode = "coarse")

if(length(args)>=4) {
 cgi <- import(args[4],format='BED')
} else {
 cgi <- annotatr::build_annotations(genome='mm10',annotations='mm10_cpg_islands')
}
body_features<-c(exons_tx,utr5,utr3,intr)
rows<-list()
for(i in seq_len(nrow(files))) {
 d<-read.delim(files$path[i],stringsAsFactors=FALSE)
 gr<-makeGRangesFromDataFrame(d,seqnames.field='chr',keep.extra.columns=TRUE)
 # CGI overlaps are counted separately from the three genomic classes.
 cgi_hits<-unique(queryHits(findOverlaps(gr,cgi,ignore.strand=TRUE)))
 n_cgi<-length(unique(paste0(as.character(seqnames(gr))[cgi_hits],':',start(gr)[cgi_hits],'-',end(gr)[cgi_hits])))
 gr<-keepStandardChromosomes(gr,pruning.mode='coarse')
 common<-intersect(seqlevelsInUse(gr),seqlevelsInUse(gene_features))
 gr<-keepSeqlevels(gr,common,pruning.mode='coarse')
 classes<-rep('intergenic_dmrs',length(gr))
 classes[unique(queryHits(findOverlaps(gr,body_features,ignore.strand=TRUE)))]<-'gene_body_dmrs'
 classes[unique(queryHits(findOverlaps(gr,prom,ignore.strand=TRUE)))]<-'promoter_dmrs'
 counts<-as.list(setNames(as.integer(table(factor(classes,levels=c('promoter_dmrs','gene_body_dmrs','intergenic_dmrs')))),c('promoter_dmrs','gene_body_dmrs','intergenic_dmrs')))
 counts$CGI<-n_cgi;total<-sum(unlist(counts))
 rows[[i]]<-list(timepoint=files$timepoint[i],direction=files$direction[i],unique_genomic_DMRs=length(gr),CGI_overlap_count=n_cgi,sum_of_four_categories=total,displayed_total=total,counts=counts,percentages=as.list(round(100*unlist(counts)/total,2)))
}
write_json(rows,file.path(out,'original_DMR_annotation_pies.json'),auto_unbox=TRUE,pretty=TRUE,digits=15)
