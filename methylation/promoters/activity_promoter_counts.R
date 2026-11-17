#!/usr/bin/env Rscript
# Measure activity-gene promoter methylation from per-sample CpG counts.
suppressPackageStartupMessages({library(methylKit);library(GenomicRanges);library(rtracklayer);library(data.table)})
setDTthreads(2)
a<-commandArgs(trailingOnly=TRUE)
if(length(a)<5)stop('Usage: activity_promoter_counts.R gene_layers_dir activity_membership.tsv samples.tsv output_dir annotation_helpers.R')
source(a[5]);models<-load_collapsed_gene_models(file.path(a[1],'genes.bed'),file.path(a[1],'tss.bed'))
classes<-read.delim(a[2],stringsAsFactors=FALSE);samples<-read.delim(a[3],stringsAsFactors=FALSE);out<-a[4]
stopifnot(all(c('gene_name','gene_set')%in%names(classes)),all(c('sample','report')%in%names(samples)),!anyDuplicated(samples$sample))
classes<-data.frame(gene=classes$gene_name,gene_class=classes$gene_set)
gr<-build_promoter_windows(models$tss,window=c(-1000,100));gr<-gr[gr$gene_symbol%in%classes$gene]
stopifnot(length(gr)==173L,!anyDuplicated(gr$gene_symbol))
# The activity-gene display uses the defined Bdnf exon IV interval.
bdnf<-which(gr$gene_symbol=='Bdnf');stopifnot(length(bdnf)==1,as.character(seqnames(gr))[bdnf]=='chr2')
ranges(gr)[bdnf]<-IRanges(109692289,109692729)
annotation<-data.frame(gene=gr$gene_symbol,gene_class=classes$gene_class[match(gr$gene_symbol,classes$gene)],chr=as.character(seqnames(gr)),start=start(gr),end=end(gr),strand=as.character(strand(gr)))
dir.create(out,recursive=TRUE,showWarnings=FALSE);write.table(annotation,file.path(out,'promoter_regions_1based.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
result<-list()
for(i in seq_len(nrow(samples))) {
 id<-samples$sample[i];cat('Reading',id,'\n');flush.console()
 x<-methRead(samples$report[i],sample.id=id,assembly='mm10',pipeline='bismarkCytosineReport',mincov=4)
 x<-filterByCoverage(x,hi.perc=99.9)
 rr<-regionCounts(x,regions=gr,cov.bases=0,strand.aware=FALSE);d<-getData(rr)
 hit<-findOverlaps(GRanges(d$chr,IRanges(d$start,d$end)),gr,type='equal',select='first')
 stopifnot(!anyNA(hit),!anyDuplicated(hit))
 d$gene<-gr$gene_symbol[hit];d$gene_class<-classes$gene_class[match(d$gene,classes$gene)];d$sample<-id;d$mCG_pct<-100*d$numCs/d$coverage
 result[[id]]<-d;write.table(d,file.path(out,paste0(id,'_promoter_counts.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
 rm(x,rr,d);gc()
}
write.table(do.call(rbind,result),file.path(out,'activity_promoter_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
