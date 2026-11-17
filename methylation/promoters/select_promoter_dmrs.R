#!/usr/bin/env Rscript
# Assign one representative promoter DMR per gene, then apply q01 and 25-point thresholds.
suppressPackageStartupMessages({library(GenomicRanges);library(dplyr);library(data.table)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: select_promoter_dmrs.R features.rds merged_DMR_manifest.tsv output_directory")
features<-readRDS(args[1]);prom<-features$promoters_plus1kb
manifest<-fread(args[2]);out<-args[3]
stopifnot(all(c('timepoint','direction','path')%in%names(manifest)),all(width(prom)==1000),nrow(manifest)==6L)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
for(tp in c(1,6,24)) {
 parts<-list()
 for(dr in c('hypo','hyper')) {
  path<-manifest[timepoint==tp & direction==dr,path]
  stopifnot(length(path)==1L)
  d<-fread(path);gr<-GRanges(d$chr,IRanges(d$start,d$end));h<-findOverlaps(gr,prom,minoverlap=1)
  parts[[dr]]<-data.frame(gene_name=as.character(mcols(prom)$gene_id[subjectHits(h)]),
      meth_diff=d$meth_diff[queryHits(h)],qvalue=d$qvalue[queryHits(h)],DMR_type=dr,
      prom_chr=d$chr[queryHits(h)],prom_start=d$start[queryHits(h)],prom_end=d$end[queryHits(h)])
 }
 selected<-bind_rows(parts)%>%filter(!is.na(gene_name))%>%group_by(gene_name)%>%slice_min(qvalue,with_ties=FALSE)%>%ungroup()
 stopifnot(!anyDuplicated(selected$gene_name))
 for(direction in c('hyper','hypo')) {
  x<-selected[selected$DMR_type==direction,];stem<-paste0(tp,'M_vs_',tp,'K_',direction,'_Promoter_DMRs_ForVolcano')
  fwrite(x,file.path(out,paste0(stem,'.txt')),sep='\t')
  x<-x[x$qvalue<=.01 & abs(x$meth_diff)>=25,]
  fwrite(x,file.path(out,paste0(stem,'_meth25_q01.txt')),sep='\t')
 }
}
