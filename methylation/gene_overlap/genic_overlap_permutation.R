#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({library(dplyr);library(GenomicRanges)})
W <- gb_path("wgbs_package")
D <- gb_path("wgbs/RNA-seq_integration/All_DMR_q05_meth5_merged_overlap_DEGs")
RUN <- file.path(W,'runs/2026-10-14_RNA_DMR_main_supp');TAB<-file.path(RUN,'tables')
rd<-function(n)read.delim(file.path(TAB,n),check.names=FALSE,stringsAsFactors=FALSE)
wr<-function(d,n)write.table(d,file.path(TAB,n),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
cv<-rd('67_covariates_verified_loci.tsv')
cv$bin_length<-ntile(cv$gene_length,5)
cv$bin_cpg<-ntile(cv$true_cg_density_per_kb,5)
cv$bin_expr<-ntile(cv$baseline_fpkm,5)
cv$stratum<-paste(cv$bin_length,cv$bin_cpg,cv$bin_expr,sep='_')
stopifnot(!anyDuplicated(cv$gene_symbol),!anyNA(cv))
wr(cv,'67_matched_covariates_and_strata.tsv')
feat<-rd('67_promoter_body_union_intervals.tsv')
gr<-makeGRangesFromDataFrame(feat,keep.extra.columns=TRUE)
members<-rd('67_permutation_gene_membership.tsv')
for(cl in setdiff(names(members),c('timepoint','gene_symbol')))members[[cl]]<-members[[cl]]%in%c('True','TRUE',TRUE)
N<-10000L;RNGkind('Mersenne-Twister','Inversion','Rejection')
summaries<-list();nulls<-list();stratas<-list();examples<-list();verify<-list();balance<-list()
for(t in c(1L,6L,24L)){
 tp<-paste0(t,'hr');seed<-20261018L+t;set.seed(seed)
 m<-members[members$timepoint==tp,];ix<-match(cv$gene_symbol,m$gene_symbol)
 stopifnot(!anyNA(ix),all(m$permutation_eligible[ix]))

 dm<-do.call(rbind,lapply(c('hyper','hypo'),function(di)read.delim(file.path(D,paste0(t,'M_vs_',t,'K_Merged_DMRs_cov4_',di,'.txt')))))
 dr<-makeGRangesFromDataFrame(dm)
 hits<-unique(feat$gene_symbol[queryHits(findOverlaps(gr,dr,ignore.strand=TRUE))])
 stopifnot(identical(sort(intersect(cv$gene_symbol,hits)),sort(cv$gene_symbol[m$is_promoter_body_DMR[ix]])))
 y<-cv$gene_symbol%in%hits
 deg<-m$is_DEG[ix]
 obs<-sum(deg&y)
 stopifnot(obs==sum(m$eligible_observed),obs==sum(m$corrected_main&m$permutation_eligible))
 pools<-split(seq_len(nrow(cv)),cv$stratum);take<-vapply(pools,function(v)sum(deg[v]),integer(1))
 st<-data.frame(timepoint=tp,stratum=names(pools),N=vapply(pools,length,integer(1)),K=vapply(pools,function(v)sum(y[v]),integer(1)),n=take,row.names=NULL)
 st$expected<-st$n*st$K/st$N
 st$variance<-with(st,ifelse(N>1,n*(K/N)*(1-K/N)*(N-n)/(N-1),0))
 stopifnot(sum(st$n)==sum(deg),all(st$n<=st$N))
 stratas[[tp]]<-st
 values<-integer(N);first<-NULL;sum_counts<-numeric(nrow(cv))
 for(i in seq_len(N)){
  draw<-unlist(lapply(names(pools)[take>0],function(s){v<-pools[[s]];v[sample.int(length(v),take[[s]],replace=FALSE)]}),use.names=FALSE)
  stopifnot(length(draw)==sum(deg),!anyDuplicated(draw))
  values[i]<-sum(y[draw]);sum_counts[draw]<-sum_counts[draw]+1
  if(i==1L)first<-draw
 }
 p<-(sum(values>=obs)+1)/(N+1)
 summaries[[tp]]<-data.frame(timepoint=tp,background_n=nrow(cv),DEG_n=sum(deg),observed=obs,null_mean=mean(values),null_sd=sd(values),empirical_p=p,seed=seed,draws=N,expected_exact=sum(st$expected),sd_exact=sqrt(sum(st$variance)),observed_over_expected=obs/sum(st$expected))
 nulls[[tp]]<-data.frame(timepoint=tp,draw=seq_len(N),overlap=values)
 examples[[tp]]<-data.frame(timepoint=tp,gene_symbol=cv$gene_symbol[first],stratum=cv$stratum[first],DMR_hit=y[first])

 for(col in c('gene_length','true_cg_density_per_kb','baseline_fpkm')){
  z<-log1p(cv[[col]])
  balance[[paste(tp,col)]]<-data.frame(timepoint=tp,covariate=col,DEG_log1p_mean=mean(z[deg]),null_log1p_mean=sum(z*sum_counts)/(N*sum(deg)),background_log1p_mean=mean(z),DEG_raw_median=median(cv[[col]][deg]),first_null_raw_median=median(cv[[col]][first]))
 }
 verify[[tp]]<-data.frame(timepoint=tp,independent_GenomicRanges_membership=TRUE,observed_eligible_n=obs,all_null_draws_same_n=TRUE,all_draws_without_replacement=TRUE,exact_stratum_counts=TRUE)
 print(summaries[[tp]]);flush.console()
}
wr(bind_rows(summaries),'67_matched_permutation_summary.tsv')
wr(bind_rows(nulls),'67_matched_permutation_draws.tsv')
wr(bind_rows(stratas),'67_permutation_stratum_accounting.tsv')
wr(bind_rows(examples),'67_first_null_draw_genes.tsv')
wr(bind_rows(verify),'67_independent_R_endpoint_QC.tsv')
wr(bind_rows(balance),'67_covariate_balance.tsv')
writeLines(capture.output(sessionInfo()),file.path(RUN,'67_R_sessionInfo.txt'))
