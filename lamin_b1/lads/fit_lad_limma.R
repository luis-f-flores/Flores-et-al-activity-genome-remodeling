#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}

suppressPackageStartupMessages(library(limma))
run<-commandArgs(TRUE)[1];tab<-file.path(run,'tables')
dir.create(file.path(run,'qc'),recursive=TRUE,showWarnings=FALSE)
mat<-function(n) as.matrix(read.table(file.path(tab,n),sep='\t',header=FALSE))
wr<-function(x,n) write.table(x,file.path(tab,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
group<-factor(c('TTX','TTX','Bic','Bic','Bic','Nimo','Nimo','Nimo'),levels=c('TTX','Bic','Nimo'))
design<-model.matrix(~0+group);colnames(design)<-levels(group)
contrast<-makeContrasts(Bic_vs_TTX=Bic-TTX,Bic_Nimo_vs_Bic=Nimo-Bic,Bic_Nimo_vs_TTX=Nimo-TTX,levels=design)
lib<-as.numeric(mat('64_library_sizes.tsv'));summaries<-list();n<-0
for(con in colnames(contrast)) {
prefix<-paste0('64_',con)
counts<-mat(paste0(prefix,'_counts.tsv'));stopifnot(ncol(counts)==8,all(counts>0))
y<-sweep(log2(counts+.5),2,log2(lib+1),'-')+log2(1e6)
cf<-lmFit(y,design)
fitted_count<-sweep(2^(cf$coefficients%*%t(design)),2,(lib+1)/1e6,'*')
sx<-rowMeans(y)+mean(log2(lib+1))-log2(1e6)
for(scheme in 'condition_depth') {
 value<-mat(paste0(prefix,'_',scheme,'.tsv'));initial<-lmFit(value,design)
 stopifnot(all(initial$df.residual==5))
 sy<-sqrt(initial$sigma);curve<-lowess(sx,sy,f=.7);fn<-approxfun(curve,rule=2)
 predicted<-fn(log2(fitted_count));stopifnot(all(is.finite(predicted)),all(predicted>0))
 w<-matrix(1/predicted^4,nrow=nrow(counts),ncol=8)
 fit<-eBayes(contrasts.fit(lmFit(value,design,weights=w),contrast))
 stopifnot(all(is.finite(fit$p.value)),all(is.finite(fit$coefficients)))
 gg<-list(1:2,3:5,6:8)
 means<-sapply(gg,function(g) rowSums(value[,g,drop=FALSE]*w[,g,drop=FALSE])/rowSums(w[,g,drop=FALSE]))
 s2<-rowSums(w*(value-means[,c(1,1,2,2,2,3,3,3)])^2)/5
 post<-if(is.infinite(fit$df.prior)) rep(fit$s2.prior,nrow(value)) else (fit$df.prior*fit$s2.prior+5*s2)/(fit$df.prior+5)
 stopifnot(max(abs(s2-fit$sigma^2))<1e-10,max(abs(post-fit$s2.post))<1e-10)
 obj<-list(counts=counts,E=value,weights=w,design=design,contrasts=contrast,fit=fit,sx=sx,sy=sy,curve=curve,fitted_count=fitted_count)
 saveRDS(obj,file.path(tab,paste0(prefix,'_',scheme,'_model.rds')))
 wr(data.frame(log2_count=sx,sqrt_SD=sy,smoothed=fn(sx)),paste0(prefix,'_',scheme,'_trend.tsv'))
 write.table(w,file.path(tab,paste0(prefix,'_',scheme,'_weights.tsv')),sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
 comps<-list(c(1,2),c(2,3),c(1,3));maxse<-0
 for(j in match(con,colnames(contrast))) {
  se<-fit$stdev.unscaled[,j]*sqrt(fit$s2.post);ci<-qt(.975,df=fit$df.total)*se
  g1<-gg[[comps[[j]][1]]];g2<-gg[[comps[[j]][2]]]
  manualSE<-sqrt(post*(1/rowSums(w[,g1,drop=FALSE])+1/rowSums(w[,g2,drop=FALSE])))
  maxse<-max(maxse,max(abs(manualSE-se)));stopifnot(maxse<1e-10)
  out<-data.frame(row_index=seq_len(nrow(value))-1,mean_z=rowMeans(value),delta_z=fit$coefficients[,j],SE=se,CI_low=fit$coefficients[,j]-ci,CI_high=fit$coefficients[,j]+ci,t=fit$t[,j],df=fit$df.total,P=fit$p.value[,j],q=p.adjust(fit$p.value[,j],'BH'))
  wr(out,paste0(prefix,'_',scheme,'_tests.tsv'))
 }
 n<-n+1;summaries[[n]]<-data.frame(unit='pairwise_segment',contrast=con,scheme=scheme,features=nrow(value),residual_df=5,prior_df=fit$df.prior,weight_min=min(w),weight_median=median(w),weight_max=max(w),manual_SE_max_error=maxse,manual_variance_max_error=max(abs(s2-fit$sigma^2)),manual_posterior_max_error=max(abs(post-fit$s2.post)))
 cat(scheme,'complete\n');flush.console()
}
}
wr(do.call(rbind,summaries),'64_model_summary.tsv')
writeLines(capture.output(sessionInfo()),file.path(run,'qc','64_R_session.txt'))
