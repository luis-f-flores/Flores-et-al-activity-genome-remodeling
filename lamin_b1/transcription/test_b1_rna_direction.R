#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}
root<-commandArgs(TRUE)[1]
d<-read.delim(file.path(root,'76_gene_measurements.tsv'));reps<-read.delim(file.path(root,'76_gene_replicates.tsv'));out<-list()
cons<-c('Bic_vs_TTX','Bic_Nimo_vs_Bic','Bic_Nimo_vs_TTX');pairs<-list(c(1,2),c(2,3),c(1,3))
for(s in c('condition_depth')) for(k in seq_along(cons)) {
 con<-cons[k];g<-d[d$contrast==con & d$scheme==s,]
 rr<-reps[reps$scheme==s,];rr<-rr[match(g$gene_id,rr$gene_id),];stopifnot(identical(rr$gene_id,g$gene_id))
 mu<-cbind(rowMeans(rr[c('L9','L12')]),rowMeans(rr[c('L11','L14','L16')]),rowMeans(rr[c('L10','L13','L15')]))
 stopifnot(max(abs(mu-as.matrix(g[c('B1_TTX','B1_Bic','B1_Nimo')])))<1e-13)
 delta<-mu[,pairs[[k]][2]]-mu[,pairs[[k]][1]];stopifnot(max(abs(delta-g$B1_delta_log2))<1e-13)
 for(direction in c('decreased','increased')) for(rna in c('up','down')) {
  x<-g$B1_direction==direction;y<-g$RNA_classification==rna
  stopifnot(all(x==(if(direction=='decreased') delta<0 else delta>0)))
  cells<-c(sum(x&y),sum(x&!y),sum(!x&y),sum(!x&!y));ft<-fisher.test(matrix(cells,2,byrow=TRUE))
  out[[length(out)+1]]<-data.frame(contrast=con,scheme=s,B1_direction=direction,DEG_direction=rna,a=cells[1],b=cells[2],c=cells[3],d=cells[4],
   universe=nrow(g),B1_genes=sum(x),DEGs=sum(y),OR=unname(ft$estimate),CI_low=ft$conf.int[1],CI_high=ft$conf.int[2],P=ft$p.value)
 }
}
res<-do.call(rbind,out);res$BH_q<-NA_real_
for(s in unique(res$scheme)) {
 i<-res$scheme==s;res$BH_q[i]<-p.adjust(res$P[i],'BH')
}
write.table(res,file.path(root,'76_OR_results.tsv'),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
writeLines(capture.output(sessionInfo()),file.path(root,'76_R_session.txt'))
print(res)
