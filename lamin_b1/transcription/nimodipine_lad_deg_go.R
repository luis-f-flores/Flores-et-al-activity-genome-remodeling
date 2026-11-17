#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}

args<-commandArgs(trailingOnly=TRUE);D<-args[1];options(stringsAsFactors=FALSE);set.seed(47)
suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db);library(GO.db);library(AnnotationDbi);library(jsonlite)})
stopifnot(as.character(packageVersion('clusterProfiler'))=='4.6.2',as.character(packageVersion('org.Mm.eg.db'))=='3.16.0',as.character(packageVersion('GO.db'))=='3.16.0')
g<-read.delim(file.path(D,'92_GO_membership_input.tsv'),check.names=FALSE,quote='')
for(n in c('DEG','baseline_LAD')){s<-tolower(as.character(g[[n]]));stopifnot(all(s %in% c('true','false')));g[[n]]<-s=='true'}
wt<-function(x,n)write.table(x,file.path(D,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
mm<-AnnotationDbi::select(org.Mm.eg.db,keys=unique(g$gene_id),keytype='ENSEMBL',columns='ENTREZID');mm<-unique(mm[!is.na(mm$ENTREZID),]);mm$ENTREZID<-as.character(mm$ENTREZID);maps<-split(mm$ENTREZID,mm$ENSEMBL)
g$Entrez<-vapply(g$gene_id,function(x)paste(sort(unique(maps[[x]])),collapse=';'),'');ng<-vapply(g$gene_id,function(x)length(unique(maps[[x]])),1L)
g$mapping_status<-ifelse(ng==0,'no_Entrez_for_Ensembl',ifelse(ng>1,'multiple_Entrez','mapped_one_to_one'))
for(con in unique(g$contrast)){
 ix<-which(g$contrast==con & ng==1);dup<-names(which(table(g$Entrez[ix])>1));g$mapping_status[ix[g$Entrez[ix] %in% dup]]<-'multiple_RNA_genes_for_Entrez'
}
g$GO_mapping_eligible<-g$mapping_status=='mapped_one_to_one';wt(g,'92_GO_mapping.tsv')
go<-AnnotationDbi::select(org.Mm.eg.db,keys=unique(g$Entrez[g$GO_mapping_eligible]),keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL'))
go<-unique(go[!is.na(go$GOALL)&go$ONTOLOGYALL=='CC',c('ENTREZID','GOALL')]);ont<-AnnotationDbi::Ontology(GO.db::GOTERM);go<-go[go$GOALL %in% names(ont),];stopifnot(all(ont[go$GOALL]=='CC'));names(go)<-c('Entrez','GO_ID');wt(go,'92_GO_annotation.tsv')
all<-list();sums<-list();objs<-list()
for(con in c('Bic_Nimo_vs_Bic','Bic_Nimo_vs_TTX')){
 d<-g[g$contrast==con,];b<-d[d$baseline_LAD & d$GO_mapping_eligible,];B<-b$Entrez;stopifnot(!anyDuplicated(B))
 for(direction in c('up','down')){
  sel<-d$baseline_LAD&d$DEG&d$classification==direction;F<-d$Entrez[sel&d$GO_mapping_eligible];key<-paste(con,direction,sep='__')
  er<-enrichGO(gene=F,universe=B,OrgDb=org.Mm.eg.db,keyType='ENTREZID',ont='CC',pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.05,minGSSize=10,maxGSSize=500)
  z<-if(is.null(er))data.frame() else er@result;sig<-if(is.null(er))data.frame() else as.data.frame(er);sp<-if(nrow(sig)>1)simplify(er,cutoff=.7,by='p.adjust',select_fun=min) else er;simp<-if(is.null(sp))data.frame() else as.data.frame(sp)
  if(nrow(z)){
   z$contrast<-con;z$direction<-direction;z$passes_original_cutoffs<-z$ID %in% sig$ID;z$retained_after_simplify<-z$ID %in% simp$ID
   z$hit_gene_ids<-vapply(strsplit(z$geneID,'/'),function(x)paste(sort(b$gene_id[match(x,b$Entrez)]),collapse=';'),'');z$hit_genes<-vapply(strsplit(z$geneID,'/'),function(x)paste(sort(b$gene_name[match(x,b$Entrez)]),collapse=';'),'');all[[key]]<-z
  }
  objs[[key]]<-list(full=er,simplified=sp)
  sums[[key]]<-data.frame(contrast=con,direction=direction,foreground_RNA=sum(sel),foreground_mapped=length(F),foreground_CC=length(intersect(F,go$Entrez)),background_RNA=sum(d$baseline_LAD),background_mapped=length(B),background_CC=length(intersect(B,go$Entrez)),tested_terms=nrow(z),significant=nrow(sig),simplified=nrow(simp),top_BH=if(nrow(z))min(z$p.adjust) else NA,unmapped_foreground=paste(d$gene_name[sel&!d$GO_mapping_eligible],collapse=';'))
 }
}
wt(do.call(rbind,all),'92_GO_full.tsv');wt(do.call(rbind,sums),'92_GO_summary.tsv');saveRDS(objs,file.path(D,'92_GO_objects.rds'));writeLines(capture.output(sessionInfo()),file.path(D,'92_GO_sessionInfo.txt'));print(do.call(rbind,sums))
