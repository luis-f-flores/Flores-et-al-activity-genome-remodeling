#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


options(stringsAsFactors=FALSE,warn=1);Sys.setenv(LANGUAGE='en',TZ='UTC');set.seed(47)
suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db);library(GO.db);library(AnnotationDbi);library(jsonlite);library(digest)})
stopifnot(as.character(packageVersion('clusterProfiler'))=='4.6.2',as.character(packageVersion('org.Mm.eg.db'))=='3.16.0',as.character(packageVersion('GO.db'))=='3.16.0')
args<-commandArgs(TRUE);IN<-args[1];OUT<-args[2];SOURCE<-args[3]
stopifnot(!file.exists(file.path(OUT,'tables/64_GO_full.tsv')))
for(d in c('tables','logs','qc'))dir.create(file.path(OUT,d),recursive=TRUE,showWarnings=FALSE)
wt<-function(x,n)write.table(x,file.path(OUT,n),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
g<-read.delim(IN,check.names=FALSE,quote='');stopifnot(all(g$scheme=='condition_depth'),!anyDuplicated(paste(g$contrast,g$gene_id)))
u<-unique(g[,c('gene_id','gene_name','gene_biotype')]);stopifnot(!anyDuplicated(u$gene_id))
mm<-AnnotationDbi::select(org.Mm.eg.db,keys=u$gene_id,keytype='ENSEMBL',columns='ENTREZID');mm<-unique(mm[!is.na(mm$ENTREZID),]);mm$ENTREZID<-as.character(mm$ENTREZID)
maps<-split(mm$ENTREZID,mm$ENSEMBL);ng<-vapply(u$gene_id,function(x)length(unique(maps[[x]])),integer(1));u$Entrez<-vapply(u$gene_id,function(x)paste(sort(unique(maps[[x]])),collapse=';'),character(1))
single<-mm[mm$ENSEMBL %in% u$gene_id[ng==1],];reverse<-names(which(table(single$ENTREZID)>1))
u$mapping_status<-ifelse(ng==0,'no_Entrez_for_Ensembl',ifelse(ng>1,'multiple_Entrez',ifelse(u$Entrez %in% reverse,'multiple_genes_for_Entrez','mapped_one_to_one')))
g$Entrez<-u$Entrez[match(g$gene_id,u$gene_id)];g$mapping_status<-u$mapping_status[match(g$gene_id,u$gene_id)];g$GO_mapping_eligible<-g$mapping_status=='mapped_one_to_one'
wt(u,'tables/64_GO_ID_mapping.tsv');wt(g,'tables/64_GO_gene_sets.tsv');valid<-u[u$mapping_status=='mapped_one_to_one',]
go<-AnnotationDbi::select(org.Mm.eg.db,keys=valid$Entrez,keytype='ENTREZID',columns=c('GOALL','ONTOLOGYALL'));go<-unique(go[!is.na(go$GOALL)&go$ONTOLOGYALL %in% c('BP','CC','MF'),c('ENTREZID','GOALL','ONTOLOGYALL')])
ontologies<-AnnotationDbi::Ontology(GO.db::GOTERM);go<-go[go$GOALL %in% names(ontologies),];stopifnot(all(go$ONTOLOGYALL==unname(ontologies[go$GOALL])))
names(go)<-c('Entrez','GO_ID','ontology');go$gene_id<-valid$gene_id[match(go$Entrez,valid$Entrez)];wt(go,'tables/64_gene_GO.tsv')
wt(rbind(data.frame(database='org.Mm.eg.db',metadata(org.Mm.eg.db)),data.frame(database='GO.db',metadata(GO.db))),'tables/64_GO_annotation_metadata.tsv')
allresults<-list();summaries<-list();backgrounds<-list();objects<-list();simplified<-list()
invisible(clusterProfiler:::get_GO_data(org.Mm.eg.db,'ALL','ENTREZID'))
for(con in c('Bic_vs_TTX','Bic_Nimo_vs_Bic','Bic_Nimo_vs_TTX'))for(direction in c('Lower B1','Higher B1')){
 sel<-g$contrast==con;fg<-sel & g$B1_group==direction;B<-g$Entrez[sel & g$GO_mapping_eligible];F<-g$Entrez[fg & g$GO_mapping_eligible]
 stopifnot(all(F %in% B),!anyDuplicated(B),!anyDuplicated(F),all(g$B1_q[fg]<.05))
 for(ont in c('BP','CC','MF')){
  key<-paste(con,direction,ont,sep='__');ann<-unique(go$Entrez[go$ontology==ont]);ba<-intersect(B,ann);fa<-intersect(F,ann)
  cat(key,': foreground',length(fa),'of',length(ba),'annotated background\n');flush.console()
  backgrounds[[key]]<-data.frame(contrast=con,direction=direction,ontology=ont,background_TSS=sum(sel),foreground_TSS=sum(fg),background_mapped=length(B),foreground_mapped=length(F),background_annotated=length(ba),foreground_annotated=length(fa),unmapped_foreground=paste(g$gene_name[fg & !g$GO_mapping_eligible],collapse=';'))
  er<-if(length(fa))enrichGO(gene=F,universe=B,OrgDb=org.Mm.eg.db,keyType='ENTREZID',ont=ont,pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.05,minGSSize=10,maxGSSize=500) else NULL
  if(is.null(er)||!nrow(er@result)){
   summaries[[key]]<-data.frame(contrast=con,direction=direction,ontology=ont,N=length(ba),n=length(fa),tested_terms=0,BH_005_terms=0,original_cutoff_terms=0,simplified_terms=0,top_term='',top_BH=NA_real_,top_k=NA_integer_,top_genes='',status=if(length(fa))'no_testable_terms' else 'no_annotated_foreground');next
  }
  objects[[key]]<-er;stopifnot(setequal(er@universe,ba));z<-er@result;z$contrast<-con;z$direction<-direction;z$ontology<-ont
  ratio<-function(x,i)as.integer(vapply(strsplit(x,'/',fixed=TRUE),'[',character(1),i))
  z$k<-ratio(z$GeneRatio,1);z$n<-ratio(z$GeneRatio,2);z$K<-ratio(z$BgRatio,1);z$N<-ratio(z$BgRatio,2);z$expected<-z$n*z$K/z$N;z$fold_enrichment<-z$k/z$expected
  stopifnot(all(z$n==length(fa)),all(z$N==length(ba)))
  z$hit_genes<-vapply(strsplit(z$geneID,'/',fixed=TRUE),function(ids)paste(sort(valid$gene_name[match(ids,valid$Entrez)]),collapse=';'),character(1))
  z$hit_gene_ids<-vapply(strsplit(z$geneID,'/',fixed=TRUE),function(ids)paste(sort(valid$gene_id[match(ids,valid$Entrez)]),collapse=';'),character(1))
  sig<-as.data.frame(er);z$passes_original_cutoffs<-z$ID %in% sig$ID;sp<-if(nrow(sig)>1)simplify(er,cutoff=.7,by='p.adjust',select_fun=min)else er;simplified[[key]]<-sp;spdf<-as.data.frame(sp);z$retained_after_simplify<-z$ID %in% spdf$ID;allresults[[key]]<-z
  best<-z[order(z$p.adjust,z$pvalue,z$ID)[1],];summaries[[key]]<-data.frame(contrast=con,direction=direction,ontology=ont,N=length(ba),n=length(fa),tested_terms=nrow(z),BH_005_terms=sum(z$p.adjust<=.05),original_cutoff_terms=nrow(sig),simplified_terms=nrow(spdf),top_term=best$Description,top_BH=best$p.adjust,top_k=best$k,top_genes=best$hit_genes,status='tested')
 }
}
full<-do.call(rbind,allresults);stopifnot(nrow(full)>0);rownames(full)<-NULL
full<-full[order(full$contrast,full$direction,full$ontology,full$p.adjust,full$pvalue,full$ID),]
wt(full,'tables/64_GO_full.tsv');wt(full[full$passes_original_cutoffs,],'tables/64_GO_significant.tsv');wt(full[full$retained_after_simplify,],'tables/64_GO_simplified.tsv');wt(do.call(rbind,summaries),'tables/64_GO_summary.tsv');wt(do.call(rbind,backgrounds),'tables/64_GO_backgrounds.tsv')
saveRDS(list(unsimplified=objects,simplified=simplified),file.path(OUT,'tables/64_GO_objects.rds'));writeLines(capture.output(sessionInfo()),file.path(OUT,'logs/64_GO_sessionInfo.txt'))
stamp<-function(p)list(path=p,bytes=unname(file.info(p)$size),sha256=digest(file=p,algo='sha256'))
packages<-c('clusterProfiler','org.Mm.eg.db','GO.db','AnnotationDbi','DOSE','GOSemSim','qvalue')
receipt<-list(status='computed',R=R.version.string,packages=setNames(lapply(packages,function(x)as.character(packageVersion(x))),packages),inputs=lapply(c(IN,SOURCE,AnnotationDbi::dbfile(org.Mm.eg.db),AnnotationDbi::dbfile(GO.db)),stamp),settings=list(minGSSize=10,maxGSSize=500,pvalueCutoff=.05,qvalueCutoff=.05,simplify=.7,correction='BH per direction/contrast/ontology',background='All uniquely mapped assessed RefGene-linked TSS genes in each pairwise LAD universe; no RNA test/DEG filter; one-to-one Entrez mapping'))
write_json(receipt,file.path(OUT,'qc/64_GO_QC.json'),pretty=TRUE,auto_unbox=TRUE,na='null');print(do.call(rbind,summaries),row.names=FALSE)
