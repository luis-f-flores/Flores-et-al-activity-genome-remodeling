#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}
args <- commandArgs(trailingOnly=TRUE)
root <- args[[1]]
suppressPackageStartupMessages({library(clusterProfiler);library(org.Mm.eg.db);library(AnnotationDbi)})
stopifnot(as.character(packageVersion('clusterProfiler'))=='4.6.2',as.character(packageVersion('org.Mm.eg.db'))=='3.16.0',as.character(packageVersion('GO.db'))=='3.16.0')
run <- file.path(root,'runs/2026-10-14_RNA_DMR_main_supp')
out <- file.path(run,'68_GO_refresh');dir.create(out,showWarnings=FALSE)
mappings <- list()
for(tp in c('6hr','24hr')) {
  before <- readLines(file.path(out,paste0(tp,'_Up_original.txt')))
  after <- readLines(file.path(out,paste0(tp,'_Up_corrected.txt')))
  ids <- mapIds(org.Mm.eg.db,keys=before,column='ENTREZID',keytype='SYMBOL',multiVals='first')
  mappings[[tp]] <- data.frame(timepoint=tp,symbol=names(ids),ENTREZID=as.character(ids),retained=names(ids)%in%after)
  e <- enrichGO(gene=unname(na.omit(ids[after])),OrgDb=org.Mm.eg.db,keyType='ENTREZID',ont='BP',pAdjustMethod='BH',pvalueCutoff=.05,qvalueCutoff=.05)
  saveRDS(e,file.path(out,paste0(tp,'_BP_corrected_unsimplified.rds')))
  write.csv(e@result,file.path(out,paste0(tp,'_BP_corrected_all_tested.csv')),row.names=FALSE)
  if(nrow(as.data.frame(e))>0)e <- simplify(e,cutoff=.7,by='p.adjust',select_fun=min)
  saveRDS(e,file.path(out,paste0(tp,'_BP_corrected_simplified.rds')))
  write.csv(as.data.frame(e),file.path(out,paste0('GO_BP_',tp,'_Up_corrected.csv')),row.names=FALSE)
}
write.table(do.call(rbind,mappings),file.path(out,'symbol_Entrez_mapping.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
