#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(DESeq2)
  library(jsonlite)
})

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: salmonte_cross_paradigm.R KCl_root Bic_quant_root output_directory")
root <- args[1]
bic_quant <- args[2]
out <- args[3]
dir.create(out, recursive=TRUE, showWarnings=FALSE)

read_quant <- function(dir) {
  qfile <- file.path(dir, "quant.sf")
  mfile <- file.path(dir, "aux_info/meta_info.json")
  stopifnot(file.exists(qfile), file.exists(mfile))
  q <- read.delim(qfile, check.names=FALSE)
  m <- jsonlite::fromJSON(mfile)
  stopifnot(identical(q$Name, unique(q$Name)), nrow(q)==301L,
            is.finite(m$num_processed), m$num_processed>0,
            all(is.finite(q$NumReads)), all(q$NumReads>=0))
  list(counts=setNames(q$NumReads, q$Name), processed=m$num_processed,
       target_hash=m$index_name_hash, sequence_hash=m$index_seq_hash,
       qfile=qfile, mfile=mfile)
}

old_ids <- unlist(lapply(c("1h","6h","24h"), function(t)
  as.vector(outer(LETTERS[1:3], c("K","M"), function(a,k)
    paste0(a, sub("h$","",t), k)))))
new_ids <- c(paste0("Basal_", LETTERS[1:3]),
             paste0("Bic_", LETTERS[1:3]),
             paste0("Bic_nimo_", LETTERS[1:3]))
old_time <- sub("^[ABC](1|6|24)[KM]$", "\\1h", old_ids)
stopifnot(all(old_time %in% c("1h","6h","24h")))
old <- setNames(lapply(seq_along(old_ids), function(i)
  read_quant(file.path(root, paste0("SalmonTE_out_",old_time[i]), old_ids[i]))), old_ids)
new <- setNames(lapply(new_ids, function(s)
  read_quant(file.path(bic_quant,s))), new_ids)
allq <- c(old,new)
ref_names <- names(allq[[1]]$counts)
stopifnot(all(vapply(allq, function(x) setequal(names(x$counts),ref_names), logical(1))),
          length(unique(vapply(allq, `[[`, "", "target_hash")))==1L,
          length(unique(vapply(allq, `[[`, "", "sequence_hash")))==1L)
make_matrix <- function(qs) {
  m <- vapply(qs, function(x) x$counts[ref_names], numeric(length(ref_names)))
  rownames(m) <- ref_names
  m
}
counts <- make_matrix(allq)
depth <- vapply(allq, `[[`, numeric(1), "processed")
cpm <- sweep(counts, 2, depth[colnames(counts)]/1e6, "/")
effect <- function(treatment, control) {
  log2((rowMeans(cpm[,treatment,drop=FALSE])+1)/
       (rowMeans(cpm[,control,drop=FALSE])+1))
}

contrasts <- list()
for (t in c("1h","6h","24h")) {
  treatment <- paste0(LETTERS[1:3], sub("h$","",t), "K")
  control <- paste0(LETTERS[1:3], sub("h$","",t), "M")
  contrasts[[paste0("KCl_Mann_",t)]] <- effect(treatment,control)
}
contrasts$Bic_TTX <- effect(paste0("Bic_",LETTERS[1:3]),paste0("Basal_",LETTERS[1:3]))
contrasts$Nimo_TTX <- effect(paste0("Bic_nimo_",LETTERS[1:3]),paste0("Basal_",LETTERS[1:3]))
contrasts$Nimo_Bic <- effect(paste0("Bic_nimo_",LETTERS[1:3]),paste0("Bic_",LETTERS[1:3]))
new_filter <- rowSums(counts[,new_ids[1:6],drop=FALSE]>=10)>=3
family <- data.frame(Name=ref_names, do.call(cbind,contrasts),
                     new_measured=new_filter,
                     mean_CPM_TTX=rowMeans(cpm[,new_ids[1:3],drop=FALSE]),
                     mean_CPM_Bic=rowMeans(cpm[,new_ids[4:6],drop=FALSE]),
                     mean_CPM_Nimo=rowMeans(cpm[,new_ids[7:9],drop=FALSE]),
                     check.names=FALSE)
kcl_de <- setNames(lapply(c("1h","6h","24h"), function(t) {
  p <- file.path(root,paste0("SalmonTE_DE_",t),"TE_DE_results.csv")
  d <- read.csv(p,row.names=1,check.names=FALSE)
  stopifnot(setequal(rownames(d),ref_names),
            all(c("log2FoldChange","padj") %in% names(d)))
  d
}),c("1h","6h","24h"))
for (t in names(kcl_de)) {
  family[[paste0("KCl_",t,"_DE_lfc")]] <- kcl_de[[t]][ref_names,"log2FoldChange"]
  family[[paste0("KCl_",t,"_DE_padj")]] <- kcl_de[[t]][ref_names,"padj"]
  ids <- old_ids[old_time==t]
  family[[paste0("KCl_",t,"_measured")]] <- rowSums(counts[,ids,drop=FALSE]>=10)>=3
}
write.table(family,file.path(out,"CROSS_PARADIGM_FAMILY_EFFECTS.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE,na="NA")

concord <- lapply(c("1h","6h","24h"), function(t) {
  ids <- old_ids[old_time==t]
  eligible <- new_filter & rowSums(counts[,ids,drop=FALSE]>=10)>=3
  x <- contrasts[[paste0("KCl_Mann_",t)]][eligible]
  y <- contrasts$Bic_TTX[eligible]
  z <- contrasts$Nimo_TTX[eligible]
  sig_up <- !is.na(kcl_de[[t]][ref_names,"padj"]) &
    kcl_de[[t]][ref_names,"padj"]<0.05 &
    kcl_de[[t]][ref_names,"log2FoldChange"]>0
  data.frame(time=t,n=sum(eligible),spearman_Bic=cor(x,y,method="spearman"),
             pearson_Bic=cor(x,y),
             spearman_Nimo=cor(x,z,method="spearman"),
             KCl_FDR_up=sum(sig_up),KCl_FDR_up_Bic_same_direction=sum(sig_up & contrasts$Bic_TTX>0))
})
concord <- do.call(rbind,concord)
write.table(concord,file.path(out,"CROSS_PARADIGM_CONCORDANCE.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)

iap <- c("IAPA_MM","IAPEZI","IAPEYI","IAPEY3_I","IAPEY4_I","IAPLTR1a_I_MM")
stopifnot(all(iap %in% rownames(cpm)))
iap_samples <- data.frame(Name=rep(iap,each=length(new_ids)),
                          sample_id=rep(new_ids,length(iap)),
                          CPM=as.vector(t(cpm[iap,new_ids,drop=FALSE])),
                          NumReads=as.vector(t(counts[iap,new_ids,drop=FALSE])))
write.table(iap_samples,file.path(out,"IAP_REPLICATE_VALUES.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)

model_ids <- new_ids
model_keep <- rowSums(counts[,model_ids,drop=FALSE]>=10)>=3
model_counts <- round(counts[model_keep,model_ids,drop=FALSE])
condition <- factor(c(rep("TTX",3),rep("Bic",3),rep("Bic_Nimo",3)),
                    levels=c("TTX","Bic","Bic_Nimo"))
coldata <- data.frame(row.names=model_ids,condition=condition)
dds <- DESeqDataSetFromMatrix(model_counts,coldata,design=~condition)
sf <- depth[model_ids]/exp(mean(log(depth[model_ids])))
sizeFactors(dds) <- sf
dds <- DESeq(dds,quiet=TRUE,minReplicatesForReplace=Inf)
model_contrasts <- list(Bic_TTX=c("condition","Bic","TTX"),
                        Nimo_TTX=c("condition","Bic_Nimo","TTX"),
                        Nimo_Bic=c("condition","Bic_Nimo","Bic"))
model <- do.call(rbind,lapply(names(model_contrasts),function(nm) {
  r <- results(dds,contrast=model_contrasts[[nm]],independentFiltering=FALSE)
  data.frame(Name=rownames(r),contrast=nm,baseMean=r$baseMean,
             log2FoldChange=r$log2FoldChange,lfcSE=r$lfcSE,
             pvalue=r$pvalue,padj=r$padj)
}))
write.table(model,file.path(out,"P2609_TE_DESEQ2_INPUT_DEPTH.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE,na="NA")


cat("Saved family effects, IAP replicate counts and DESeq2 contrasts.\n")
