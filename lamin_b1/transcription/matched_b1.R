args <- commandArgs(trailingOnly=TRUE)
input <- args[[1]]
run <- args[[2]]
dir.create(run,recursive=TRUE,showWarnings=FALSE)
d <- read.delim(file.path(input,'fixed_observations.tsv'),check.names=FALSE)
rows <- lapply(c('down','up'),function(direction) {
  a <- d[d$RNA_direction==direction,]
  a$matched_set <- factor(a$matched_set)
  fit <- lm(B1_delta ~ 0 + matched_set + is_DEG,data=a,weights=analysis_weight)
  clusters <- sort(unique(a$LAD_cluster))
  g <- length(clusters)
  beta <- unname(coef(fit)['is_DEG'])
  jack <- vapply(clusters,function(k) {
    sub <- droplevels(a[a$LAD_cluster!=k,])
    f <- lm(B1_delta ~ 0 + matched_set + is_DEG,data=sub,weights=analysis_weight)
    unname(coef(f)['is_DEG'])
  },numeric(1))
  stopifnot(all(is.finite(jack)))
  se <- sqrt((g-1)/g*sum((jack-mean(jack))^2))
  write.table(data.frame(LAD_cluster=clusters,estimate_without_LAD=jack),
              file.path(run,paste0('R_jackknife_',direction,'.tsv')),
              sep='\t',row.names=FALSE,quote=FALSE)
  t <- beta/se
  df <- g-1
  critical <- qt(.975,df)
  data.frame(RNA_direction=direction,n_DEGs=nlevels(a$matched_set),
             observations=nrow(a),unique_genes=length(unique(a$gene_id)),
             LAD_clusters=g,mean_excess_log2=beta,SE_LAD_cluster=se,t_statistic=t,df=df,
             CI95_low=beta-critical*se,CI95_high=beta+critical*se,
             P_LAD_cluster=2*pt(-abs(t),df),model_rank=fit$rank,
             jackknife_mean=mean(jack),jackknife_min=min(jack),jackknife_max=max(jack))
})
out <- do.call(rbind,rows)
out$Holm_P_LAD_cluster <- p.adjust(out$P_LAD_cluster,method='holm')
write.table(out,file.path(run,'R_LAD_adjusted.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(out,digits=15)
sessionInfo()
