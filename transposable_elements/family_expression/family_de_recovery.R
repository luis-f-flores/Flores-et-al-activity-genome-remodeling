#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(pheatmap); library(matrixStats)
})
expr <- read.csv("../SalmonTE_out_24h/EXPR.csv", check.names=FALSE)
rownames(expr) <- expr$TE
expr <- expr[, setdiff(colnames(expr), "TE"), drop=FALSE]
expr <- expr[, c("A24M","B24M","C24M","A24K","B24K","C24K")]
samples <- data.frame(
  row.names=c("A24M","B24M","C24M","A24K","B24K","C24K"),
  condition=factor(c("control","control","control","KCl","KCl","KCl"),
                   levels=c("control","KCl"))
)
dds <- DESeqDataSetFromMatrix(round(expr), colData=samples, design=~condition)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef="condition_KCl_vs_control", type="apeglm")
outdir <- "../SalmonTE_DE_24h"; dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
write.csv(as.data.frame(res), file=file.path(outdir, "TE_DE_results.csv"))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
pca <- plotPCA(vsd, intgroup="condition") + ggtitle("TE families PCA")
ggsave(file.path(outdir, "PCA.pdf"), pca, width=6, height=5)
vst_mat <- assay(vsd); rv <- matrixStats::rowVars(vst_mat)
sel <- order(rv, decreasing=TRUE)[seq_len(min(50, length(rv)))]
pdf(file.path(outdir, "heatmap.pdf"), width=6, height=8)
pheatmap(vst_mat[sel, ], show_rownames=FALSE,
         annotation_col=as.data.frame(colData(vsd)))
dev.off()
