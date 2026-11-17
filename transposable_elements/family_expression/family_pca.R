#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}
suppressPackageStartupMessages({
  library(DESeq2); library(ggplot2); library(dplyr); library(readr)
})


outdir <- "SalmonTE_all_timepoints"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

expr_paths <- c(
  "SalmonTE_out_1h/EXPR.csv",
  "SalmonTE_out_6h/EXPR.csv",
  "SalmonTE_out_24h/EXPR.csv"
)


read_expr <- function(p) {
  x <- read.csv(p, check.names = FALSE)
  stopifnot("TE" %in% names(x))
  rownames(x) <- x$TE
  x[, setdiff(colnames(x), "TE"), drop = FALSE]
}
lst <- lapply(expr_paths, read_expr)
common_tes <- Reduce(intersect, lapply(lst, rownames))
lst <- lapply(lst, function(x) x[common_tes, , drop = FALSE])
counts <- do.call(cbind, lst)


samples <- colnames(counts)

parse_one <- function(s) {

  m <- regexec("^([A-Z])(\\d+)([MK])$", s)
  r <- regmatches(s, m)[[1]]
  if (length(r) != 4) stop(paste("Unrecognized sample name:", s))
  data.frame(sample=s,
             rep=r[2],
             time=as.integer(r[3]),
             condition=ifelse(r[4]=="M","control","KCl"),
             stringsAsFactors=FALSE)
}
meta <- do.call(rbind, lapply(samples, parse_one))
row.names(meta) <- meta$sample
meta$time <- factor(meta$time, levels=sort(unique(meta$time)))
meta$condition <- factor(meta$condition, levels=c("control","KCl"))


dds <- DESeqDataSetFromMatrix(round(counts), colData=meta, design=~1)
dds <- estimateSizeFactors(dds)
vst <- varianceStabilizingTransformation(dds, blind=TRUE)


pdat <- DESeq2::plotPCA(vst, intgroup=c("condition","time"), returnData=TRUE)
percentVar <- round(100 * attr(pdat, "percentVar"), 1)

p <- ggplot(pdat, aes(PC1, PC2, color=condition, shape=time)) +
  geom_point(size=3, alpha=0.9) +
  labs(title="TE families PCA (1h + 6h + 24h)",
       x=paste0("PC1 (", percentVar[1], "%)"),
       y=paste0("PC2 (", percentVar[2], "%)")) +
  coord_equal() +
  theme_minimal(base_size = 12)


ggsave(file.path(outdir, "TE_PCA_all_timepoints.pdf"), p, width=7, height=5)
write_csv(pdat, file.path(outdir, "TE_PCA_all_timepoints_scores.csv"))
write_csv(as.data.frame(colData(vst)), file.path(outdir, "TE_sample_metadata_used.csv"))


cat("\nWrote:\n",
    file.path(outdir, "TE_PCA_all_timepoints.pdf"), "\n",
    file.path(outdir, "TE_PCA_all_timepoints_scores.csv"), "\n",
    file.path(outdir, "TE_sample_metadata_used.csv"), "\n")
