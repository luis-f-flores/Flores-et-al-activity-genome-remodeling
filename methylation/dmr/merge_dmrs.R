#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


library(GenomicRanges)
library(data.table)


args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L) stop("Usage: merge_dmrs.R filtered_DMRs.tsv comparison output_directory")
comparison <- args[2]
base_dir <- args[3]
dir.create(base_dir,recursive=TRUE,showWarnings=FALSE)
dmrs <- fread(args[1])
stopifnot(all(c("chr","start","end","pvalue","qvalue","meth.diff") %in% names(dmrs)))
hyper_dmrs <- dmrs[meth.diff>0]
hypo_dmrs <- dmrs[meth.diff<0]

hyper_gr <- GRanges(
  seqnames = hyper_dmrs$chr,
  ranges   = IRanges(start = hyper_dmrs$start, end = hyper_dmrs$end),
  strand   = "*",
  pvalue   = hyper_dmrs$pvalue,
  qvalue   = hyper_dmrs$qvalue,
  meth_diff = hyper_dmrs$meth.diff
)


hypo_gr <- GRanges(
  seqnames = hypo_dmrs$chr,
  ranges   = IRanges(start = hypo_dmrs$start, end = hypo_dmrs$end),
  strand   = "*",
  pvalue   = hypo_dmrs$pvalue,
  qvalue   = hypo_dmrs$qvalue,
  meth_diff = hypo_dmrs$meth.diff
)


hyper_merged <- reduce(hyper_gr)
hyper_overlaps <- findOverlaps(hyper_gr, hyper_merged)

hypo_merged <- reduce(hypo_gr)
hypo_overlaps <- findOverlaps(hypo_gr, hypo_merged)


merged_hyper <- data.frame(
  chr      = seqnames(hyper_merged),
  start    = start(hyper_merged),
  end      = end(hyper_merged),
  meth_diff = tapply(hyper_dmrs$meth.diff[queryHits(hyper_overlaps)], subjectHits(hyper_overlaps), mean),
  pvalue   = tapply(hyper_dmrs$pvalue[queryHits(hyper_overlaps)], subjectHits(hyper_overlaps), min),
  qvalue   = tapply(hyper_dmrs$qvalue[queryHits(hyper_overlaps)], subjectHits(hyper_overlaps), min)
)

merged_hypo <- data.frame(
  chr      = seqnames(hypo_merged),
  start    = start(hypo_merged),
  end      = end(hypo_merged),
  meth_diff = tapply(hypo_dmrs$meth.diff[queryHits(hypo_overlaps)], subjectHits(hypo_overlaps), mean),
  pvalue   = tapply(hypo_dmrs$pvalue[queryHits(hypo_overlaps)], subjectHits(hypo_overlaps), min),
  qvalue   = tapply(hypo_dmrs$qvalue[queryHits(hypo_overlaps)], subjectHits(hypo_overlaps), min)
)


write.table(merged_hyper,
            file = file.path(base_dir, paste0(comparison, "_Merged_DMRs_cov4_hyper.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(merged_hypo,
            file = file.path(base_dir, paste0(comparison, "_Merged_DMRs_cov4_hypo.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)


cat(paste0("\nMerged DMRs: ", comparison, " → Merged DMRs saved in:\n", base_dir, "\n"))
