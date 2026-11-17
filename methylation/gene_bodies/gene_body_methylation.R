#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({library(methylKit);library(GenomicRanges);library(rtracklayer)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5L)stop("Usage: gene_body_methylation.R timepoint samples.tsv gene_layers_dir output_directory annotation_helpers.R")
tp<-args[1];stopifnot(tp%in%c('1hr','6hr','24hr'))
source(args[5]);samples<-read.delim(args[2],stringsAsFactors=FALSE);stopifnot(!anyDuplicated(samples$sample))
ids<-c(paste0(LETTERS[1:3],sub('hr$','',tp),'M'),paste0(LETTERS[1:3],sub('hr$','',tp),'K'))
selected<-samples[match(ids,samples$sample),];stopifnot(!anyNA(selected$sample))
raw_cpg_report_files<-setNames(list(setNames(selected$report,ids)),tp)
gene_body_bed_file<-file.path(args[3],'genes.bed');tss_bed_file<-file.path(args[3],'tss.bed')
OUTPUT_TABLES_DIR<-args[4];dir.create(OUTPUT_TABLES_DIR,recursive=TRUE,showWarnings=FALSE)
DMR_DELTA<-5;DMR_QVAL<-.05

cat("Running methylKit region-level gene-body test for timepoint:", tp, "\n")

raw_cpg_files <- raw_cpg_report_files[[tp]]
stopifnot(!is.null(raw_cpg_files), length(raw_cpg_files) == 6)

sample_ids <- as.list(names(raw_cpg_files))
treatment  <- ifelse(grepl("K$", names(raw_cpg_files)), 1, 0)

cat("Reading", length(raw_cpg_files), "Bismark CpG reports...\n")
myobj <- methRead(as.list(raw_cpg_files),
                   sample.id = sample_ids,
                   assembly = "mm10",
                   treatment = treatment,
                   pipeline = "bismarkCytosineReport",
                   mincov = 4)

cat("Filtering high coverage (hi.perc=99.9) and normalizing...\n")
myobj <- filterByCoverage(myobj, hi.perc = 99.9)
myobj <- normalizeCoverage(myobj)

gene_models <- load_collapsed_gene_models(gene_body_bed_file, tss_bed_file)
gene_body_windows <- get_gene_body_windows(gene_models)

cat("Aggregating per-CpG counts into per-gene-body regions (", length(gene_body_windows), "genes)...\n")
region_obj <- regionCounts(myobj, regions = gene_body_windows)

cat("Uniting samples (min.per.group=2)...\n")
meth_regions <- unite(region_obj, destrand = FALSE, min.per.group = 2L)

cat("Running calculateDiffMeth (one test per gene body)...\n")
region_diff <- calculateDiffMeth(meth_regions, test = "Chisq", adjust = "SLIM",
                                  effect = "wmean", mc.cores = 4, save.db = FALSE)

diff_df <- getData(region_diff)
diff_df$gene_symbol <- gene_body_windows$gene_symbol[GenomicRanges::findOverlaps(
  GenomicRanges::GRanges(diff_df$chr, IRanges::IRanges(diff_df$start, diff_df$end)),
  gene_body_windows, type = "equal", select = "first"
)]
diff_df$gene_id <- diff_df$gene_symbol
diff_df$timepoint <- tp

write.csv(diff_df,
          file.path(OUTPUT_TABLES_DIR, sprintf("08_genebody_methylkit_direct_%s_all.csv", tp)),
          row.names = FALSE)

genebody_significant <- diff_df[abs(diff_df$meth.diff) >= DMR_DELTA & diff_df$qvalue <= DMR_QVAL, ]

write.csv(genebody_significant,
          file.path(OUTPUT_TABLES_DIR, sprintf("08_genebody_significant_%s.csv", tp)),
          row.names = FALSE)

cat("Gene-body SIGNIFICANT count (DMR-style q<=", DMR_QVAL, "/|delta|>=", DMR_DELTA, "%,", tp, "):",
    nrow(genebody_significant), "\n")

