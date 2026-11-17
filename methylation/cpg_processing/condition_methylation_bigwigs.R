#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({library(methylKit);library(GenomicRanges);library(Biostrings);library(rtracklayer)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop("Usage: condition_methylation_bigwigs.R timepoint samples.tsv mm10.fa output_directory")
tp<-args[1];stopifnot(tp%in%c('1hr','6hr','24hr'))
samples<-read.delim(args[2],stringsAsFactors=FALSE);stopifnot(!anyDuplicated(samples$sample))
ids<-c(paste0(LETTERS[1:3],sub('hr$','',tp),'M'),paste0(LETTERS[1:3],sub('hr$','',tp),'K'))
selected<-samples[match(ids,samples$sample),];stopifnot(!anyNA(selected$sample))
raw_cpg_report_files<-setNames(list(setNames(selected$report,ids)),tp)
MM10_FASTA<-args[3];OUTPUT_FIGURES_DIR<-args[4];dir.create(OUTPUT_FIGURES_DIR,recursive=TRUE,showWarnings=FALSE)

cat("Building genome-wide mean-%methylation bigwigs for timepoint:", tp, "\n")

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

chrom_lengths <- Biostrings::fasta.seqlengths(MM10_FASTA)

make_condition_bw <- function(idx, label) {
  cat("Uniting", length(idx), label, "replicates genome-wide (per-CpG, no region aggregation)...\n")
  sub_obj <- methylKit::reorganize(myobj, sample.ids = sample_ids[idx], treatment = rep(0, length(idx)))
  united <- unite(sub_obj, destrand = FALSE, min.per.group = 1L)

  pm <- percMethylation(united)
  mean_pct <- rowMeans(pm, na.rm = TRUE)

  d <- getData(united)
  gr <- GenomicRanges::GRanges(d$chr, IRanges::IRanges(d$start, d$end), score = mean_pct)
  gr <- gr[stats::complete.cases(data.frame(score = mean_pct)), ]

  seqlengths(gr) <- chrom_lengths[as.character(unique(GenomicRanges::seqnames(gr)))]

  out_path <- file.path(OUTPUT_FIGURES_DIR, sprintf("16_genomewide_%s_%s_mean_pct_meth.bw", label, tp))
  rtracklayer::export.bw(gr, out_path)
  cat(label, "bigwig written:", out_path, "(", length(gr), "CpGs )\n")
}

k_idx <- which(treatment == 1)
m_idx <- which(treatment == 0)

make_condition_bw(k_idx, "KCl")
make_condition_bw(m_idx, "Mannitol")

cat("\nCondition-mean methylation bigWigs complete for", tp, ". Output in", OUTPUT_FIGURES_DIR, "\n")

