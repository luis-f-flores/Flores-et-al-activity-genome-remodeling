#!/usr/bin/env Rscript
# Shared WGBS input names, thresholds and analysis configuration.
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(readxl)
})

PROJECT_ROOT <- gb_path("wgbs_package")

RUN_ID <- "2026-06-29_full_rerun"
RUN_DIR <- file.path(PROJECT_ROOT, "runs", RUN_ID)

OUTPUT_TABLES_DIR  <- file.path(RUN_DIR, "tables")
source(file.path(PROJECT_ROOT, "R", "utils_figure_paths.R"))
OUTPUT_FIGURES_DIR <- revision_figure_dir(RUN_ID)
dir.create(OUTPUT_TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

TIMEPOINTS <- c("1hr", "6hr", "24hr")

SNAPSHOT_ROOT <- PROJECT_ROOT
WGBS_ROOT     <- file.path(SNAPSHOT_ROOT, "wgbs")
RNA_ROOT      <- file.path(SNAPSHOT_ROOT, "rna_seq")
GENESET_ROOT  <- file.path(SNAPSHOT_ROOT, "genesets")

cpg_baseline_0hr_bw <- NULL

dmc_files <- list(
  "1hr"  = file.path(WGBS_ROOT, "DMC_1M_vs_1K_analysis_cov4/all_DMCs_1M_vs_1K.txt"),
  "6hr"  = file.path(WGBS_ROOT, "DMC_6M_vs_6K_analysis_cov4/all_DMCs_6M_vs_6K.txt"),
  "24hr" = file.path(WGBS_ROOT, "DMC_24M_vs_24hrKCl_analysis_cov4/all_DMCs_24M_vs_24hrKCl.txt")
)

dmc_baseline_files <- list(
  "1hr"  = file.path(WGBS_ROOT, "DMC_0hr_vs_1hK_analysis/all_DMCs_0hr_vs_1hK.txt"),
  "24hr" = file.path(WGBS_ROOT, "DMC_0hr_vs_24hrKCl_analysis/all_DMCs_0hr_vs_24hrKCl.txt")
)

dmr_files <- list(
  "1hr"  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp/DMRs_diff5_q005.txt"),
  "6hr"  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp/DMRs_diff5_q005.txt"),
  "24hr" = file.path(WGBS_ROOT, "methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp/DMRs_diff5_q005.txt")
)

dmr_merged_files <- list(
  "1hr"  = list(hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp/Merged_q05_meth5/1M_vs_1K_Merged_DMRs_cov4_hyper.txt"),
                hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp/Merged_q05_meth5/1M_vs_1K_Merged_DMRs_cov4_hypo.txt")),
  "6hr"  = list(hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp/Merged_q05_meth5/6M_vs_6K_Merged_DMRs_cov4_hyper.txt"),
                hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp/Merged_q05_meth5/6M_vs_6K_Merged_DMRs_cov4_hypo.txt")),
  "24hr" = list(hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp/Merged_q05_meth5/24M_vs_24K_Merged_DMRs_cov4_hyper.txt"),
                hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp/Merged_q05_meth5/24M_vs_24K_Merged_DMRs_cov4_hypo.txt"))
)

dmr_unfiltered_files <- list(
  "1hr"  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp/DMRs_all_unfiltered_win_500bp.txt"),
  "6hr"  = file.path(WGBS_ROOT, "methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp/DMRs_all_unfiltered_win_500bp.txt"),
  "24hr" = file.path(WGBS_ROOT, "methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp/DMRs_all_unfiltered_win_500bp.txt")
)

RAW_CPG_DIR <- file.path(WGBS_ROOT, "raw_cpg_reports")
raw_cpg_report_files <- list(
  "0hr" = c(
    A0 = file.path(RAW_CPG_DIR, "A0_CKDL240037499-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    B0 = file.path(RAW_CPG_DIR, "B0_CKDL240037506-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    C0 = file.path(RAW_CPG_DIR, "C0_CKDL240037513-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt")
  ),
  "1hr" = c(
    A1K = file.path(RAW_CPG_DIR, "A1K_CKDL240037501-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    A1M = file.path(RAW_CPG_DIR, "A1M_CKDL240037500-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    B1K = file.path(RAW_CPG_DIR, "B1K_CKDL240037508-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    B1M = file.path(RAW_CPG_DIR, "B1M_CKDL240037507-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    C1K = file.path(RAW_CPG_DIR, "C1K_CKDL240037515-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    C1M = file.path(RAW_CPG_DIR, "C1M_CKDL240037514-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt")
  ),
  "6hr" = c(
    A6K = file.path(RAW_CPG_DIR, "A6K_CKDL240037503-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    A6M = file.path(RAW_CPG_DIR, "A6M_CKDL240037502-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    B6K = file.path(RAW_CPG_DIR, "B6K_CKDL240037510-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    B6M = file.path(RAW_CPG_DIR, "B6M_CKDL240037509-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    C6K = file.path(RAW_CPG_DIR, "C6K_CKDL240037517-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    C6M = file.path(RAW_CPG_DIR, "C6M_CKDL240037516-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt")
  ),
  "24hr" = c(
    A24K = file.path(RAW_CPG_DIR, "A24K_CKDL240037505-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    A24M = file.path(RAW_CPG_DIR, "A24M_CKDL240037504-1A_22J3WNLT4_L3_1_val_1_bismark_hisat2.CpG_report.txt"),
    B24K = file.path(RAW_CPG_DIR, "B24K_CKDL240037512-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    B24M = file.path(RAW_CPG_DIR, "B24M_CKDL240037511-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    C24K = file.path(RAW_CPG_DIR, "C24K_CKDL240037519-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt"),
    C24M = file.path(RAW_CPG_DIR, "C24M_CKDL240037518-1A_22J3WNLT4_L4_1_val_1_bismark_hisat2.CpG_report.txt")
  )
)

dmg_promoter_files <- list(
  "1hr" = list(
    hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/1M_vs_1K_hyper_Promoter_DMRs_ForVolcano_meth25_q01.txt"),
    hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/1M_vs_1K_hypo_Promoter_DMRs_ForVolcano_meth25_q01.txt")
  ),
  "6hr" = list(
    hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/6M_vs_6K_hyper_Promoter_DMRs_ForVolcano_meth25_q01.txt"),
    hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/6M_vs_6K_hypo_Promoter_DMRs_ForVolcano_meth25_q01.txt")
  ),
  "24hr" = list(
    hyper = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/24M_vs_24K_hyper_Promoter_DMRs_ForVolcano_meth25_q01.txt"),
    hypo  = file.path(WGBS_ROOT, "methylKit_DMR/DMG_overlap/Promoter_DMGs/Results_with_5meth_q05/DMGs_q01_meth25_from_DMRs_q05meth5/24M_vs_24K_hypo_Promoter_DMRs_ForVolcano_meth25_q01.txt")
  )
)

coverage_summary_file <- file.path(WGBS_ROOT, "methylKit_QC/methylKit_analysis/coverage_summary.csv")

deg_files <- list(
  "1hr"  = file.path(RNA_ROOT, "Hour_1_KClvsHour_1_Mann_deg.xls"),
  "6hr"  = file.path(RNA_ROOT, "Hour_6_KClvsHour_6_Mann_deg.xls"),
  "24hr" = file.path(RNA_ROOT, "Hour_24_KClvsHour_24_Mann_deg.xls")
)

fpkm_matrix_file <- NULL

master_fpkm_file <- file.path(RNA_ROOT, "gene_fpkm.xls")

gene_body_bed_file <- file.path(RNA_ROOT, "gene_layers_mm10/genes.bed")
exon_bed_file       <- file.path(RNA_ROOT, "gene_layers_mm10/exons.bed")
intron_bed_file      <- file.path(RNA_ROOT, "gene_layers_mm10/introns.bed")
tss_bed_file         <- file.path(RNA_ROOT, "gene_layers_mm10/tss.bed")

te_context_promoter_bed_file <- file.path(RNA_ROOT, "gene_layers_mm10/promoters_-1kb_+100bp.bed")

gtf_file <- NULL

bdnf_exon_IV_region <- list(chr = "chr2", start = 109692501, end = 109693000)


activity_gene_lists <- list(
  dPRG = readLines(file.path(GENESET_ROOT, "dPRG_gene_list.txt")),
  SRG  = readLines(file.path(GENESET_ROOT, "SRG_gene_list.txt")),
  rPRG = readLines(file.path(GENESET_ROOT, "rPRG_gene_list.txt"))
)

.activity_alias_map <- c(Fam46a = "Tent5a", Gltscr2 = "Nop53")
activity_gene_lists <- lapply(activity_gene_lists, function(g) {
  g <- trimws(g)
  hits <- intersect(names(.activity_alias_map), g)
  unique(c(g, unname(.activity_alias_map[hits])))
})

te_catalog_bed_file  <- file.path(RNA_ROOT, "TE_catalog_mm10/TE_mm10.bed")
te_catalog_meta_file <- file.path(RNA_ROOT, "TE_catalog_mm10/TE_meta.tsv")

DMC_PVAL <- 0.001
DMC_DELTA <- 5

DMR_QVAL <- 0.05
DMR_DELTA <- 5

DMG_QVAL <- 0.01
DMG_DELTA <- 25

DEG_QVAL <- 0.01
DEG_FC   <- 1.3

FPKM_BACKGROUND_OLD <- 0
FPKM_BACKGROUND_NEW <- 5

PROMOTER_WINDOW <- c(-1000, 100)

.check_path <- function(name, value) {
  if (is.null(value)) return(NULL)
  if (is.character(value) && length(value) == 1 && !file.exists(value)) {
    return(sprintf("  [ FILE MISSING ] %s -> %s", name, value))
  }
  NULL
}

.validate_config <- function() {
  msgs <- character(0)

  msgs <- c(msgs, .check_path("cpg_baseline_0hr_bw", cpg_baseline_0hr_bw))
  for (tp in TIMEPOINTS) {
    msgs <- c(msgs, .check_path(sprintf("dmc_files[['%s']]", tp), dmc_files[[tp]]))
    msgs <- c(msgs, .check_path(sprintf("dmr_files[['%s']]", tp), dmr_files[[tp]]))
    msgs <- c(msgs, .check_path(sprintf("dmr_merged_files[['%s']]$hyper", tp), dmr_merged_files[[tp]]$hyper))
    msgs <- c(msgs, .check_path(sprintf("dmr_merged_files[['%s']]$hypo", tp), dmr_merged_files[[tp]]$hypo))
    msgs <- c(msgs, .check_path(sprintf("dmr_unfiltered_files[['%s']]", tp), dmr_unfiltered_files[[tp]]))
    msgs <- c(msgs, .check_path(sprintf("dmg_promoter_files[['%s']]$hyper", tp), dmg_promoter_files[[tp]]$hyper))
    msgs <- c(msgs, .check_path(sprintf("dmg_promoter_files[['%s']]$hypo", tp), dmg_promoter_files[[tp]]$hypo))
    msgs <- c(msgs, .check_path(sprintf("deg_files[['%s']]", tp), deg_files[[tp]]))
  }
  msgs <- c(msgs, .check_path("coverage_summary_file", coverage_summary_file))
  for (tp in TIMEPOINTS) {
    for (s in names(raw_cpg_report_files[[tp]])) {
      msgs <- c(msgs, .check_path(sprintf("raw_cpg_report_files[['%s']]['%s']", tp, s),
                                   raw_cpg_report_files[[tp]][[s]]))
    }
  }
  msgs <- c(msgs, .check_path("master_fpkm_file", master_fpkm_file))
  msgs <- c(msgs, .check_path("gene_body_bed_file", gene_body_bed_file))
  msgs <- c(msgs, .check_path("te_catalog_bed_file", te_catalog_bed_file))
  msgs <- c(msgs, .check_path("te_catalog_meta_file", te_catalog_meta_file))

  msgs <- Filter(Negate(is.null), msgs)

  if (length(msgs) == 0) {
    cat("config.R: all input paths resolve.\n")
    cat("Promoter DMGs use dmg_promoter_files; TE context uses te_context_promoter_bed_file.\n")
  } else {
    cat("config.R: still missing", length(msgs), "item(s):\n")
    cat(paste(msgs, collapse = "\n"), "\n")
  }
  invisible(length(msgs) == 0)
}

.validate_config()
