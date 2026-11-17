#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


suppressPackageStartupMessages({
  library(data.table)
})

base_dirs <- c(
  gb_path("wgbs/methylKit_DMR/DMR_1M_vs_1K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_6M_vs_6K_multiwindow/DMR_win_500bp"),
  gb_path("wgbs/methylKit_DMR/DMR_24M_vs_24K_multiwindow/DMR_win_500bp")
)

merged_subdir <- "Merged_q05_meth5"
out_subdir    <- "SingleFeature_q05_meth5"


feature_priority <- c("Promoter", "Exon", "Intron", "Intergenic", "5'UTR", "3'UTR")

map_feature_category <- function(x) {
  x <- tolower(x)
  f <- ifelse(grepl("promoter|tss", x), "Promoter",
       ifelse(grepl("fiveutr|5utr", x), "5'UTR",
       ifelse(grepl("threeutr|3utr", x), "3'UTR",
       ifelse(grepl("exon", x), "Exon",
       ifelse(grepl("intron", x), "Intron",
       ifelse(grepl("intergenic", x), "Intergenic",

       ifelse(grepl("\\bgene\\b", x), NA_character_,
              NA_character_)))))))
  factor(f, levels = feature_priority)
}


summary_all <- list()

collapse_one_contrast <- function(base_dir) {
  base_dir <- path.expand(base_dir)
  if (!dir.exists(base_dir)) {
    cat("Base dir not found:", base_dir, "\n")
    return(invisible(NULL))
  }


  comp_dir  <- basename(dirname(base_dir))
  comparison <- sub("^DMR_", "", sub("_multiwindow$", "", comp_dir))

  cat("\nContrast:", comparison, "\n")
  merged_dir <- file.path(base_dir, merged_subdir)
  if (!dir.exists(merged_dir)) {
    cat("Merged_q05_meth5 dir not found:", merged_dir, "\n")
    return(invisible(NULL))
  }

  out_dir <- file.path(merged_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  for (cls in c("hyper", "hypo")) {
    cat("  Class:", cls, "\n")

    raw_file   <- file.path(merged_dir, sprintf("%s_Merged_DMRs_cov4_%s.txt", comparison, cls))
    annot_file <- file.path(merged_dir, sprintf("GTF_Annot_%s_Merged_DMRs_cov4_%s.txt", comparison, cls))

    if (!file.exists(raw_file)) {
      cat("     Missing raw DMR file:", raw_file, "\n")
      next
    }
    if (!file.exists(annot_file)) {
      cat("     Missing annotation file:", annot_file, "\n")
      next
    }


    dmrs <- fread(raw_file)
    req_raw <- c("chr", "start", "end")
    if (!all(req_raw %in% names(dmrs))) {
      stop("Raw DMR file missing chr/start/end: ", raw_file)
    }


    if ("meth.diff" %in% names(dmrs)) setnames(dmrs, "meth.diff", "meth_diff")
    if (!"meth_diff" %in% names(dmrs)) stop("No meth_diff column in ", raw_file)
    if (!"qvalue"   %in% names(dmrs)) stop("No qvalue column in ", raw_file)

    dmrs[, dmr_id := paste(chr, start, end, sep = "_")]
    raw_n   <- nrow(dmrs)
    raw_ids <- dmrs$dmr_id


    ann <- fread(annot_file)
    req_ann <- c("chr", "start", "end", "strand", "feature")
    if (!all(req_ann %in% names(ann))) {
      stop("Annotation file missing chr/start/end/strand/feature: ", annot_file)
    }

    ann[, dmr_id := paste(chr, start, end, sep = "_")]
    ann[, feature_category := as.character(map_feature_category(feature))]


    setDT(ann)


    collapsed <- ann[!is.na(feature_category),
      {
        cats <- unique(feature_category)

        idx  <- match(feature_priority, cats, nomatch = Inf)
        best <- feature_priority[which.min(idx)]

        gn <- NA_character_
        if ("gene_name" %in% names(ann)) {
          gn_vec <- unique(gene_name[feature_category == best])
          gn_vec <- gn_vec[!is.na(gn_vec) & gn_vec != ""]
          if (length(gn_vec)) gn <- paste(gn_vec, collapse = ",")
        }
        .(strand = strand[1L],
          feature_category = best,
          gene_name = gn)
      },
      by = .(dmr_id, chr, start, end)
    ]


    assigned_ids <- collapsed$dmr_id
    missing_ids  <- setdiff(raw_ids, assigned_ids)

    if (length(missing_ids) > 0) {
      intergenic <- dmrs[dmr_id %in% missing_ids,
                         .(dmr_id, chr, start, end)]
      intergenic[, `:=`(
        strand = "*",
        feature_category = "Intergenic",
        gene_name = NA_character_
      )]
      collapsed <- rbind(collapsed, intergenic, use.names = TRUE, fill = TRUE)
    }


    collapsed <- merge(
      collapsed,
      dmrs[, .(dmr_id, meth_diff, qvalue)],
      by = "dmr_id",
      all.x = TRUE
    )


    if (nrow(collapsed) != raw_n) {
      stop("Row count mismatch for ", comparison, " ", cls,
           ": raw_n = ", raw_n, ", collapsed = ", nrow(collapsed))
    }


    collapsed[, `:=`(comparison = comparison, class = cls)]
    setcolorder(
      collapsed,
      c("comparison", "class",
        "chr", "start", "end", "strand",
        "feature_category",
        "meth_diff", "qvalue",
        "gene_name", "dmr_id")
    )


    summary_dt <- collapsed[
      , .N, by = feature_category
    ][order(factor(feature_category, levels = feature_priority))]

    summary_dt[, `:=`(comparison = comparison, class = cls)]
    setcolorder(summary_dt, c("comparison", "class", "feature_category", "N"))


    out_file <- file.path(
      out_dir,
      sprintf("%s_Merged_DMRs_cov4_%s_singleFeature_q05_meth5.tsv",
              comparison, cls)
    )
    sum_file <- file.path(
      out_dir,
      sprintf("Summary_SingleFeature_%s_%s_q05_meth5.tsv",
              comparison, cls)
    )

    fwrite(collapsed[, !c("dmr_id")], out_file,
           sep = "\t", quote = FALSE, na = "NA")
    fwrite(summary_dt, sum_file,
           sep = "\t", quote = FALSE, na = "NA")

    cat("    Collapsed file:", out_file, "\n")
    cat("    Summary file:  ", sum_file, "\n")


    summary_all[[paste(comparison, cls, sep = "_")]] <<- summary_dt
  }

  invisible(TRUE)
}

for (bd in base_dirs) {
  collapse_one_contrast(bd)
}

if (length(summary_all)) {
  all_dt <- rbindlist(summary_all, use.names = TRUE, fill = TRUE)
  out_global <- "SingleFeature_counts_all_q05_meth5.tsv"
  fwrite(all_dt, out_global, sep = "\t", quote = FALSE, na = "NA")
  cat("\nGlobal summary written:", out_global, "\n")
} else {
  cat("\n No summaries produced (check inputs).\n")
}

cat("\nSingle-feature q05/meth5 assignments written.\n")
