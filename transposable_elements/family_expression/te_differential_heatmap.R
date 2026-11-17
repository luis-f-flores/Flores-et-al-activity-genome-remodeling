#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


packages <- c("readr","dplyr","tidyr","stringr","ggplot2","forcats","scales","purrr")

invisible(lapply(packages, require, character.only = TRUE))


args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: te_differential_heatmap.R SalmonTE_results_root output_directory")
base_dir <- args[1]
baseMean_cutoff <- 30
color_cap <- 3
focus_IAP_only <- FALSE
focus_component <- "all"
limit_sigonly_topN <- NA_integer_


facet_mode <- "none"
row_order_mode <- "alpha"


out_dir <- args[2]
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


files <- sort(Sys.glob(file.path(base_dir, "SalmonTE_DE_*", "TE_DE_results.csv")))
if (length(files) == 0) stop("No TE_DE_results.csv files found under ", normalizePath(base_dir))

label_from_dir <- function(path) {
  tpt <- sub("^SalmonTE_DE_", "", basename(dirname(path)))
  paste0(tpt, " KCl vs Mann")
}
contrast_labels <- vapply(files, label_from_dir, character(1))


preferred <- c("1h KCl vs Mann","6h KCl vs Mann","24h KCl vs Mann")
ord <- order(match(contrast_labels, preferred), na.last = TRUE)
files <- files[ord]; contrast_labels <- contrast_labels[ord]

message("Using CSV result files:\n",
        paste(sprintf("  %s -> %s", basename(dirname(files)), basename(files)), collapse = "\n"))


canonize_names <- function(nms) {
  nms2 <- stringr::str_replace_all(nms, " ", " ")
  nms2 <- stringr::str_trim(nms2)
  nms2 <- stringr::str_replace_all(nms2, "[^A-Za-z0-9]+", "")
  tolower(nms2)
}


classify_type <- function(fam) {
  fam <- stringr::str_replace_all(fam, "−", "-")
  dplyr::case_when(

    stringr::str_detect(fam, "(^IAPLTR|^RLTR|^MLTR|(^|[-_])LTR([-_]|$)|\\bLTR\\b)") ~ "LTR",

    stringr::str_detect(fam, "^IAP") ~ "internal",

    stringr::str_detect(fam, "([_-]|\\b)(I|int)([_-]|\\b)") ~ "internal",
    stringr::str_detect(fam, "_I$") ~ "internal",
    TRUE ~ "other"
  )
}


star_for_padj <- function(padj) {
  dplyr::case_when(
    is.na(padj)   ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

`%||%` <- function(a,b) if (!is.null(a)) a else b


load_one <- function(f, label) {
  df <- readr::read_csv(f, show_col_types = FALSE)


  if (!"family" %in% names(df)) names(df)[1] <- "family"


  req <- c("family","baseMean","log2FoldChange","pvalue","padj")
  missing_req <- setdiff(req, names(df))
  if (length(missing_req)) {

    nm <- setNames(names(df), canonize_names(names(df)))
    revmap <- function(key) names(nm)[match(key, names(nm))]
    mapto <- list(
      family         = revmap("family") %||% names(df)[1],
      baseMean       = revmap("basemean"),
      log2FoldChange = revmap("log2foldchange"),
      pvalue         = revmap("pvalue"),
      padj           = revmap("padj")
    )
    for (k in names(mapto)) {
      if (!is.null(mapto[[k]]) && !(k %in% names(df))) {
        names(df)[names(df) == mapto[[k]]] <- k
      }
    }
    missing_req <- setdiff(req, names(df))
    if (length(missing_req)) {
      stop("Missing required columns: ", paste(missing_req, collapse = ", "),
           "\nAvailable columns: ", paste(names(df), collapse = ", "))
    }
  }

  df <- df |>
    dplyr::transmute(
      family         = as.character(.data$family),
      baseMean       = suppressWarnings(as.numeric(.data$baseMean)),
      log2FoldChange = suppressWarnings(as.numeric(.data$log2FoldChange)),
      pvalue         = suppressWarnings(as.numeric(.data$pvalue)),
      padj           = suppressWarnings(as.numeric(.data$padj)),
      contrast       = label
    )

  df$family <- stringr::str_replace_all(df$family, "−", "-")
  df$type   <- classify_type(df$family)


  df$star <- ifelse(
    !is.na(df$baseMean) & df$baseMean >= baseMean_cutoff,
    star_for_padj(df$padj),
    ""
  )


  df$sig_padj <- !is.na(df$padj) & df$padj <= 0.05 &
                 !is.na(df$baseMean) & df$baseMean >= baseMean_cutoff
  df
}


tbl <- purrr::map2_dfr(files, contrast_labels, load_one)

type_override_file <- file.path(base_dir, "TE_type_overrides.csv")
if (file.exists(type_override_file)) {
  overrides <- readr::read_csv(type_override_file, show_col_types = FALSE) %>%
    dplyr::select(family, type_override = type) %>%
    dplyr::mutate(type_override = factor(type_override, levels = c("internal","LTR","other")))

  tbl <- tbl %>%
    dplyr::left_join(overrides, by = "family") %>%
    dplyr::mutate(
      type = dplyr::if_else(
        !is.na(.data$type_override),
        as.character(.data$type_override),
        as.character(.data$type)
      )
    ) %>%
    dplyr::mutate(type = factor(type, levels = c("internal","LTR","other"))) %>%
    dplyr::select(-type_override)

  message("Applied manual type overrides from ", basename(type_override_file))
} else {
  message("No TE_type_overrides.csv found in ", normalizePath(base_dir),
          "; you can add one with columns: family,type (internal|LTR|other)")
}


if (isTRUE(focus_IAP_only)) {
  tbl <- dplyr::filter(tbl, stringr::str_detect(.data$family, "^IAP"))
  if (identical(focus_component, "LTR")) {
    tbl <- dplyr::filter(tbl, .data$type == "LTR")
  } else if (identical(focus_component, "internal")) {
    tbl <- dplyr::filter(tbl, .data$type == "internal")
  }
}


readr::write_csv(
  tbl %>%
    dplyr::distinct(family, type) %>%
    dplyr::count(type, name = "n") %>%
    dplyr::arrange(desc(n)),
  file.path(out_dir, "TE_type_summary.csv")
)
writeLines(
  tbl %>%
    dplyr::filter(type == "other") %>%
    dplyr::distinct(family) %>%
    dplyr::arrange(family) %>%
    dplyr::pull(family),
  con = file.path(out_dir, "TE_other_families.txt")
)

keepers_full <- tbl |>
  dplyr::group_by(.data$family, .data$type) |>
  dplyr::summarize(
    any_ok = any(.data$baseMean >= baseMean_cutoff, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$any_ok) |>
  dplyr::pull(.data$family)

full_tbl <- dplyr::filter(tbl, .data$family %in% keepers_full)


if (identical(row_order_mode, "alpha")) {
  full_row_order <- full_tbl |>
    dplyr::distinct(family) |>
    dplyr::arrange(family) |>
    dplyr::pull(family)
} else {
  full_row_order <- full_tbl |>
    dplyr::group_by(.data$family, .data$type) |>
    dplyr::summarize(
      max_abs = suppressWarnings(max(abs(.data$log2FoldChange), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(dplyr::coalesce(.data$max_abs, -Inf)), .data$family) |>
    dplyr::pull(.data$family)
}

full_tbl$family   <- factor(full_tbl$family, levels = unique(full_row_order))
full_tbl$contrast <- factor(full_tbl$contrast, levels = contrast_labels)
full_tbl$type     <- factor(full_tbl$type, levels = c("internal","LTR","other"))

make_prefix <- function(x) sub("[_-].*$", "", x)
full_tbl$prefix <- make_prefix(as.character(full_tbl$family))

has_sig <- tbl |>
  dplyr::group_by(.data$family) |>
  dplyr::summarize(
    any_sig = any(.data$sig_padj, na.rm = TRUE),
    .groups = "drop"
  )

sigonly_tbl <- tbl |>
  dplyr::inner_join(has_sig, by = "family") |>
  dplyr::filter(.data$any_sig) |>
  dplyr::select(-any_sig)


if (is.integer(limit_sigonly_topN) && !is.na(limit_sigonly_topN) && limit_sigonly_topN > 0) {
  top_fams <- sigonly_tbl |>
    dplyr::group_by(.data$family, .data$type) |>
    dplyr::summarize(
      max_abs = suppressWarnings(max(abs(.data$log2FoldChange), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$max_abs)) |>
    dplyr::slice_head(n = limit_sigonly_topN) |>
    dplyr::pull(.data$family)
  sigonly_tbl <- dplyr::filter(sigonly_tbl, .data$family %in% top_fams)
}

if (identical(row_order_mode, "alpha")) {
  sig_row_order <- sigonly_tbl |>
    dplyr::distinct(family) |>
    dplyr::arrange(family) |>
    dplyr::pull(family)
} else {
  sig_row_order <- sigonly_tbl |>
    dplyr::group_by(.data$family, .data$type) |>
    dplyr::summarize(
      max_abs = suppressWarnings(max(abs(.data$log2FoldChange), na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(dplyr::coalesce(.data$max_abs, -Inf)), .data$family) |>
    dplyr::pull(.data$family)
}

sigonly_tbl$family   <- factor(sigonly_tbl$family, levels = unique(sig_row_order))
sigonly_tbl$contrast <- factor(sigonly_tbl$contrast, levels = contrast_labels)
sigonly_tbl$type     <- factor(sigonly_tbl$type, levels = c("internal","LTR","other"))
sigonly_tbl$prefix   <- make_prefix(as.character(sigonly_tbl$family))


write_te_heatmap_tables <- function(df, file_stem) {
  mat_log2fc <- df |>
    dplyr::select(family, type, contrast, log2FoldChange) |>
    tidyr::pivot_wider(names_from = contrast, values_from = log2FoldChange) |>
    dplyr::arrange(match(family, levels(df$family)))
  readr::write_csv(mat_log2fc, file.path(out_dir, paste0(file_stem, "_log2FC_matrix.csv")))

  mat_star <- df |>
    dplyr::select(family, type, contrast, star) |>
    tidyr::pivot_wider(names_from = contrast, values_from = star) |>
    dplyr::arrange(match(family, levels(df$family)))
  readr::write_csv(mat_star, file.path(out_dir, paste0(file_stem, "_stars_matrix.csv")))
}


write_te_heatmap_tables(sigonly_tbl, "TE_DE_heatmap_SIGONLY")
