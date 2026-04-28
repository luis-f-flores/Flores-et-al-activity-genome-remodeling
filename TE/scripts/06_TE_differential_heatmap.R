#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_TE_differential_heatmap.R
# Package: Transposable element analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for TE catalog construction, TE differential-expression annotation, TE-gene overlap, IAP summaries, methylation metaprofiles, and heatmaps.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# 06_TE_differential_heatmap.R
#
# Generate TE-family differential-expression heatmaps from SalmonTE/DESeq2
# result tables.
#

# Inputs:
#   A parent directory containing one subdirectory per contrast:
#       SalmonTE_DE_<label>/TE_DE_results.csv
#
# Each TE_DE_results.csv must contain at least:
#   family, baseMean, log2FoldChange, pvalue, padj
#

# Outputs:
#   TE_DE_summary/TE_DE_heatmap_FULL.pdf
#   TE_DE_summary/TE_DE_heatmap_FULL.jpg
#   TE_DE_summary/TE_DE_heatmap_FULL_log2FC_matrix.csv
#   TE_DE_summary/TE_DE_heatmap_FULL_stars_matrix.csv
#   TE_DE_summary/TE_DE_heatmap_SIGONLY.pdf
#   TE_DE_summary/TE_DE_heatmap_SIGONLY.jpg
#   TE_DE_summary/TE_DE_heatmap_SIGONLY_log2FC_matrix.csv
#   TE_DE_summary/TE_DE_heatmap_SIGONLY_stars_matrix.csv
#
# Example:
#   Rscript 06_TE_differential_heatmap.R \
#     --base-dir /path/to/TE_results \
#     --out-dir /path/to/TE_results/TE_DE_summary \
#     --baseMean-cutoff 30 \
#     --color-cap 3 \
#     --row-order alpha \
#     --facet none
#
# Notes:
#   - Heatmap fill is log2FoldChange for KCl vs mannitol.
#   - Asterisks mark FDR significance only for cells passing the baseMean cutoff.
#   - FULL heatmap includes TE families passing baseMean cutoff in at least one contrast.
#   - SIGONLY heatmap includes TE families significant in at least one contrast.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(purrr)
})

# ----------------------------- argument parser -----------------------------

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx)) return(default)
  if (idx == length(args)) stop("Missing value after ", flag)
  args[[idx + 1]]
}

to_logical <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

base_dir <- get_arg("--base-dir", ".")
out_dir <- get_arg("--out-dir", file.path(base_dir, "TE_DE_summary"))
pattern <- get_arg("--pattern", "SalmonTE_DE_*/TE_DE_results.csv")

baseMean_cutoff <- as.numeric(get_arg("--baseMean-cutoff", "30"))
color_cap <- as.numeric(get_arg("--color-cap", "3"))

focus_IAP_only <- to_logical(get_arg("--focus-IAP-only", "FALSE"))
focus_component <- get_arg("--focus-component", "all")  # all | LTR | internal
limit_sigonly_topN_arg <- get_arg("--limit-sigonly-topN", "NA")
limit_sigonly_topN <- suppressWarnings(as.integer(limit_sigonly_topN_arg))

facet_mode <- get_arg("--facet", "none")       # none | type | prefix
row_order_mode <- get_arg("--row-order", "alpha") # alpha | effect

if (!dir.exists(base_dir)) stop("base_dir does not exist: ", base_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!facet_mode %in% c("none", "type", "prefix")) {
  stop("--facet must be one of: none, type, prefix")
}
if (!row_order_mode %in% c("alpha", "effect")) {
  stop("--row-order must be one of: alpha, effect")
}
if (!focus_component %in% c("all", "LTR", "internal")) {
  stop("--focus-component must be one of: all, LTR, internal")
}

# ----------------------------- helpers -------------------------------------

canonize_names <- function(nms) {
  nms |>
    stringr::str_replace_all("\u00A0", " ") |>
    stringr::str_trim() |>
    stringr::str_replace_all("[^A-Za-z0-9]+", "") |>
    tolower()
}

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

label_from_dir <- function(path) {
  label <- sub("^SalmonTE_DE_", "", basename(dirname(path)))
  paste0(label, " KCl vs Mannitol")
}

classify_type <- function(fam) {
  fam <- stringr::str_replace_all(fam, "\u2212", "-")
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
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01 ~ "**",
    padj < 0.05 ~ "*",
    TRUE ~ ""
  )
}

make_prefix <- function(x) sub("[_-].*$", "", x)

# ----------------------------- discover files ------------------------------

files <- sort(Sys.glob(file.path(base_dir, pattern)))
if (length(files) == 0) {
  stop("No TE_DE_results.csv files found using pattern: ", file.path(base_dir, pattern))
}

contrast_labels <- vapply(files, label_from_dir, character(1))

# Preferred ordering if the usual timepoints are present.
preferred <- c("1h KCl vs Mannitol", "6h KCl vs Mannitol", "24h KCl vs Mannitol",
               "1hr KCl vs Mannitol", "6hr KCl vs Mannitol", "24hr KCl vs Mannitol")
ord <- order(match(contrast_labels, preferred), contrast_labels, na.last = TRUE)
files <- files[ord]
contrast_labels <- contrast_labels[ord]

message("Using TE DE result files:")
message(paste(sprintf("  %s -> %s", dirname(files), basename(files)), collapse = "\n"))

# ----------------------------- load data -----------------------------------

load_one <- function(f, label) {
  df <- readr::read_csv(f, show_col_types = FALSE)

  if (!"family" %in% names(df)) names(df)[1] <- "family"

  required <- c("family", "baseMean", "log2FoldChange", "pvalue", "padj")
  missing_required <- setdiff(required, names(df))

  if (length(missing_required) > 0) {
    canonical <- setNames(names(df), canonize_names(names(df)))

    find_name <- function(key) {
      out <- names(canonical)[match(key, names(canonical))]
      if (is.na(out)) NULL else out
    }

    mapping <- list(
      family = find_name("family") %||% names(df)[1],
      baseMean = find_name("basemean"),
      log2FoldChange = find_name("log2foldchange"),
      pvalue = find_name("pvalue"),
      padj = find_name("padj")
    )

    for (target in names(mapping)) {
      old <- mapping[[target]]
      if (!is.null(old) && !(target %in% names(df))) {
        names(df)[names(df) == old] <- target
      }
    }

    missing_required <- setdiff(required, names(df))
    if (length(missing_required) > 0) {
      stop(
        "Missing required columns in ", f, ": ",
        paste(missing_required, collapse = ", "),
        "\nAvailable columns: ", paste(names(df), collapse = ", ")
      )
    }
  }

  df |>
    transmute(
      family = as.character(.data$family),
      baseMean = suppressWarnings(as.numeric(.data$baseMean)),
      log2FoldChange = suppressWarnings(as.numeric(.data$log2FoldChange)),
      pvalue = suppressWarnings(as.numeric(.data$pvalue)),
      padj = suppressWarnings(as.numeric(.data$padj)),
      contrast = label
    ) |>
    mutate(
      family = stringr::str_replace_all(.data$family, "\u2212", "-"),
      type = classify_type(.data$family),
      star = ifelse(
        !is.na(.data$baseMean) & .data$baseMean >= baseMean_cutoff,
        star_for_padj(.data$padj),
        ""
      ),
      sig_padj = !is.na(.data$padj) &
        .data$padj <= 0.05 &
        !is.na(.data$baseMean) &
        .data$baseMean >= baseMean_cutoff
    )
}

tbl <- purrr::map2_dfr(files, contrast_labels, load_one)

# Optional manual family type overrides.
override_file <- file.path(base_dir, "TE_type_overrides.csv")
if (file.exists(override_file)) {
  overrides <- readr::read_csv(override_file, show_col_types = FALSE) |>
    select(family, type_override = type)

  tbl <- tbl |>
    left_join(overrides, by = "family") |>
    mutate(
      type = if_else(!is.na(.data$type_override), .data$type_override, .data$type),
      type = factor(.data$type, levels = c("internal", "LTR", "other"))
    ) |>
    select(-type_override)

  message("Applied TE type overrides from: ", override_file)
} else {
  tbl <- tbl |> mutate(type = factor(.data$type, levels = c("internal", "LTR", "other")))
}

if (focus_IAP_only) {
  tbl <- tbl |> filter(stringr::str_detect(.data$family, "^IAP"))

  if (focus_component == "LTR") {
    tbl <- tbl |> filter(.data$type == "LTR")
  } else if (focus_component == "internal") {
    tbl <- tbl |> filter(.data$type == "internal")
  }
}

# ----------------------------- summaries -----------------------------------

readr::write_csv(
  tbl |>
    distinct(family, type) |>
    count(type, name = "n") |>
    arrange(desc(n)),
  file.path(out_dir, "TE_type_summary.csv")
)

writeLines(
  tbl |>
    filter(.data$type == "other") |>
    distinct(family) |>
    arrange(family) |>
    pull(family),
  con = file.path(out_dir, "TE_other_families.txt")
)

# ----------------------------- row sets -------------------------------------

keepers_full <- tbl |>
  group_by(family) |>
  summarize(any_ok = any(.data$baseMean >= baseMean_cutoff, na.rm = TRUE), .groups = "drop") |>
  filter(.data$any_ok) |>
  pull(family)

full_tbl <- tbl |> filter(.data$family %in% keepers_full)

has_sig <- tbl |>
  group_by(family) |>
  summarize(any_sig = any(.data$sig_padj, na.rm = TRUE), .groups = "drop")

sigonly_tbl <- tbl |>
  inner_join(has_sig, by = "family") |>
  filter(.data$any_sig) |>
  select(-any_sig)

if (!is.na(limit_sigonly_topN) && limit_sigonly_topN > 0) {
  top_families <- sigonly_tbl |>
    group_by(family) |>
    summarize(max_abs = max(abs(.data$log2FoldChange), na.rm = TRUE), .groups = "drop") |>
    arrange(desc(.data$max_abs)) |>
    slice_head(n = limit_sigonly_topN) |>
    pull(family)

  sigonly_tbl <- sigonly_tbl |> filter(.data$family %in% top_families)
}

order_families <- function(df) {
  if (row_order_mode == "alpha") {
    df |> distinct(family) |> arrange(family) |> pull(family)
  } else {
    df |>
      group_by(family) |>
      summarize(max_abs = max(abs(.data$log2FoldChange), na.rm = TRUE), .groups = "drop") |>
      arrange(desc(.data$max_abs), family) |>
      pull(family)
  }
}

prepare_plot_df <- function(df) {
  fam_order <- order_families(df)

  df |>
    mutate(
      family = factor(.data$family, levels = unique(fam_order)),
      contrast = factor(.data$contrast, levels = contrast_labels),
      type = factor(.data$type, levels = c("internal", "LTR", "other")),
      prefix = make_prefix(as.character(.data$family))
    )
}

full_tbl <- prepare_plot_df(full_tbl)
sigonly_tbl <- prepare_plot_df(sigonly_tbl)

# ----------------------------- plotting -------------------------------------

plot_te_heatmap <- function(df, title, file_stem) {
  if (nrow(df) == 0) {
    warning("No rows to plot for: ", file_stem)
    return(invisible(NULL))
  }

  n_rows <- length(levels(df$family))
  fig_height <- max(6, 0.25 * n_rows + 2)
  fig_width <- 8.5

  p <- ggplot(df, aes(x = contrast, y = family, fill = log2FoldChange)) +
    geom_tile(color = NA, na.rm = TRUE) +
    scale_fill_gradient2(
      low = "#3b70f2",
      mid = "white",
      high = "#d73027",
      midpoint = 0,
      limits = c(-color_cap, color_cap),
      oob = scales::squish,
      na.value = "grey90",
      name = "log2FC\n(KCl vs Mannitol)"
    ) +
    geom_text(
      data = filter(df, .data$star != ""),
      aes(label = star),
      size = 3,
      color = "black",
      na.rm = TRUE
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL,
      caption = paste0(
        "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001. ",
        "Asterisks shown only where baseMean >= ", baseMean_cutoff, "."
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      plot.title = element_text(face = "bold", size = 14)
    )

  if (facet_mode == "type") {
    p <- p + facet_grid(rows = vars(type), scales = "free_y", space = "free_y", switch = "y")
  } else if (facet_mode == "prefix") {
    p <- p + facet_grid(rows = vars(prefix), scales = "free_y", space = "free_y", switch = "y")
  }

  pdf_file <- file.path(out_dir, paste0(file_stem, ".pdf"))
  jpg_file <- file.path(out_dir, paste0(file_stem, ".jpg"))

  ggsave(pdf_file, p, width = fig_width, height = fig_height, units = "in")
  ggsave(jpg_file, p, width = fig_width, height = fig_height, units = "in", dpi = 300)

  log2fc_matrix <- df |>
    select(family, type, contrast, log2FoldChange) |>
    pivot_wider(names_from = contrast, values_from = log2FoldChange) |>
    arrange(match(family, levels(df$family)))

  stars_matrix <- df |>
    select(family, type, contrast, star) |>
    pivot_wider(names_from = contrast, values_from = star) |>
    arrange(match(family, levels(df$family)))

  readr::write_csv(log2fc_matrix, file.path(out_dir, paste0(file_stem, "_log2FC_matrix.csv")))
  readr::write_csv(stars_matrix, file.path(out_dir, paste0(file_stem, "_stars_matrix.csv")))

  message("Saved: ", pdf_file)
  message("Saved: ", jpg_file)
}

plot_te_heatmap(
  full_tbl,
  title = "TE differential expression - full context",
  file_stem = "TE_DE_heatmap_FULL"
)

plot_te_heatmap(
  sigonly_tbl,
  title = "TE differential expression - FDR-significant families",
  file_stem = "TE_DE_heatmap_SIGONLY"
)

message("Done. Outputs written to: ", normalizePath(out_dir))
