#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}


minmax_scale_pm0.5 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng) - 0.5
}

fit_methylation_expression_regression <- function(df, fpkm_filter = NULL) {
  timepoints <- unique(df$timepoint)
  results <- lapply(timepoints, function(tp) {
    sub <- df[df$timepoint == tp, ]
    if (!is.null(fpkm_filter)) {
      sub <- sub[sub$gene_id %in% fpkm_filter, ]
    }
    sub <- sub[stats::complete.cases(sub[, c("log2FoldChange", "delta_mCG")]), ]

    if (nrow(sub) < 3) {
      return(data.frame(
        timepoint = tp, n = nrow(sub), slope = NA, intercept = NA,
        r_squared = NA, adj_r_squared = NA, p_value = NA
      ))
    }

    sub$delta_mCG_scaled <- minmax_scale_pm0.5(sub$delta_mCG)
    fit <- stats::lm(log2FoldChange ~ delta_mCG_scaled, data = sub)
    s <- summary(fit)

    data.frame(
      timepoint     = tp,
      n             = nrow(sub),
      slope         = stats::coef(fit)[["delta_mCG_scaled"]],
      intercept     = stats::coef(fit)[["(Intercept)"]],
      r_squared     = s$r.squared,
      adj_r_squared = s$adj.r.squared,
      p_value       = stats::coef(s)["delta_mCG_scaled", "Pr(>|t|)"]
    )
  })
  do.call(rbind, results)
}

compare_regressions <- function(result_a, result_b, label_a = "A", label_b = "B") {
  merged <- merge(result_a, result_b, by = "timepoint", suffixes = paste0("_", c(label_a, label_b)))
  merged
}
