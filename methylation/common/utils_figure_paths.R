#!/usr/bin/env Rscript
GB_WORK_ROOT <- Sys.getenv("GB_WORK_ROOT", unset = file.path(getwd(), "work"))
gb_path <- function(name) {
  if (startsWith(name, "@tool/")) return(basename(name))
  if (startsWith(name, "@fonts")) return(file.path(Sys.getenv("GB_FONT_DIR", unset="/usr/share/fonts/truetype/msttcorefonts"), sub("^@fonts/?", "", name)))
  file.path(GB_WORK_ROOT, name)
}

FIGURE_PACKAGE <- "WGBS"
FIGURE_LINUX_ROOT <- gb_path("figures")
FIGURE_SOURCE_LINUX_ROOT <- gb_path("wgbs_package")
FIGURE_WORK_ROOT <- FIGURE_LINUX_ROOT
revision_figure_dir <- function(run_id) {
  stopifnot(length(run_id)==1L, grepl("^[A-Za-z0-9_-]+$",run_id))
  target <- file.path(PROJECT_ROOT,"runs",run_id,"figures")
  dir.create(target,recursive=TRUE,showWarnings=FALSE)
  target
}
