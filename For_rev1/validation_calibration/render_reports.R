#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(profile = "smoke", reports = character(0))
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--profile")) {
      i <- i + 1L
      out$profile <- args[[i]]
    } else if (identical(arg, "--reports")) {
      i <- i + 1L
      out$reports <- strsplit(args[[i]], ",", fixed = TRUE)[[1]]
    } else if (arg %in% c("-h", "--help")) {
      cat("Usage: Rscript analysis/validation_calibration/render_reports.R --profile smoke\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
if (is.na(file_arg) || !nzchar(file_arg)) {
  file_arg <- "analysis/validation_calibration/render_reports.R"
}
root <- normalizePath(dirname(file_arg), mustWork = TRUE)

source(file.path(root, "R", "00_paths.R"))
Sys.setenv(VALIDATION_CALIBRATION_ROOT = root)
activate_local_lib()
ensure_vc_dirs()

if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render validation calibration reports.", call. = FALSE)
}

all_reports <- c(
  "01_smoke_and_runtime.Rmd",
  "02_type1_global_null.Rmd",
  "03_partial_null_false_attribution.Rmd",
  "04_split_stability.Rmd",
  "05_gaussian_sensitivity.Rmd",
  "06_runtime_projection.Rmd",
  "07_combined_report.Rmd"
)
reports <- if (length(args$reports)) args$reports else all_reports

for (report in reports) {
  input <- file.path(root, "reports", report)
  output_file <- sub("\\.Rmd$", paste0("_", args$profile, ".html"), report)
  message("Rendering ", report)
  rmarkdown::render(
    input = input,
    output_file = output_file,
    output_dir = file.path(root, "output", "html"),
    params = list(profile = args$profile, root = root),
    envir = new.env(parent = globalenv()),
    quiet = TRUE
  )
}
