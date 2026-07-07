#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(profile = "smoke")
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--profile")) {
      i <- i + 1L
      out$profile <- args[[i]]
    } else if (arg %in% c("-h", "--help")) {
      cat("Usage: Rscript validation_calibration/render_non_split_reports.R --profile production\n")
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
  file_arg <- "validation_calibration/render_non_split_reports.R"
}
root <- normalizePath(dirname(file_arg), mustWork = TRUE)

reports <- c(
  "01_smoke_and_runtime.Rmd",
  "02_type1_global_null.Rmd",
  "03_partial_null_false_attribution.Rmd",
  "05_gaussian_sensitivity.Rmd",
  "06_runtime_projection.Rmd",
  "07_non_split_combined_report.Rmd"
)

for (report in reports) {
  input <- file.path(root, "reports", report)
  output_file <- paste0(tools::file_path_sans_ext(report), "_", args$profile, ".html")
  rmarkdown::render(
    input,
    params = list(profile = args$profile, root = root),
    output_dir = file.path(root, "output", "html"),
    output_file = output_file,
    quiet = FALSE
  )
}

