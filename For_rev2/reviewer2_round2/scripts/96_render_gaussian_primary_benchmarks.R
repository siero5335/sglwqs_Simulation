if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("The rmarkdown package is required to render the Gaussian primary benchmark report.", call. = FALSE)
}
output_name <- Sys.getenv("R2R2_PRIMARY_REPORT_DIR", "gaussian_primary_benchmarks")
output_file <- if (identical(output_name, "gaussian_primary_benchmarks")) {
  "gaussian_primary_benchmarks.html"
} else {
  paste0(output_name, ".html")
}
rmarkdown::render(
  "reviewer2_round2/gaussian_primary_benchmarks.Rmd",
  output_file = output_file,
  output_dir = file.path("reviewer2_round2/results/summaries", output_name),
  knit_root_dir = getwd(),
  envir = new.env(parent = globalenv()),
  quiet = TRUE
)
