source("reviewer2_round2/R/design_helpers.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("method_wrappers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

manifest_path <- r2r2_file("config", "scenario_manifest.csv")
manifest <- if (file.exists(manifest_path)) {
  readr::read_csv(manifest_path, show_col_types = FALSE)
} else {
  make_production_manifest()
}
status <- summarize_manifest_status(manifest)
cat("Valid completed jobs:", sum(status$valid_complete), "/", nrow(status), "\n")
print(status |>
  filter(.data$valid_complete) |>
  count(.data$battery, .data$method, name = "done") |>
  arrange(.data$battery, .data$method))

cat("\nIncomplete by battery:\n")
print(status |>
  group_by(.data$battery) |>
  summarise(
    total = n(),
    done = sum(.data$valid_complete),
    remaining = .data$total - .data$done,
    .groups = "drop"
  ))
