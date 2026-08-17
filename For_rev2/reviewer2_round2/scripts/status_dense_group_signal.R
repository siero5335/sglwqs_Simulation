source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("dense_group_signal_helpers.R")

selected_families <- dense_group_signal_selected_families()
manifest <- make_dense_group_signal_manifest()
manifest <- manifest[manifest$family %in% selected_families, , drop = FALSE]
paths <- dense_group_signal_manifest_paths(manifest)
valid <- vapply(paths, is_complete_job_file, logical(1))
status <- manifest |>
  dplyr::mutate(valid_complete = valid) |>
  dplyr::group_by(.data$family, .data$n, .data$method) |>
  dplyr::summarise(done = sum(.data$valid_complete), expected = dplyr::n(), .groups = "drop")

cat(sprintf(
  "Battery H (%s) valid jobs: %d/%d; remaining: %d\n",
  paste(selected_families, collapse = ", "), sum(valid), length(valid), sum(!valid)
))
print(status, n = Inf)
