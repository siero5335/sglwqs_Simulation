if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("gaussian_primary_benchmark_generators.R")
source_r2r2("gaussian_primary_benchmark_helpers.R")
install_gaussian_primary_benchmark_dispatch(.GlobalEnv)

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("gaussian_primary_benchmark_smoke")

manifest <- make_gaussian_primary_benchmark_smoke_manifest()
write_manifest(manifest, "gaussian_primary_benchmark_smoke_manifest.csv")

cor_job <- manifest[grepl("correlation", manifest$battery) & manifest$method == "qgcomp", , drop = FALSE][1, ]
active_job <- manifest[grepl("active_component", manifest$battery) & manifest$method == "qgcomp", , drop = FALSE][1, ]
cor_a <- generate_gaussian_primary_benchmark_job_data(cor_job)
cor_b <- generate_gaussian_primary_benchmark_job_data(cor_job)
active <- generate_gaussian_primary_benchmark_job_data(active_job)

assert_or_stop(identical(cor_a, cor_b), "Identical seed did not reproduce the correlation dataset.")
assert_or_stop(identical(attr(cor_a, "truth"), attr(cor_b, "truth")), "Identical seed did not reproduce correlation truth metadata.")
assert_or_stop(all(is.finite(cor_a$Y)) && all(is.finite(active$Y)), "Gaussian outcome contains non-finite values.")
assert_or_stop(abs(stats::sd(cor_a$Y - attr(cor_a, "eta_linear")) - 1) < 0.15, "Correlation residual SD is inconsistent with 1.")
assert_or_stop(abs(stats::sd(active$Y - attr(active, "eta_linear")) - 1) < 0.15, "Active-component residual SD is inconsistent with 1.")
assert_or_stop(identical(attr(cor_a, "beta"), gaussian_correlation_benchmark_beta()), "Correlation truth differs from the GitHub design.")
assert_or_stop(identical(attr(active, "beta"), gaussian_active_component_beta()), "Active-component truth differs from the supplied design.")
assert_or_stop(sum(lengths(attr(cor_a, "groups"))) == 13L, "Correlation group sizes do not sum to 13.")
assert_or_stop(sum(lengths(attr(active, "groups"))) == 13L, "Active-component group sizes do not sum to 13.")

split_a <- make_train_validation_split(nrow(cor_a), cor_job$split_seed, 0.6)
split_b <- make_train_validation_split(nrow(cor_a), cor_job$split_seed, 0.6)
assert_or_stop(identical(split_a, split_b), "Split seed was not reproducible.")
assert_or_stop(!length(intersect(split_a$train, split_a$validation)), "Train and validation sets overlap.")
assert_or_stop(length(union(split_a$train, split_a$validation)) == nrow(cor_a), "Train and validation sets are not exhaustive.")

settings <- r2r2_settings(n_boot = as.integer(Sys.getenv("R2R2_SMOKE_N_BOOT", "10")), smoke = TRUE)
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
results <- run_gaussian_primary_benchmark_jobs(
  manifest,
  settings = settings,
  workers = as.integer(Sys.getenv("R2R2_SMOKE_WORKERS", "1")),
  force = force
)

assert_or_stop(length(results) == nrow(manifest), "Smoke result count does not match the manifest.")
fit_rows <- dplyr::bind_rows(lapply(results, `[[`, "method_metrics"))
assert_or_stop(all(fit_rows$fit_success), "At least one Gaussian primary benchmark smoke fit failed.")
assert_or_stop(all(fit_rows$family == "gaussian"), "A smoke output used a non-Gaussian family.")
assert_or_stop(all(vapply(results, function(x) identical(x$schema_version, "reviewer2_round2_v1"), logical(1))), "Smoke schema version mismatch.")

component_columns <- lapply(results, function(x) sort(names(x$component_metrics)))
assert_or_stop(length(unique(vapply(component_columns, paste, collapse = "|", character(1)))) == 1L, "Component schemas differ across methods.")

paths <- gaussian_primary_manifest_paths(manifest)
mtime_before <- file.info(paths)$mtime
invisible(run_gaussian_primary_benchmark_jobs(manifest, settings = settings, workers = 1L, force = FALSE))
mtime_after <- file.info(paths)$mtime
assert_or_stop(identical(mtime_before, mtime_after), "Resume check rewrote complete smoke outputs.")

status <- summarize_gaussian_primary_benchmark_status(
  manifest,
  "gaussian_primary_benchmark_smoke_completion_status.csv"
)
checks <- data.frame(
  check = c(
    "data_seed_reproducibility", "gaussian_residual_sd", "published_truth_match",
    "group_size_sum", "split_reproducibility", "split_disjoint_exhaustive",
    "all_method_wrappers", "schema_compatibility", "resume_without_rewrite"
  ),
  passed = TRUE,
  stringsAsFactors = FALSE
)
atomic_write_csv(checks, r2r2_result_file("summaries", "gaussian_primary_benchmark_smoke_checks.csv"))

report <- c(
  "# Gaussian Primary Benchmark Smoke Test",
  "",
  sprintf("- Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- Jobs: %d/%d valid", sum(status$valid_complete), nrow(status)),
  sprintf("- Fit successes: %d/%d", sum(fit_rows$fit_success), nrow(fit_rows)),
  sprintf("- Bootstrap iterations per WQS job: %d", settings$n_boot),
  sprintf("- SGL-WQS source: `%s`", settings$sglwqs_source),
  sprintf("- SGL-WQS commit: `%s`", git_commit_for_path(settings$sglwqs_source)),
  "",
  "All generator, truth, split, wrapper, schema, and resume checks passed."
)
writeLines(report, r2r2_result_file("summaries", "GAUSSIAN_PRIMARY_BENCHMARK_SMOKE_REPORT.md"))
cat(paste(report, collapse = "\n"), "\n")
