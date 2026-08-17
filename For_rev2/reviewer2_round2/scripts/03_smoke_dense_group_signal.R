if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("dense_group_signal_generator.R")
source_r2r2("dense_group_signal_helpers.R")
install_dense_group_signal_dispatch(.GlobalEnv)

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("dense_group_signal_smoke")

selected_families <- dense_group_signal_selected_families()
manifest <- make_dense_group_signal_smoke_manifest()
manifest <- manifest[manifest$family %in% selected_families, , drop = FALSE]
write_manifest(manifest, "dense_group_signal_smoke_manifest.csv")
seed <- unique(manifest$data_seed)[[1L]]

dense_gaussian_a <- generate_dense_group_signal_data(500L, seed, "gaussian")
dense_gaussian_b <- generate_dense_group_signal_data(500L, seed, "gaussian")
sparse_gaussian <- generate_factorial_data(
  500L, seed, "gaussian", "unbalanced", "weak_heterogeneous", 0.45, 0.05
)
dense_binomial <- generate_dense_group_signal_data(500L, seed, "binomial")

vars <- unlist(attr(dense_gaussian_a, "groups"), use.names = FALSE)
covars <- legacy13_covariates()
truth <- attr(dense_gaussian_a, "truth")
signal_spec <- attr(dense_gaussian_a, "signal_spec")

assert_or_stop(identical(dense_gaussian_a, dense_gaussian_b), "Dense generator is not seed-reproducible.")
assert_or_stop(identical(dense_gaussian_a[c(vars, covars)], sparse_gaussian[c(vars, covars)]),
  "Dense and sparse matched datasets do not share X and covariates."
)
assert_or_stop(isTRUE(all.equal(
  dense_gaussian_a$Y - attr(dense_gaussian_a, "eta_linear"),
  sparse_gaussian$Y - attr(sparse_gaussian, "eta_linear"),
  tolerance = 1e-12
)), "Dense and sparse Gaussian datasets do not share residual draws.")
assert_or_stop(abs(stats::sd(dense_gaussian_a$Y - attr(dense_gaussian_a, "eta_linear")) - 1) < 0.15,
  "Dense Gaussian residual SD is inconsistent with 1."
)
assert_or_stop(length(unique(dense_binomial$Y)) == 2L, "Dense binomial outcome does not contain both classes.")
assert_or_stop(sum(lengths(attr(dense_gaussian_a, "groups"))) == 100L, "Dense group sizes do not sum to 100.")
assert_or_stop(sum(truth$IsActive) == 30L, "Dense truth does not contain 30 active components.")
assert_or_stop(all(truth$True_Beta[!truth$IsActiveGroup] == 0), "A null group contains nonzero truth.")
assert_or_stop(all(abs(signal_spec$Target_Signal_Variance - signal_spec$Dense_Signal_Variance) < 1e-12),
  "Sparse and dense group signal variances differ."
)

job <- manifest[manifest$method == "qgcomp" & manifest$family == "gaussian", , drop = FALSE][1, ]
split_a <- make_train_validation_split(job$n, job$split_seed, 0.6)
split_b <- make_train_validation_split(job$n, job$split_seed, 0.6)
assert_or_stop(identical(split_a, split_b), "Split seed was not reproducible.")
assert_or_stop(!length(intersect(split_a$train, split_a$validation)), "Train and validation sets overlap.")
assert_or_stop(length(union(split_a$train, split_a$validation)) == job$n, "Split is not exhaustive.")

settings <- r2r2_settings(n_boot = as.integer(Sys.getenv("R2R2_SMOKE_N_BOOT", "10")), smoke = TRUE)
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
results <- run_dense_group_signal_jobs(
  manifest,
  settings = settings,
  workers = as.integer(Sys.getenv("R2R2_SMOKE_WORKERS", "1")),
  force = force
)

assert_or_stop(length(results) == nrow(manifest), "Smoke result count does not match the manifest.")
fit_rows <- dplyr::bind_rows(lapply(results, `[[`, "method_metrics"))
assert_or_stop(all(fit_rows$fit_success), "At least one dense group-signal smoke fit failed.")
assert_or_stop(all(vapply(results, function(x) identical(x$schema_version, "reviewer2_round2_v1"), logical(1))),
  "Dense smoke schema version mismatch."
)
component_columns <- lapply(results, function(x) sort(names(x$component_metrics)))
assert_or_stop(length(unique(vapply(component_columns, paste, collapse = "|", character(1)))) == 1L,
  "Dense component schemas differ across methods."
)

paths <- dense_group_signal_manifest_paths(manifest)
mtime_before <- file.info(paths)$mtime
invisible(run_dense_group_signal_jobs(manifest, settings = settings, workers = 1L, force = FALSE))
mtime_after <- file.info(paths)$mtime
assert_or_stop(identical(mtime_before, mtime_after), "Resume check rewrote complete dense smoke outputs.")

status <- summarize_dense_group_signal_status(manifest, "dense_group_signal_smoke_completion_status.csv")
checks <- data.frame(
  check = c(
    "data_seed_reproducibility", "matched_X_and_covariates", "matched_gaussian_residuals",
    "gaussian_residual_sd", "binary_both_classes", "group_size_sum", "active_component_count",
    "null_groups_zero", "group_signal_variance_match", "split_reproducibility",
    "split_disjoint_exhaustive", "all_method_wrappers", "schema_compatibility",
    "resume_without_rewrite"
  ),
  passed = TRUE,
  stringsAsFactors = FALSE
)
atomic_write_csv(checks, r2r2_result_file("summaries", "dense_group_signal_smoke_checks.csv"))
atomic_write_csv(signal_spec, r2r2_result_file("summaries", "dense_group_signal_specification.csv"))

report <- c(
  "# Dense Group-Signal Smoke Test",
  "",
  sprintf("- Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("- Families: %s", paste(selected_families, collapse = ", ")),
  sprintf("- Jobs: %d/%d valid", sum(status$valid_complete), nrow(status)),
  sprintf("- Fit successes: %d/%d", sum(fit_rows$fit_success), nrow(fit_rows)),
  sprintf("- Bootstrap iterations per WQS job: %d", settings$n_boot),
  sprintf("- Active components: %d/100", sum(truth$IsActive)),
  sprintf("- SGL-WQS source: `%s`", settings$sglwqs_source),
  sprintf("- SGL-WQS commit: `%s`", git_commit_for_path(settings$sglwqs_source)),
  "",
  "All matched-data, signal-variance, truth, split, wrapper, schema, and resume checks passed."
)
writeLines(report, r2r2_result_file("summaries", "DENSE_GROUP_SIGNAL_SMOKE_REPORT.md"))
cat(paste(report, collapse = "\n"), "\n")
