if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("metric_helpers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("smoke")

workers <- as.integer(Sys.getenv("R2R2_SMOKE_WORKERS", "2"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings(n_boot = as.integer(Sys.getenv("R2R2_SMOKE_N_BOOT", "10")), smoke = TRUE)

smoke_manifest <- make_smoke_manifest()
write_manifest(smoke_manifest, "smoke_manifest.csv")

seed_dat_1 <- generate_legacy13_data(60, 71, "gaussian", "baseline")
seed_dat_2 <- generate_legacy13_data(60, 71, "gaussian", "baseline")
weak_truth <- legacy13_truth("weak")
split_1 <- make_train_validation_split(60, 170071)
split_2 <- make_train_validation_split(60, 170071)
seed_checks <- data.frame(
  check = c(
    "identical_generated_data", "identical_split", "unbalanced_sum_p100", "p200_group_sum",
    "weak_support_preserves_ten_active_components", "weak_mg_is_active_negative",
    "weak_support_matches_nonzero_beta"
  ),
  passed = c(
    identical(digest::digest(seed_dat_1), digest::digest(seed_dat_2)),
    identical(split_1, split_2),
    sum(factorial_group_sizes("unbalanced")) == 100L,
    sum(lengths(make_group_definitions_equal(200L, 10L))) == 200L,
    sum(weak_truth$IsActive) == 10L,
    weak_truth$IsActive[weak_truth$Variable == "mg"] &&
      weak_truth$True_Direction[weak_truth$Variable == "mg"] == "Negative",
    identical(weak_truth$IsActive, abs(weak_truth$True_Beta) > r2r2_truth_tolerance())
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(seed_checks, r2r2_result_file("summaries", "smoke_seed_and_design_checks.csv"))

message("Running smoke jobs with n_boot=", settings$n_boot, "; workers=", workers)
invisible(run_job_table(smoke_manifest, settings = settings, workers = workers, force = force))

fallback_job <- make_base_job(
  battery = "smoke_failure_fallback",
  scenario_id = "smoke_failure_bad_method",
  family = "gaussian",
  n = 50L,
  p = 13L,
  data_seed = 71L,
  method = "bad_method"
)
fallback_res <- run_method_job_resumable(fallback_job, settings = settings, force = TRUE)
fallback_check <- data.frame(
  check = "failure_fallback_creates_failure_metric",
  passed = isFALSE(fallback_res$method_metrics$fit_success[[1]]) &&
    identical(fallback_res$method_metrics$failure_stage[[1]], "fit"),
  stringsAsFactors = FALSE
)
atomic_write_csv(fallback_check, r2r2_result_file("summaries", "smoke_failure_fallback_check.csv"))

all_checks <- dplyr::bind_rows(seed_checks, fallback_check)
results <- aggregate_and_write_outputs(smoke_manifest, include_smoke = TRUE)
method_metrics <- results$method_metrics |> dplyr::filter(grepl("^smoke_", .data$scenario_id))
schema_checks <- data.frame(
  check = c(
    "method_metrics_available",
    "component_metrics_available",
    "sglwqs_direction_schema_available",
    "all_method_wrappers_attempted",
    "gaussian_identity_branch_finite",
    "resume_paths_valid"
  ),
  passed = c(
    nrow(results$method_metrics) > 0L,
    nrow(results$component_metrics) > 0L,
    nrow(results$sglwqs_direction_results) > 0L,
    all(r2r2_methods() %in% unique(method_metrics$method)),
    all(is.finite(results$method_metrics$runtime_sec[results$method_metrics$family == "gaussian" & results$method_metrics$fit_success])),
    all(file.exists(manifest_job_paths(smoke_manifest[smoke_manifest$run_this, , drop = FALSE])))
  ),
  stringsAsFactors = FALSE
)
all_checks <- dplyr::bind_rows(all_checks, schema_checks)
atomic_write_csv(all_checks, r2r2_result_file("summaries", "smoke_test_checks.csv"))

runtime_est <- results$method_metrics |>
  dplyr::filter(grepl("^smoke_", .data$scenario_id), .data$fit_success) |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(
    smoke_jobs = dplyr::n(),
    median_runtime_sec = stats::median(.data$runtime_sec, na.rm = TRUE),
    .groups = "drop"
  )
prod_manifest <- make_production_manifest()
prod_counts <- prod_manifest |> dplyr::filter(.data$run_this) |> dplyr::count(.data$method, name = "production_jobs")
runtime_projection <- dplyr::left_join(prod_counts, runtime_est, by = "method") |>
  dplyr::mutate(
    projected_core_hours = .data$production_jobs * .data$median_runtime_sec / 3600,
    projected_wall_hours_15_workers = .data$projected_core_hours / 15
  )
atomic_write_csv(runtime_projection, r2r2_result_file("summaries", "smoke_runtime_projection.csv"))

smoke_md <- c(
  "# Smoke Test Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Checks",
  capture.output(knitr::kable(all_checks, format = "markdown")),
  "",
  "## Runtime Projection From Smoke Jobs",
  capture.output(knitr::kable(runtime_projection, format = "markdown", digits = 3)),
  "",
  "## Notes",
  "- Smoke jobs use reduced bootstrap count and are not interpreted as production results.",
  "- Failure fallback is checked using a deliberately invalid method label; the corresponding result is stored as a failure RDS."
)
writeLines(smoke_md, r2r2_result_file("summaries", "smoke_test_report.md"))

if (!all(all_checks$passed)) {
  stop("One or more smoke checks failed. See results/summaries/smoke_test_checks.csv", call. = FALSE)
}

message("Smoke test passed.")
