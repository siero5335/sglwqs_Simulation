source("reviewer2_round2/R/design_helpers.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")

unset_namespace <- Sys.getenv("R2R2_RESULT_NAMESPACE", unset = "")
if (nzchar(unset_namespace)) {
  stop("Production validation must run with R2R2_RESULT_NAMESPACE unset.", call. = FALSE)
}

battery_key <- Sys.getenv("R2R2_VALIDATE_BATTERY", unset = "")
spec <- switch(battery_key,
  C = list(name = "C_gaussian_split_stability", expected = 400L),
  D = list(name = "D_unbalanced_effect_factorial", expected = 1440L),
  stop("Set R2R2_VALIDATE_BATTERY to C or D.", call. = FALSE)
)

manifest <- make_production_manifest() |>
  dplyr::filter(.data$battery == spec$name, .data$run_this)
paths <- manifest_job_paths(manifest)
valid <- vapply(paths, is_complete_job_file, logical(1))
results <- lapply(paths[valid], readRDS)
method_metrics <- dplyr::bind_rows(lapply(results, `[[`, "method_metrics"))
data_diagnostics <- dplyr::bind_rows(lapply(results, `[[`, "data_diagnostics"))
seed_metadata <- dplyr::bind_rows(lapply(results, `[[`, "seed_metadata"))

prevalence_check <- if (identical(battery_key, "D")) {
  binomial_rows <- method_metrics$family == "binomial"
  all(c(
    "binomial_target_prevalence",
    "binomial_expected_prevalence",
    "y_positive_rate",
    "y_positive_count"
  ) %in% names(data_diagnostics)) &&
    all(abs(data_diagnostics$binomial_target_prevalence[binomial_rows] - 0.20) < 1e-12) &&
    all(abs(
      data_diagnostics$binomial_expected_prevalence[binomial_rows] -
        data_diagnostics$binomial_target_prevalence[binomial_rows]
    ) < 1e-10) &&
    all(data_diagnostics$y_positive_rate[binomial_rows] >= 0.10) &&
    all(data_diagnostics$y_positive_rate[binomial_rows] <= 0.30)
} else {
  TRUE
}

boot_methods <- method_metrics$method %in% c("SGL-WQS", "gWQS", "groupWQS")
boot_setting_ok <- nrow(method_metrics) == spec$expected &&
  all(method_metrics$bootstrap_attempted[boot_methods] == 200L, na.rm = FALSE) &&
  all(method_metrics$bootstrap_attempted[!boot_methods] == 0L, na.rm = FALSE)

hash_check <- if (identical(battery_key, "C")) {
  hash_table <- data.frame(
    scenario_id = manifest$scenario_id[valid],
    data_hash = data_diagnostics$data_hash,
    split_hash = data_diagnostics$split_hash,
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$scenario_id) |>
    dplyr::summarise(
      data_hashes = dplyr::n_distinct(.data$data_hash),
      split_hashes = dplyr::n_distinct(.data$split_hash),
      .groups = "drop"
    )
  nrow(hash_table) == 4L && all(hash_table$data_hashes == 1L) && all(hash_table$split_hashes == 100L)
} else {
  hash_table <- data.frame(
    scenario_id = manifest$scenario_id[valid],
    data_seed = manifest$data_seed[valid],
    data_hash = data_diagnostics$data_hash,
    split_hash = data_diagnostics$split_hash,
    stringsAsFactors = FALSE
  ) |>
    dplyr::group_by(.data$scenario_id, .data$data_seed) |>
    dplyr::summarise(
      data_hashes = dplyr::n_distinct(.data$data_hash),
      split_hashes = dplyr::n_distinct(.data$split_hash),
      methods = dplyr::n(),
      .groups = "drop"
    )
  nrow(hash_table) == 360L && all(hash_table$data_hashes == 1L) &&
    all(hash_table$split_hashes == 1L) && all(hash_table$methods == 4L)
}

checks <- data.frame(
  check = c(
    "manifest_expected_jobs",
    "all_atomic_outputs_valid",
    "one_metric_row_per_job",
    "production_bootstrap_setting_200",
    "seed_metadata_saved",
    "data_and_split_hash_design",
    "factorial_binomial_prevalence_target_20pct"
  ),
  passed = c(
    nrow(manifest) == spec$expected,
    sum(valid) == spec$expected,
    nrow(method_metrics) == spec$expected,
    boot_setting_ok,
    nrow(seed_metadata) == spec$expected &&
      all(is.finite(seed_metadata$data_seed)) && all(is.finite(seed_metadata$split_seed)),
    hash_check,
    prevalence_check
  ),
  detail = c(
    sprintf("%d/%d", nrow(manifest), spec$expected),
    sprintf("%d/%d", sum(valid), spec$expected),
    sprintf("%d/%d", nrow(method_metrics), spec$expected),
    "n_boot=200 for SGL-WQS/gWQS/groupWQS; qgcomp=0",
    sprintf("%d rows", nrow(seed_metadata)),
    if (identical(battery_key, "C")) "one fixed dataset and 100 distinct splits per scenario" else "identical data/split across four methods per scenario-seed",
    if (identical(battery_key, "D")) "expected prevalence calibrated to 0.20; realized prevalence required in [0.10, 0.30]" else "not applicable to Battery C"
  ),
  stringsAsFactors = FALSE
)

fit_summary <- method_metrics |>
  dplyr::group_by(.data$scenario_id, .data$family, .data$method) |>
  dplyr::summarise(
    completed = dplyr::n(),
    fit_success = sum(.data$fit_success),
    recorded_failure = sum(!.data$fit_success),
    runtime_median_sec = stats::median(.data$runtime_sec, na.rm = TRUE),
    runtime_max_sec = safe_max(.data$runtime_sec),
    .groups = "drop"
  )

failure_summary <- method_metrics |>
  dplyr::filter(!.data$fit_success) |>
  dplyr::count(
    .data$scenario_id, .data$family, .data$method,
    .data$failure_stage, .data$error_class, .data$error_message,
    name = "jobs"
  )

atomic_write_csv(checks, r2r2_result_file("summaries", paste0("delegated_", battery_key, "_production_gate_checks.csv")))
atomic_write_csv(fit_summary, r2r2_result_file("summaries", paste0("delegated_", battery_key, "_production_fit_summary.csv")))
atomic_write_csv(failure_summary, r2r2_result_file("summaries", paste0("delegated_", battery_key, "_production_failure_summary.csv")))
atomic_write_csv(hash_table, r2r2_result_file("summaries", paste0("delegated_", battery_key, "_production_hash_checks.csv")))

if (!all(checks$passed)) {
  stop("Delegated Battery ", battery_key, " production gate failed.", call. = FALSE)
}

message(
  "Delegated Battery ", battery_key, " production gate passed: ",
  sum(method_metrics$fit_success), " fit successes; ",
  sum(!method_metrics$fit_success), " recorded failures."
)
