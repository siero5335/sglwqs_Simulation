bind_direction_results <- function(results) {
  do.call(rbind, lapply(results, `[[`, "direction_results"))
}

bind_diagnostics <- function(results) {
  do.call(rbind, lapply(results, `[[`, "diagnostics"))
}

binom_ci <- function(x, n, conf.level = 0.95) {
  if (!is.finite(x) || !is.finite(n) || n <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  stats::binom.test(as.integer(x), as.integer(n), conf.level = conf.level)$conf.int
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

summarize_calibration <- function(direction_results) {
  if (nrow(direction_results) == 0L) {
    return(data.frame())
  }
  dplyr::group_by(
    direction_results,
    profile, job_type, scenario_family, scenario_id, scenario_label,
    n, outcome_family, effect, within_r, between_r, group, direction, is_true_null
  ) |>
    dplyr::summarise(
      attempted = dplyr::n(),
      fit_success_n = sum(fit_success, na.rm = TRUE),
      fit_failure_rate = 1 - fit_success_n / attempted,
      retained_n = sum(retained, na.rm = TRUE),
      retention_rate = retained_n / fit_success_n,
      rejection_count_unconditional = sum(fit_success & is_true_null & rejected_0_05, na.rm = TRUE),
      rejection_denom_unconditional = sum(fit_success & is_true_null, na.rm = TRUE),
      rejection_rate_unconditional = rejection_count_unconditional / rejection_denom_unconditional,
      rejection_count_retained = sum(retained & is_true_null & rejected_0_05, na.rm = TRUE),
      rejection_denom_retained = sum(retained & is_true_null, na.rm = TRUE),
      rejection_rate_retained = rejection_count_retained / rejection_denom_retained,
      estimate_mean = mean(estimate, na.rm = TRUE),
      estimate_median = stats::median(estimate, na.rm = TRUE),
      empirical_sd = stats::sd(estimate, na.rm = TRUE),
      model_se_mean = mean(std_error, na.rm = TRUE),
      se_to_empirical_sd = model_se_mean / empirical_sd,
      null_coverage_count = sum(retained & is_true_null & conf_low <= 0 & conf_high >= 0, na.rm = TRUE),
      null_coverage_denom = sum(retained & is_true_null, na.rm = TRUE),
      null_coverage = null_coverage_count / null_coverage_denom,
      median_runtime_sec = stats::median(runtime_sec, na.rm = TRUE),
      median_boot_success_rate = stats::median(boot_success_rate, na.rm = TRUE),
      lambda_boundary_rate = mean(selected_lambda_at_path_boundary, na.rm = TRUE),
      all_zero_rate = mean(all_zero_exposure, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      rejection_unconditional_lcl = binom_ci(rejection_count_unconditional, rejection_denom_unconditional)[1],
      rejection_unconditional_ucl = binom_ci(rejection_count_unconditional, rejection_denom_unconditional)[2],
      rejection_retained_lcl = binom_ci(rejection_count_retained, rejection_denom_retained)[1],
      rejection_retained_ucl = binom_ci(rejection_count_retained, rejection_denom_retained)[2],
      null_coverage_lcl = binom_ci(null_coverage_count, null_coverage_denom)[1],
      null_coverage_ucl = binom_ci(null_coverage_count, null_coverage_denom)[2]
    ) |>
    dplyr::ungroup()
}

summarize_split_stability <- function(direction_results) {
  split <- direction_results[direction_results$job_type == "split_stability", , drop = FALSE]
  if (nrow(split) == 0L) {
    return(data.frame())
  }
  dplyr::group_by(
    split,
    profile, scenario_family, scenario_id, scenario_label,
    n, outcome_family, effect, within_r, between_r, group, direction, is_true_null
  ) |>
    dplyr::summarise(
      split_count = dplyr::n(),
      estimable_split_count = sum(retained & is.finite(p_value), na.rm = TRUE),
      fit_success_count = sum(fit_success, na.rm = TRUE),
      p_value_mean = mean(p_value, na.rm = TRUE),
      p_value_sd = stats::sd(p_value, na.rm = TRUE),
      p_value_cv = p_value_sd / p_value_mean,
      p_value_iqr = stats::IQR(p_value, na.rm = TRUE),
      p_value_min = safe_min(p_value),
      p_value_max = safe_max(p_value),
      estimate_mean = mean(estimate, na.rm = TRUE),
      estimate_sd = stats::sd(estimate, na.rm = TRUE),
      positive_sign_rate = mean(sign == "positive", na.rm = TRUE),
      negative_sign_rate = mean(sign == "negative", na.rm = TRUE),
      rejection_rate = mean(rejected_0_05, na.rm = TRUE),
      retention_rate = mean(retained, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_runtime <- function(diagnostics, profile) {
  if (nrow(diagnostics) == 0L) {
    return(data.frame())
  }
  observed <- dplyr::group_by(diagnostics, profile, job_type) |>
    dplyr::summarise(
      observed_jobs = dplyr::n(),
      successful_jobs = sum(fit_success, na.rm = TRUE),
      median_runtime_sec = stats::median(runtime_sec, na.rm = TRUE),
      p90_runtime_sec = as.numeric(stats::quantile(runtime_sec, 0.9, na.rm = TRUE)),
      total_observed_core_sec = sum(runtime_sec, na.rm = TRUE),
      .groups = "drop"
    )

  production_counts <- data.frame(
    job_type = c("replicate", "split_stability"),
    production_jobs = c(
      4 * 500 + 4 * 300 + 2 * 200,
      4 * 100
    ),
    stringsAsFactors = FALSE
  )

  dplyr::left_join(observed, production_counts, by = "job_type") |>
    dplyr::mutate(
      projected_core_hours = median_runtime_sec * production_jobs / 3600,
      projected_wall_hours_8_workers = projected_core_hours / 8
    )
}

write_summary_tables <- function(results, profile) {
  ensure_vc_dirs()
  directions <- bind_direction_results(results)
  diagnostics <- bind_diagnostics(results)
  calibration <- summarize_calibration(directions)
  split <- summarize_split_stability(directions)
  runtime <- summarize_runtime(diagnostics, profile)

  prefix <- profile$name
  atomic_write_csv(directions, vc_file("output", "tables", paste0("validation_replicate_results_", prefix, ".csv")))
  atomic_write_csv(diagnostics, vc_file("output", "tables", paste0("validation_diagnostics_", prefix, ".csv")))
  atomic_write_csv(calibration, vc_file("output", "tables", paste0("validation_calibration_summary_", prefix, ".csv")))
  atomic_write_csv(split, vc_file("output", "tables", paste0("validation_split_stability_", prefix, ".csv")))
  atomic_write_csv(runtime, vc_file("output", "tables", paste0("validation_runtime_projection_", prefix, ".csv")))

  list(
    directions = directions,
    diagnostics = diagnostics,
    calibration = calibration,
    split_stability = split,
    runtime = runtime
  )
}
