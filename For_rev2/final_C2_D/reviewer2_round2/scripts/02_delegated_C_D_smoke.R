source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("metric_helpers.R")

# Keep delegated smoke artifacts physically separate from production raw files.
Sys.setenv(R2R2_RESULT_NAMESPACE = "delegated_C_D_smoke")
r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("delegated_C_D_smoke")

workers <- as.integer(Sys.getenv("R2R2_SMOKE_WORKERS", "8"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings(
  n_boot = as.integer(Sys.getenv("R2R2_SMOKE_N_BOOT", "10")),
  smoke = TRUE
)

paired_seed <- 42L
paired_jobs <- purrr::map_dfr(r2r2_methods(), function(method) {
  make_base_job(
    battery = "D_unbalanced_effect_factorial",
    scenario_id = "delegated_smoke_D_binomial_balanced_heterogeneous",
    family = "binomial",
    n = 180L,
    p = 100L,
    data_seed = paired_seed,
    method = method,
    signal_profile = "heterogeneous",
    effect_profile = "heterogeneous",
    group_structure = "balanced",
    within_rho = 0.45,
    cross_rho = 0.05
  )
})

smoke_manifest <- dplyr::bind_rows(
  make_base_job(
    battery = "C_gaussian_split_stability",
    scenario_id = "delegated_smoke_C_gaussian_global_null_n500",
    family = "gaussian",
    n = 500L,
    p = 13L,
    data_seed = 991001L,
    method = "SGL-WQS",
    signal_profile = "global_null",
    effect_profile = "global_null",
    group_structure = "legacy13",
    within_rho = 0.20,
    cross_rho = 0,
    split_replicate = 1L,
    split_seed = 1061002L,
    method_seed = 1081002L
  ),
  paired_jobs,
  make_base_job(
    "D_unbalanced_effect_factorial",
    "delegated_smoke_D_binomial_unbalanced_heterogeneous",
    "binomial", 180L, 100L, paired_seed, "qgcomp",
    signal_profile = "heterogeneous", effect_profile = "heterogeneous",
    group_structure = "unbalanced", within_rho = 0.45, cross_rho = 0.05
  ),
  make_base_job(
    "D_unbalanced_effect_factorial",
    "delegated_smoke_D_gaussian_balanced_heterogeneous",
    "gaussian", 180L, 100L, paired_seed, "qgcomp",
    signal_profile = "heterogeneous", effect_profile = "heterogeneous",
    group_structure = "balanced", within_rho = 0.45, cross_rho = 0.05
  ),
  make_base_job(
    "D_unbalanced_effect_factorial",
    "delegated_smoke_D_gaussian_unbalanced_heterogeneous",
    "gaussian", 180L, 100L, paired_seed, "SGL-WQS",
    signal_profile = "heterogeneous", effect_profile = "heterogeneous",
    group_structure = "unbalanced", within_rho = 0.45, cross_rho = 0.05
  )
) |>
  dplyr::mutate(job_id = paste0("delegated_smoke_", .data$job_id))

atomic_write_csv(
  smoke_manifest,
  r2r2_result_file("summaries", "delegated_C_D_smoke_manifest.csv")
)

# Deterministic design checks independent of model convergence.
truth_unbalanced <- make_factorial_truth("unbalanced", "heterogeneous")
truth_weak <- make_factorial_truth("unbalanced", "weak_heterogeneous")
null_groups <- setdiff(names(truth_unbalanced$groups), c("G01", "G05", "G10"))
null_vars <- unlist(truth_unbalanced$groups[null_groups], use.names = FALSE)

cor_groups <- make_group_definitions_from_sizes(factorial_group_sizes("unbalanced"))
cor_x <- generate_latent_exposure_matrix(5000L, cor_groups, seed = 8675309L, within_rho = 0.45, cross_rho = 0.05)
cor_mat <- stats::cor(log(cor_x))
group_id <- rep(seq_along(cor_groups), lengths(cor_groups))
upper <- upper.tri(cor_mat)
within_mean <- mean(cor_mat[upper & outer(group_id, group_id, `==`)])
cross_mean <- mean(cor_mat[upper & outer(group_id, group_id, `!=`)])

gaussian_dat <- generate_factorial_data(
  2000L, 777L, "gaussian", "unbalanced", "heterogeneous",
  within_rho = 0.45, cross_rho = 0.05
)
gaussian_residual <- gaussian_dat$Y - attr(gaussian_dat, "eta_linear")
binomial_dat <- generate_factorial_data(
  2000L, 778L, "binomial", "unbalanced", "heterogeneous",
  within_rho = 0.45, cross_rho = 0.05
)
binomial_target <- attr(binomial_dat, "binomial_target_prevalence")
binomial_expected <- attr(binomial_dat, "binomial_expected_prevalence")
binomial_realized <- mean(binomial_dat$Y)

paired_data_hashes <- vapply(
  split(paired_jobs, seq_len(nrow(paired_jobs))),
  function(job) digest::digest(generate_job_data(job), algo = "sha256"),
  character(1)
)
paired_split_hashes <- vapply(
  split(paired_jobs, seq_len(nrow(paired_jobs))),
  function(job) {
    dat <- generate_job_data(job)
    split <- make_train_validation_split(
      nrow(dat), job$split_seed, settings$train_prop,
      y = dat$Y, family = job$family
    )
    digest::digest(split, algo = "sha256")
  },
  character(1)
)

r2r2_load_local_sglwqs()
local_source_version <- as.character(utils::packageVersion("sglwqs"))
design_checks <- data.frame(
  check = c(
    "local_source_is_requested_version",
    "unbalanced_group_sizes_exact",
    "active_groups_have_three_components",
    "null_groups_are_exactly_zero",
    "weak_profile_is_half_heterogeneous",
    "direction_plan_positive_negative_mixed",
    "target_within_correlation",
    "target_cross_correlation",
    "gaussian_identity_residual_sd_one",
    "binomial_expected_prevalence_target_20pct",
    "binomial_realized_prevalence_near_target",
    "identical_data_across_methods",
    "identical_split_across_methods",
    "smoke_namespace_not_production_raw"
  ),
  passed = c(
    identical(local_source_version, "0.8.13.9001"),
    identical(factorial_group_sizes("unbalanced"), c(4L, 6L, 8L, 8L, 10L, 10L, 12L, 12L, 14L, 16L)),
    all(vapply(truth_unbalanced$groups[c("G01", "G05", "G10")], function(g) sum(truth_unbalanced$beta[g] != 0), integer(1)) == 3L),
    all(truth_unbalanced$beta[null_vars] == 0),
    isTRUE(all.equal(abs(truth_weak$beta), 0.5 * abs(truth_unbalanced$beta))),
    all(truth_unbalanced$beta[truth_unbalanced$groups$G01] >= 0) &&
      all(truth_unbalanced$beta[truth_unbalanced$groups$G05] <= 0) &&
      any(truth_unbalanced$beta[truth_unbalanced$groups$G10] > 0) &&
      any(truth_unbalanced$beta[truth_unbalanced$groups$G10] < 0),
    abs(within_mean - 0.45) < 0.03,
    abs(cross_mean - 0.05) < 0.03,
    abs(stats::sd(gaussian_residual) - 1) < 0.06,
    identical(binomial_target, factorial_binomial_target_prevalence()) &&
      abs(binomial_expected - binomial_target) < 1e-10,
    abs(binomial_realized - binomial_target) < 0.05,
    length(unique(paired_data_hashes)) == 1L,
    length(unique(paired_split_hashes)) == 1L,
    all(startsWith(
      manifest_job_paths(smoke_manifest),
      file.path(r2r2_root(), "results", "delegated_C_D_smoke", "raw")
    ))
  ),
  detail = c(
    local_source_version,
    paste(factorial_group_sizes("unbalanced"), collapse = ","),
    "G01/G05/G10",
    as.character(max(abs(truth_unbalanced$beta[null_vars]))),
    "0.5 multiplier",
    "positive/negative/mixed",
    sprintf("%.4f", within_mean),
    sprintf("%.4f", cross_mean),
    sprintf("%.4f", stats::sd(gaussian_residual)),
    sprintf("target=%.4f; expected=%.4f", binomial_target, binomial_expected),
    sprintf("target=%.4f; realized=%.4f", binomial_target, binomial_realized),
    paste(unique(paired_data_hashes), collapse = ","),
    paste(unique(paired_split_hashes), collapse = ","),
    r2r2_result_root()
  ),
  stringsAsFactors = FALSE
)
atomic_write_csv(design_checks, r2r2_result_file("summaries", "delegated_C_D_design_checks.csv"))

message("Running delegated C/D smoke jobs with n_boot=", settings$n_boot, "; workers=", workers)
invisible(run_job_table(smoke_manifest, settings = settings, workers = workers, force = force))

# A second pass must skip all completed atomic outputs without changing them.
smoke_paths <- manifest_job_paths(smoke_manifest)
mtime_before <- file.info(smoke_paths)$mtime
invisible(run_job_table(smoke_manifest, settings = settings, workers = workers, force = FALSE))
mtime_after <- file.info(smoke_paths)$mtime
resume_ok <- all(!is.na(mtime_before)) && identical(mtime_before, mtime_after)

fallback_job <- make_base_job(
  battery = "D_unbalanced_effect_factorial",
  scenario_id = "delegated_smoke_failure_bad_method",
  family = "gaussian",
  n = 80L,
  p = 100L,
  data_seed = 71L,
  method = "bad_method",
  signal_profile = "heterogeneous",
  effect_profile = "heterogeneous",
  group_structure = "unbalanced",
  within_rho = 0.45,
  cross_rho = 0.05
)
fallback_res <- run_method_job_resumable(fallback_job, settings = settings, force = TRUE)
fallback_ok <- isFALSE(fallback_res$method_metrics$fit_success[[1]]) &&
  identical(fallback_res$method_metrics$failure_stage[[1]], "dispatch")

smoke_results <- read_all_job_results()
method_metrics <- bind_result_slot(smoke_results, "method_metrics") |>
  dplyr::filter(grepl("^delegated_smoke_", .data$scenario_id))
component_metrics <- bind_result_slot(smoke_results, "component_metrics")
sglwqs_direction_results <- bind_result_slot(smoke_results, "sglwqs_direction_results")
atomic_write_csv(method_metrics, r2r2_result_file("tables", "delegated_C_D_atomic_method_metrics.csv"))
atomic_write_csv(component_metrics, r2r2_result_file("tables", "delegated_C_D_atomic_component_metrics.csv"))
atomic_write_csv(sglwqs_direction_results, r2r2_result_file("tables", "delegated_C_D_sglwqs_direction_results.csv"))
output_checks <- data.frame(
  check = c(
    "all_atomic_outputs_valid",
    "all_four_wrappers_attempted",
    "battery_C_sglwqs_completed",
    "sglwqs_p100_completed",
    "qgcomp_p100_completed",
    "resume_skips_completed_outputs",
    "failure_fallback_recorded"
  ),
  passed = c(
    all(vapply(smoke_paths, is_complete_job_file, logical(1))),
    all(r2r2_methods() %in% unique(method_metrics$method)),
    any(method_metrics$battery == "C_gaussian_split_stability" & method_metrics$method == "SGL-WQS"),
    any(method_metrics$p == 100L & method_metrics$method == "SGL-WQS"),
    any(method_metrics$p == 100L & method_metrics$method == "qgcomp"),
    resume_ok,
    fallback_ok
  ),
  detail = c(
    sprintf("%d/%d", sum(vapply(smoke_paths, is_complete_job_file, logical(1))), length(smoke_paths)),
    paste(r2r2_methods()[r2r2_methods() %in% unique(method_metrics$method)], collapse = ", "),
    "C Gaussian global-null n=500",
    "D p=100",
    "D p=100",
    as.character(resume_ok),
    paste(fallback_res$method_metrics$failure_stage, fallback_res$method_metrics$error_message, sep = ": ")
  ),
  stringsAsFactors = FALSE
)
all_checks <- dplyr::bind_rows(design_checks, output_checks)
atomic_write_csv(all_checks, r2r2_result_file("summaries", "delegated_C_D_smoke_checks.csv"))

failure_table <- method_metrics |>
  dplyr::filter(!.data$fit_success) |>
  dplyr::count(.data$method, .data$failure_stage, .data$error_class, .data$error_message, name = "jobs")
runtime_table <- method_metrics |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(
    jobs = dplyr::n(),
    success = sum(.data$fit_success),
    failure = sum(!.data$fit_success),
    median_runtime_sec = stats::median(.data$runtime_sec, na.rm = TRUE),
    max_runtime_sec = safe_max(.data$runtime_sec),
    .groups = "drop"
  )

report_path <- file.path(r2r2_root(), "results", "summaries", "delegated_C_D_smoke_report.md")
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "# Delegated Battery C/D Smoke Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Local SGL-WQS source: `", r2r2_sglwqs_source_path(), "`"),
  paste0("SGL-WQS version: `", local_source_version, "`"),
  paste0("Smoke settings: n_boot=", settings$n_boot, ", workers=", workers),
  paste0("Namespaced raw output: `", r2r2_result_file("raw"), "`"),
  "",
  "## Checks",
  "",
  capture.output(knitr::kable(all_checks, format = "markdown")),
  "",
  "## Wrapper Runtime and Completion",
  "",
  capture.output(knitr::kable(runtime_table, format = "markdown", digits = 3)),
  "",
  "## Recorded Fit Failures",
  "",
  if (nrow(failure_table)) capture.output(knitr::kable(failure_table, format = "markdown")) else "No fit failures were recorded.",
  "",
  "## Interpretation Boundary",
  "",
  "- These reduced-bootstrap jobs are implementation checks only and are not production evidence.",
  "- Validation-stage SGL-WQS p-values are conditional/exploratory outputs, not formal selective inference.",
  "- qgcomp attribution remains labelled coefficient-derived and is not treated as a constrained WQS weight."
), report_path)

if (!all(all_checks$passed)) {
  stop("One or more delegated C/D smoke checks failed; see ", report_path, call. = FALSE)
}

message("Delegated C/D smoke test passed: ", report_path)
