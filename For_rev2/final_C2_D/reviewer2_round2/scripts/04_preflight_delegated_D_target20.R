source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
r2r2_load_local_sglwqs()

if (nzchar(Sys.getenv("R2R2_RESULT_NAMESPACE", unset = ""))) {
  stop("Target-20 production preflight requires R2R2_RESULT_NAMESPACE to be unset.", call. = FALSE)
}

manifest <- make_production_manifest() |>
  dplyr::filter(.data$battery == "D_unbalanced_effect_factorial", .data$run_this)
binomial_manifest <- manifest |>
  dplyr::filter(.data$family == "binomial")
dataset_jobs <- binomial_manifest |>
  dplyr::distinct(.data$scenario_id, .data$data_seed, .keep_all = TRUE)

dataset_diagnostics <- purrr::pmap_dfr(dataset_jobs, function(...) {
  job <- data.frame(..., stringsAsFactors = FALSE)
  dat <- generate_job_data(job)
  qa_generated_data(dat, job)
  data.frame(
    scenario_id = job$scenario_id,
    data_seed = as.integer(job$data_seed),
    target_prevalence = attr(dat, "binomial_target_prevalence"),
    expected_prevalence = attr(dat, "binomial_expected_prevalence"),
    realized_prevalence = mean(dat$Y),
    positive_count = sum(dat$Y),
    calibrated_intercept = attr(dat, "binomial_intercept"),
    data_hash = digest::digest(dat, algo = "sha256"),
    stringsAsFactors = FALSE
  )
})

paired_jobs <- binomial_manifest |>
  dplyr::filter(
    .data$scenario_id == dplyr::first(.data$scenario_id),
    .data$data_seed == dplyr::first(.data$data_seed)
  )
paired_hashes <- vapply(
  split(paired_jobs, seq_len(nrow(paired_jobs))),
  function(job) digest::digest(generate_job_data(job), algo = "sha256"),
  character(1)
)

checks <- data.frame(
  check = c(
    "manifest_expected_jobs",
    "manifest_binomial_target_20pct",
    "expected_binomial_datasets",
    "expected_prevalence_calibrated",
    "realized_prevalence_in_preflight_range",
    "identical_data_across_methods",
    "local_source_is_requested_version"
  ),
  passed = c(
    nrow(manifest) == 1440L,
    all(abs(binomial_manifest$binomial_target_prevalence - 0.20) < 1e-12),
    nrow(dataset_diagnostics) == 180L,
    all(abs(dataset_diagnostics$expected_prevalence - 0.20) < 1e-10),
    all(dataset_diagnostics$realized_prevalence >= 0.15) &&
      all(dataset_diagnostics$realized_prevalence <= 0.25),
    length(unique(paired_hashes)) == 1L,
    identical(as.character(utils::packageVersion("sglwqs")), "0.8.13.9001")
  ),
  detail = c(
    sprintf("%d/1440", nrow(manifest)),
    sprintf("range %.4f--%.4f", min(binomial_manifest$binomial_target_prevalence), max(binomial_manifest$binomial_target_prevalence)),
    sprintf("%d/180", nrow(dataset_diagnostics)),
    sprintf("range %.12f--%.12f", min(dataset_diagnostics$expected_prevalence), max(dataset_diagnostics$expected_prevalence)),
    sprintf("range %.4f--%.4f", min(dataset_diagnostics$realized_prevalence), max(dataset_diagnostics$realized_prevalence)),
    paste(unique(paired_hashes), collapse = ","),
    as.character(utils::packageVersion("sglwqs"))
  ),
  stringsAsFactors = FALSE
)

atomic_write_csv(
  dataset_diagnostics,
  r2r2_result_file("summaries", "delegated_D_target20_preflight_datasets.csv")
)
atomic_write_csv(
  checks,
  r2r2_result_file("summaries", "delegated_D_target20_preflight_checks.csv")
)

if (!all(checks$passed)) {
  stop("Delegated Battery D target-20 preflight failed.", call. = FALSE)
}

message(
  "Delegated Battery D target-20 preflight passed: realized prevalence range ",
  sprintf("%.3f--%.3f", min(dataset_diagnostics$realized_prevalence), max(dataset_diagnostics$realized_prevalence)),
  "."
)
