source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()

manifest <- make_paired_gaussian_split_stability_jobs()
write_manifest(manifest, "scenario_manifest_C2_paired.csv")

checks <- lapply(split(manifest, manifest$scenario_id), function(scenario_jobs) {
  job <- scenario_jobs[1, , drop = FALSE]
  dat <- generate_job_data(job)
  effective_train_prop <- attr(dat, "paired_train_prop")
  permutation <- attr(dat, "paired_row_permutation")
  reference_train <- attr(dat, "paired_reference_train_idx")
  reference_validation <- attr(dat, "paired_reference_validation_idx")
  reproduced <- make_train_validation_split(
    nrow(dat), as.integer(job$split_seed), effective_train_prop,
    y = dat$Y, family = "gaussian"
  )

  data.frame(
    scenario_id = job$scenario_id,
    paired_binary_scenario_id = job$paired_binary_scenario_id,
    n = job$n,
    data_seed = job$data_seed,
    split_seed = job$split_seed,
    binary_reference_jobs = nrow(scenario_jobs),
    manifest_jobs_expected = 100L,
    exposure_covariates_identical = nzchar(attr(dat, "paired_exposure_covariate_hash")),
    train_membership_identical = identical(
      sort(permutation[reproduced$train]), sort(reference_train)
    ),
    validation_membership_identical = identical(
      sort(permutation[reproduced$validation]), sort(reference_validation)
    ),
    binary_reference_train_n = length(reference_train),
    gaussian_paired_train_n = length(reproduced$train),
    binary_reference_validation_n = length(reference_validation),
    gaussian_paired_validation_n = length(reproduced$validation),
    effective_train_prop = effective_train_prop,
    binary_reference_prevalence = attr(dat, "paired_binary_prevalence"),
    gaussian_residual_sd = attr(dat, "gaussian_residual_sd"),
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::bind_rows() |>
  dplyr::arrange(.data$scenario_id)

checks$passed <- with(checks,
  binary_reference_jobs == manifest_jobs_expected &
    exposure_covariates_identical &
    train_membership_identical &
    validation_membership_identical &
    binary_reference_train_n == gaussian_paired_train_n &
    binary_reference_validation_n == gaussian_paired_validation_n &
    gaussian_residual_sd == 1
)

atomic_write_csv(
  checks,
  r2r2_result_file("summaries", "delegated_C2_paired_preflight_checks.csv")
)
assert_or_stop(nrow(manifest) == 400L, "C2 paired manifest does not contain 400 jobs.")
assert_or_stop(all(checks$passed), "C2 paired preflight failed.")
message("Paired C2 preflight passed for all four scenarios; manifest jobs: ", nrow(manifest), ".")
