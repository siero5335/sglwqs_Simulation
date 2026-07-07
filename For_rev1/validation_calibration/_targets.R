root <- normalizePath(Sys.getenv("VALIDATION_CALIBRATION_ROOT", getwd()), mustWork = TRUE)
Sys.setenv(VALIDATION_CALIBRATION_ROOT = root)
source(file.path(root, "R", "00_paths.R"))
activate_local_lib()
source_vc("01_profiles.R")
source_vc("02_dgp.R")
source_vc("03_fit_sglwqs.R")
source_vc("04_summaries.R")

targets::tar_option_set(
  packages = c(
    "MASS", "Matrix", "dplyr", "readr", "pkgload"
  ),
  format = "rds"
)

list(
  targets::tar_target(
    profile,
    get_profile(Sys.getenv("VALIDATION_CALIBRATION_PROFILE", "smoke"))
  ),
  targets::tar_target(
    replicate_jobs,
    make_replicate_jobs(profile)
  ),
  targets::tar_target(
    replicate_job,
    split(replicate_jobs, seq_len(nrow(replicate_jobs))),
    iteration = "list"
  ),
  targets::tar_target(
    replicate_result,
    run_sglwqs_job(replicate_job, profile),
    pattern = map(replicate_job),
    iteration = "list"
  ),
  targets::tar_target(
    split_jobs,
    make_split_stability_jobs(profile)
  ),
  targets::tar_target(
    split_job,
    split(split_jobs, seq_len(nrow(split_jobs))),
    iteration = "list"
  ),
  targets::tar_target(
    split_result,
    run_sglwqs_job(split_job, profile),
    pattern = map(split_job),
    iteration = "list"
  ),
  targets::tar_target(
    all_results,
    c(replicate_result, split_result),
    iteration = "list"
  ),
  targets::tar_target(
    summary_tables,
    write_summary_tables(all_results, profile)
  )
)
