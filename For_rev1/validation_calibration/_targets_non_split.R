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

write_non_split_summary_tables <- function(results, profile) {
  ensure_vc_dirs()
  directions <- bind_direction_results(results)
  diagnostics <- bind_diagnostics(results)
  calibration <- summarize_calibration(directions)
  runtime <- summarize_runtime(diagnostics, profile)

  prefix <- profile$name
  atomic_write_csv(
    directions,
    vc_file("output", "tables", paste0("validation_replicate_results_", prefix, ".csv"))
  )
  atomic_write_csv(
    diagnostics,
    vc_file("output", "tables", paste0("validation_diagnostics_", prefix, ".csv"))
  )
  atomic_write_csv(
    calibration,
    vc_file("output", "tables", paste0("validation_calibration_summary_", prefix, ".csv"))
  )
  atomic_write_csv(
    runtime,
    vc_file("output", "tables", paste0("validation_runtime_projection_", prefix, ".csv"))
  )

  list(
    directions = directions,
    diagnostics = diagnostics,
    calibration = calibration,
    runtime = runtime
  )
}

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
    summary_tables,
    write_non_split_summary_tables(replicate_result, profile)
  )
)

