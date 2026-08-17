source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("battery_C2_paired")

workers <- as.integer(Sys.getenv("R2R2_WORKERS", "1"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings()
manifest <- make_paired_gaussian_split_stability_jobs()
write_manifest(manifest, "scenario_manifest_C2_paired.csv")

invisible(run_job_table(manifest, settings = settings, workers = workers, force = force))

paths <- manifest_job_paths(manifest)
status <- manifest |>
  dplyr::mutate(
    output_path = paths,
    completed = file.exists(paths),
    valid_complete = vapply(paths, is_complete_job_file, logical(1)),
    size_bytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_)
  )
atomic_write_csv(
  status,
  r2r2_result_file("summaries", "delegated_C2_paired_completion_status.csv")
)
assert_or_stop(all(status$valid_complete), "C2 paired production ended with incomplete jobs.")
message("Paired C2 production complete: ", sum(status$valid_complete), "/", nrow(status), ".")
