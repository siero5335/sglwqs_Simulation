if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("battery_B")

workers <- as.integer(Sys.getenv("R2R2_WORKERS", "15"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings()
manifest <- make_production_manifest() |> dplyr::filter(.data$battery %in% c("B1_gaussian_highdim", "B2_weak_gaussian_highdim"))
invisible(run_job_table(manifest, settings = settings, workers = workers, force = force))
summarize_manifest_status(make_production_manifest())
