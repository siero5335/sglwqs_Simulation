if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("gaussian_primary_benchmark_generators.R")
source_r2r2("gaussian_primary_benchmark_helpers.R")
install_gaussian_primary_benchmark_dispatch(.GlobalEnv)

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("gaussian_primary_benchmarks")

workers <- as.integer(Sys.getenv("R2R2_WORKERS", "15"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings()
manifest <- make_gaussian_primary_benchmark_manifest()
write_manifest(manifest, "gaussian_primary_benchmark_manifest.csv")

invisible(run_gaussian_primary_benchmark_jobs(
  manifest,
  settings = settings,
  workers = workers,
  force = force
))
status <- summarize_gaussian_primary_benchmark_status(manifest)
cat(sprintf(
  "Gaussian primary benchmarks complete: %d/%d valid atomic jobs.\n",
  sum(status$valid_complete),
  nrow(status)
))
