if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("dense_group_signal_generator.R")
source_r2r2("dense_group_signal_helpers.R")
install_dense_group_signal_dispatch(.GlobalEnv)

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("dense_group_signal")

workers <- as.integer(Sys.getenv("R2R2_WORKERS", "15"))
force <- Sys.getenv("R2R2_FORCE", "false") %in% c("1", "true", "TRUE")
settings <- r2r2_settings()
selected_families <- dense_group_signal_selected_families()
full_manifest <- make_dense_group_signal_manifest()
write_manifest(full_manifest, "dense_group_signal_manifest.csv")
manifest <- full_manifest[full_manifest$family %in% selected_families, , drop = FALSE]

invisible(run_dense_group_signal_jobs(
  manifest,
  settings = settings,
  workers = workers,
  force = force
))
family_label <- paste(selected_families, collapse = "_")
status <- summarize_dense_group_signal_status(
  manifest,
  sprintf("dense_group_signal_completion_status_%s.csv", family_label)
)
cat(sprintf(
  "Dense group-signal battery (%s) complete: %d/%d valid atomic jobs.\n",
  paste(selected_families, collapse = ", "),
  sum(status$valid_complete),
  nrow(status)
))
