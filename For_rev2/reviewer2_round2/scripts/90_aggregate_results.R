source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("metric_helpers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()
write_session_artifacts("aggregate")

manifest_path <- r2r2_file("config", "scenario_manifest.csv")
manifest <- if (file.exists(manifest_path)) {
  readr::read_csv(manifest_path, show_col_types = FALSE)
} else {
  make_production_manifest()
}
aggregate_and_write_outputs(manifest)
message("Aggregation complete: ", r2r2_file("REVIEWER2_ROUND2_SIMULATION_REPORT.md"))
