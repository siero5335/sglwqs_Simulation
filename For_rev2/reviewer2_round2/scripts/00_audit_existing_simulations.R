if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("method_wrappers.R")
source_r2r2("metric_helpers.R")

r2r2_set_thread_env()
r2r2_load_packages()
r2r2_make_dirs()

manifest <- make_production_manifest()
write_manifest(manifest, "scenario_manifest.csv")
write_session_artifacts("audit")

existing_files <- data.frame(
  role = c(
    "binary_sample_size_rmd",
    "highdim_p_checkpointed_rmd",
    "validation_dgp",
    "validation_split_stability_rmd",
    "old_sglwqs_source",
    "envint_revision_sglwqs_source",
    "private_sglwqs_source"
  ),
  path = c(
    "results/final_res/compare_mixture_methods_by_samplesize_30seeds.Rmd",
    "analysis/compare_mixture_methods_highdim_p_checkpointed.Rmd",
    "analysis/validation_calibration/R/02_dgp.R",
    "analysis/validation_split_stability_standalone.Rmd",
    "source/sglwqs_old_analysis/sglwqs-main",
    "source/sglwqs-envint-revision-docs",
    "source/sglwqs-private-main"
  ),
  stringsAsFactors = FALSE
)
existing_files$exists <- file.exists(file.path(r2r2_repo_root(), existing_files$path)) |
  dir.exists(file.path(r2r2_repo_root(), existing_files$path))
atomic_write_csv(existing_files, r2r2_result_file("summaries", "existing_simulation_audit_files.csv"))

manifest_summary <- manifest |>
  dplyr::count(.data$battery, .data$run_this, name = "atomic_jobs") |>
  dplyr::arrange(.data$battery, .data$run_this)
atomic_write_csv(manifest_summary, r2r2_result_file("summaries", "scenario_manifest_summary.csv"))

audit_md <- c(
  "# Existing Simulation Audit",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Existing Code Identified",
  "",
  capture.output(knitr::kable(existing_files, format = "markdown")),
  "",
  "## New Production Manifest",
  "",
  capture.output(knitr::kable(manifest_summary, format = "markdown")),
  "",
  "## Key Audit Notes",
  "",
  "- The binary sample-size Rmd defines the original 13-exposure DGP used for matched Gaussian Battery A.",
  "- The validation calibration DGP already contains a Gaussian identity branch with residual SD 1; this definition is reused for split-stability where applicable.",
  "- The repository high-dimensional checkpointed Rmd uses a Gaussian outcome; the new family-explicit high-dimensional runner is stored under `reviewer2_round2/` and does not overwrite previous outputs.",
  "- Battery A retains the original local SGL-WQS analysis source. Audit/smoke and Batteries B-E use the pinned `paper/envint-revision-docs` source at Git commit `2fdd519e520a7dad1162810643e175cd616b1154`.",
  "- Battery E n=1,000 rows are marked as `run_this = FALSE` and are intended to be reused from Battery D's unbalanced weak-heterogeneous n=1,000 scenario.",
  "- qgcomp attribution is stored as coefficient-derived attribution. WQS-family weights are stored as WQS-type constrained weights."
)
writeLines(audit_md, r2r2_result_file("summaries", "existing_simulation_audit.md"))

message("Wrote manifest: ", r2r2_file("config", "scenario_manifest.csv"))
message("Audit complete.")
