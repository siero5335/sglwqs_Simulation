#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
profile_name <- if (length(args) >= 1L) args[[1]] else "production"

root_path <- normalizePath(
  file.path("analysis", "validation_calibration"),
  mustWork = TRUE
)

Sys.setenv(
  VALIDATION_CALIBRATION_ROOT = root_path,
  VALIDATION_CALIBRATION_PROFILE = profile_name
)

source(file.path(root_path, "R", "00_paths.R"), chdir = TRUE)
source_vc("01_profiles.R")
source_vc("02_dgp.R")

profile <- get_profile(profile_name)
split_jobs <- make_split_stability_jobs(profile)

manifest_path <- function(job) {
  vc_file(
    "output", "manifests", profile$name, job$job_type,
    paste0(job$job_id, ".rds")
  )
}

split_jobs$manifest_path <- vapply(
  split(seq_len(nrow(split_jobs)), seq_len(nrow(split_jobs))),
  function(i) manifest_path(split_jobs[i, , drop = FALSE]),
  character(1)
)
split_jobs$completed <- file.exists(split_jobs$manifest_path)

progress <- aggregate(
  completed ~ scenario_id + scenario_family + n + within_r + between_r,
  data = split_jobs,
  FUN = function(x) c(completed = sum(x), total = length(x))
)
progress <- do.call(data.frame, progress)
names(progress)[names(progress) == "completed.completed"] <- "completed"
names(progress)[names(progress) == "completed.total"] <- "total"
progress$pending <- progress$total - progress$completed
progress$completion_rate <- progress$completed / progress$total
progress <- progress[order(progress$scenario_family, progress$n), ]

out_dir <- file.path("analysis", "results", "validation_split_stability_progress")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(
  out_dir,
  paste0("validation_split_stability_progress_", profile$name, ".csv")
)
write.csv(progress, out_path, row.names = FALSE, na = "")

print(progress, row.names = FALSE)
cat("Wrote:", normalizePath(out_path, mustWork = FALSE), "\n")
