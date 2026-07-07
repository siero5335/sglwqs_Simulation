#!/usr/bin/env Rscript

cache_dir <- file.path("analysis", "analysis_cache", "highdim_correlation_checkpointed")
if (!dir.exists(cache_dir)) stop("Correlation checkpoint cache directory was not found.")

status_files <- list.files(
  cache_dir,
  pattern = "^highdim_correlation_checkpointed_status_.*\\.csv$",
  full.names = TRUE
)
if (length(status_files) == 0L) stop("No correlation checkpoint status CSV was found.")

status_file <- status_files[which.max(file.info(status_files)$mtime)]
status <- read.csv(status_file, stringsAsFactors = FALSE)

resolve_path <- function(path) {
  candidates <- c(path, file.path("analysis", path))
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) return(path)
  hit[[1]]
}
status$Resolved_File <- vapply(status$Scenario_File, resolve_path, character(1))
status$Completed <- file.exists(status$Resolved_File)

out_dir <- file.path("analysis", "results", "highdim_correlation_progress")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(status, file.path(out_dir, "checkpoint_status_latest.csv"), row.names = FALSE)

group_cols <- c("Scenario", "p", "n", "n_groups", "truth_profile")
if ("within_rho" %in% names(status)) {
  group_cols <- c(group_cols, "within_rho")
}
if ("cross_rho" %in% names(status)) {
  group_cols <- c(group_cols, "cross_rho")
}
if ("seed_count" %in% names(status)) {
  group_cols <- c(group_cols, "seed_count")
}

progress <- aggregate(
  status["Completed"],
  status[group_cols],
  sum
)
names(progress)[names(progress) == "Completed"] <- "Completed_Scenarios"
expected <- aggregate(Seed ~ Scenario, status, length)
names(expected)[names(expected) == "Seed"] <- "Expected_Scenarios"
progress <- merge(progress, expected, by = "Scenario", all.x = TRUE)
progress$Remaining_Scenarios <- progress$Expected_Scenarios - progress$Completed_Scenarios
progress <- progress[order(progress$Scenario), ]
write.csv(progress, file.path(out_dir, "checkpoint_progress_by_scenario.csv"), row.names = FALSE)
print(progress, row.names = FALSE)

completed_files <- status$Resolved_File[status$Completed]
if (length(completed_files) == 0L) {
  cat("No completed scenario RDS files are available yet.\n")
  quit(save = "no")
}

results <- lapply(completed_files, readRDS)
method_levels <- c("qgcomp", "gWQS", "groupWQS", "SGL-WQS")

timing_df <- do.call(rbind, lapply(results, function(res) {
  if (length(res$timing) == 0L) return(NULL)
  do.call(rbind, lapply(names(res$timing), function(method) {
    data.frame(
      Scenario = res$Scenario,
      p = res$p,
      n = res$n,
      n_groups = res$n_groups,
      truth_profile = res$truth_profile,
      Seed = res$seed,
      Method = method,
      Seconds = as.numeric(res$timing[[method]]),
      stringsAsFactors = FALSE
    )
  }))
}))
if (!is.null(timing_df)) {
  timing_summary <- aggregate(Seconds ~ Scenario + Method, timing_df, function(x) {
    c(Median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), N = length(x))
  })
  timing_summary <- do.call(data.frame, timing_summary)
  names(timing_summary) <- sub("^Seconds\\.", "", names(timing_summary))
  write.csv(timing_df, file.path(out_dir, "checkpoint_timing_completed.csv"), row.names = FALSE)
  write.csv(timing_summary, file.path(out_dir, "checkpoint_timing_summary_completed.csv"), row.names = FALSE)
}

accuracy_df <- do.call(rbind, lapply(results, function(res) {
  if (length(res$accuracy) == 0L) return(NULL)
  do.call(rbind, lapply(names(res$accuracy), function(method) {
    data.frame(
      Scenario = res$Scenario,
      p = res$p,
      n = res$n,
      n_groups = res$n_groups,
      truth_profile = res$truth_profile,
      Seed = res$seed,
      Method = method,
      Active_Direction_Accuracy = as.numeric(res$accuracy[[method]]),
      stringsAsFactors = FALSE
    )
  }))
}))
if (!is.null(accuracy_df)) {
  accuracy_summary <- aggregate(Active_Direction_Accuracy ~ Scenario + Method, accuracy_df, function(x) {
    c(Mean = mean(x, na.rm = TRUE), Median = median(x, na.rm = TRUE), N = sum(!is.na(x)))
  })
  accuracy_summary <- do.call(data.frame, accuracy_summary)
  names(accuracy_summary) <- sub("^Active_Direction_Accuracy\\.", "", names(accuracy_summary))
  write.csv(accuracy_df, file.path(out_dir, "checkpoint_accuracy_completed.csv"), row.names = FALSE)
  write.csv(accuracy_summary, file.path(out_dir, "checkpoint_accuracy_summary_completed.csv"), row.names = FALSE)
}

method_status <- do.call(rbind, lapply(results, function(res) {
  errors <- names(res$errors)
  do.call(rbind, lapply(method_levels, function(method) {
    data.frame(
      Scenario = res$Scenario,
      p = res$p,
      n = res$n,
      n_groups = res$n_groups,
      truth_profile = res$truth_profile,
      Seed = res$seed,
      Method = method,
      Failed = method %in% errors,
      Error = if (method %in% errors) as.character(res$errors[[method]]) else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}))
failure_summary <- aggregate(Failed ~ Scenario + Method, method_status, sum)
names(failure_summary)[names(failure_summary) == "Failed"] <- "Failures"
attempts <- aggregate(Seed ~ Scenario + Method, method_status, length)
names(attempts)[names(attempts) == "Seed"] <- "Completed_Scenarios"
failure_summary <- merge(failure_summary, attempts, by = c("Scenario", "Method"), all.x = TRUE)
failure_summary$Failure_Rate <- failure_summary$Failures / failure_summary$Completed_Scenarios

write.csv(method_status, file.path(out_dir, "checkpoint_method_status_completed.csv"), row.names = FALSE)
write.csv(failure_summary, file.path(out_dir, "checkpoint_failure_summary_completed.csv"), row.names = FALSE)

cat("Wrote correlation progress summaries to:", out_dir, "\n")
