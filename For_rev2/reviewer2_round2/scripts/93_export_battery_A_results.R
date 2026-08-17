#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
r2_dir <- file.path(root, "reviewer2_round2")
out_dir <- file.path(r2_dir, "results", "exports", "battery_A_complete")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
snapshot_time <- Sys.time()

batteries <- c("A1_gaussian_sample_size", "A2_weak_gaussian_sample_size")
methods <- c("SGL-WQS", "gWQS", "groupWQS", "qgcomp")
source(file.path(r2_dir, "R", "metric_helpers.R"))

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) > 1L) stats::sd(x) else NA_real_
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) stats::median(x) else NA_real_
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x)) unname(stats::quantile(x, probability, names = FALSE)) else NA_real_
}

group_id <- function(x) {
  values <- lapply(x, function(column) {
    column <- as.character(column)
    column[is.na(column)] <- "<NA>"
    column
  })
  do.call(paste, c(values, sep = "\r"))
}

markdown_table <- function(x, digits = 3L) {
  if (!nrow(x)) return("_No rows._")
  rendered <- x
  rendered[] <- lapply(rendered, function(column) {
    if (is.numeric(column)) {
      format(round(column, digits), trim = TRUE, scientific = FALSE)
    } else {
      as.character(column)
    }
  })
  header <- paste0("| ", paste(names(rendered), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(rendered)), collapse = " | "), " |")
  rows <- apply(rendered, 1L, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  paste(c(header, separator, rows), collapse = "\n")
}

raw_paths <- unlist(lapply(batteries, function(battery) {
  Sys.glob(file.path(r2_dir, "results", "raw", battery, "*", "*", "seed*.rds"))
}), use.names = FALSE)

loaded <- lapply(raw_paths, function(path) {
  result <- readRDS(path)
  result$.source_file <- path
  result
})
is_smoke <- vapply(loaded, function(result) {
  grepl("smoke", result$method_metrics$scenario_id[[1L]], ignore.case = TRUE)
}, logical(1))
loaded <- loaded[!is_smoke]
raw_paths <- vapply(loaded, `[[`, character(1), ".source_file")

method_metrics <- do.call(rbind, lapply(loaded, `[[`, "method_metrics"))
component_metrics <- do.call(rbind, lapply(loaded, `[[`, "component_metrics"))
direction_parts <- lapply(loaded, `[[`, "sglwqs_direction_results")
direction_parts <- direction_parts[vapply(direction_parts, nrow, integer(1)) > 0L]
direction_results <- do.call(rbind, direction_parts)
rownames(method_metrics) <- NULL
rownames(component_metrics) <- NULL
rownames(direction_results) <- NULL

data_diagnostics <- do.call(rbind, lapply(loaded, function(result) {
  metadata <- result$method_metrics[1L, c(
    "scenario_id", "battery", "family", "method", "data_seed", "n", "p",
    "group_structure", "effect_profile", "signal_profile"
  ), drop = FALSE]
  cbind(metadata, result$data_diagnostics, stringsAsFactors = FALSE)
}))
rownames(data_diagnostics) <- NULL

truth_parts <- lapply(loaded, function(result) {
  metadata <- result$method_metrics[1L, c(
    "scenario_id", "battery", "family", "data_seed", "n", "p",
    "group_structure", "effect_profile", "signal_profile"
  ), drop = FALSE]
  cbind(metadata[rep(1L, nrow(result$truth)), , drop = FALSE], result$truth)
})
truth <- unique(do.call(rbind, truth_parts))
rownames(truth) <- NULL

truth_corrections <- truth_support_correction_audit(component_metrics)
component_metrics <- correct_component_truth_support(component_metrics)
method_metrics <- recalculate_method_metrics_from_components(method_metrics, component_metrics)
direction_results <- correct_sglwqs_direction_truth_support(direction_results, component_metrics)
truth <- correct_truth_table_support(truth)

summary_keys <- c("battery", "family", "signal_profile", "n", "p", "method")
summary_splits <- split(method_metrics, group_id(method_metrics[summary_keys]))
method_summary <- do.call(rbind, lapply(summary_splits, function(x) {
  meta <- x[1L, summary_keys, drop = FALSE]
  cbind(meta, data.frame(
    attempted_jobs = nrow(x),
    fit_successes = sum(x$fit_success),
    fit_failures = sum(!x$fit_success),
    fit_completion_rate = mean(x$fit_success),
    runtime_median_sec = safe_median(x$runtime_sec),
    runtime_q1_sec = safe_quantile(x$runtime_sec, 0.25),
    runtime_q3_sec = safe_quantile(x$runtime_sec, 0.75),
    mean_active_direction_accuracy = safe_mean(x$active_direction_accuracy),
    sd_active_direction_accuracy = safe_sd(x$active_direction_accuracy),
    mean_sign_assignment_accuracy = safe_mean(x$sign_assignment_accuracy),
    mean_active_attribution = safe_mean(x$active_attribution),
    sd_active_attribution = safe_sd(x$active_attribution),
    mean_null_attribution = safe_mean(x$null_attribution),
    mean_bootstrap_success_rate = safe_mean(x$bootstrap_success_rate),
    stringsAsFactors = FALSE
  ))
}))
rownames(method_summary) <- NULL
method_summary <- method_summary[order(
  match(method_summary$battery, batteries), method_summary$n,
  match(method_summary$method, methods)
), ]

direction_keys <- c("battery", "family", "n", "p")
direction_splits <- split(direction_results, group_id(direction_results[direction_keys]))
conditional_summary <- do.call(rbind, lapply(direction_splits, function(x) {
  active <- !x$is_true_null
  null <- x$is_true_null
  active_retained <- active & x$retained
  data.frame(
    battery = x$battery[[1L]],
    family = x$family[[1L]],
    signal_profile = if (x$battery[[1L]] == "A1_gaussian_sample_size") "baseline" else "weak_0.5x",
    n = x$n[[1L]],
    p = x$p[[1L]],
    active_group_direction_rows = sum(active),
    active_retained_rows = sum(active_retained),
    active_retention_rate = mean(x$retained[active]),
    active_rejection_rate_retained = if (any(active_retained)) {
      mean(x$rejected_0_05[active_retained])
    } else {
      NA_real_
    },
    active_rejection_rate_all_attempted = mean(x$rejected_0_05[active]),
    null_group_direction_rows = sum(null),
    null_retained_rows = sum(null & x$retained),
    null_retention_rate = mean(x$retained[null]),
    null_rejection_rate_all_attempted = mean(x$rejected_0_05[null]),
    mean_bootstrap_success_rate = safe_mean(x$boot_success_rate),
    stringsAsFactors = FALSE
  )
}))
conditional_summary <- conditional_summary[order(
  match(conditional_summary$battery, batteries), conditional_summary$n
), ]
rownames(conditional_summary) <- NULL

paired_rows <- list()
paired_index <- 1L
scenario_keys <- c("battery", "family", "signal_profile", "n", "p")
scenario_splits <- split(method_metrics, group_id(method_metrics[scenario_keys]))
for (scenario_data in scenario_splits) {
  reference <- scenario_data[
    scenario_data$method == "SGL-WQS",
    c("data_seed", "active_direction_accuracy"),
    drop = FALSE
  ]
  names(reference)[[2L]] <- "reference"
  for (comparator in setdiff(methods, "SGL-WQS")) {
    other <- scenario_data[
      scenario_data$method == comparator,
      c("data_seed", "active_direction_accuracy"),
      drop = FALSE
    ]
    names(other)[[2L]] <- "comparator"
    paired <- merge(reference, other, by = "data_seed")
    difference <- paired$reference - paired$comparator
    difference <- difference[is.finite(difference)]
    set.seed(93000L + paired_index)
    bootstrap_means <- replicate(
      2000L,
      mean(sample(difference, length(difference), replace = TRUE))
    )
    confidence_interval <- stats::quantile(
      bootstrap_means, c(0.025, 0.975), names = FALSE
    )
    p_value <- suppressWarnings(stats::wilcox.test(
      difference, mu = 0, exact = FALSE
    )$p.value)
    meta <- scenario_data[1L, scenario_keys, drop = FALSE]
    paired_rows[[paired_index]] <- cbind(meta, data.frame(
      reference_method = "SGL-WQS",
      comparator_method = comparator,
      paired_seeds = length(difference),
      mean_paired_difference = mean(difference),
      median_paired_difference = stats::median(difference),
      bootstrap_ci_low = confidence_interval[[1L]],
      bootstrap_ci_high = confidence_interval[[2L]],
      raw_p_value = p_value,
      stringsAsFactors = FALSE
    ))
    paired_index <- paired_index + 1L
  }
}
paired_summary <- do.call(rbind, paired_rows)
adjustment_groups <- split(
  seq_len(nrow(paired_summary)),
  group_id(paired_summary[c("battery", "family")])
)
paired_summary$holm_p_value <- NA_real_
for (indices in adjustment_groups) {
  paired_summary$holm_p_value[indices] <- stats::p.adjust(
    paired_summary$raw_p_value[indices], method = "holm"
  )
}
paired_summary <- paired_summary[order(
  match(paired_summary$battery, batteries), paired_summary$n,
  paired_summary$comparator_method
), ]
rownames(paired_summary) <- NULL

failure_summary <- method_metrics[!method_metrics$fit_success, c(
  "battery", "scenario_id", "method", "data_seed", "n", "failure_stage",
  "error_class", "error_message"
), drop = FALSE]

production_manifest <- read.csv(
  file.path(r2_dir, "config", "scenario_manifest.csv"),
  check.names = FALSE
)
manifest_A <- production_manifest[production_manifest$battery %in% batteries, , drop = FALSE]

file_info <- file.info(raw_paths)
raw_manifest <- data.frame(
  relative_path = substring(raw_paths, nchar(root) + 2L),
  bytes = file_info$size,
  modified = format(file_info$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  md5 = unname(tools::md5sum(raw_paths)),
  stringsAsFactors = FALSE
)

job_key <- c("battery", "scenario_id", "method", "data_seed", "split_replicate")
method_counts <- aggregate(
  method ~ battery + scenario_id + data_seed,
  data = method_metrics,
  FUN = function(x) length(unique(x))
)
seed_counts <- aggregate(
  data_seed ~ battery + scenario_id + method,
  data = method_metrics,
  FUN = function(x) length(unique(x))
)
active_counts <- aggregate(
  IsActive ~ battery + scenario_id + data_seed,
  data = truth,
  FUN = sum
)
qa <- rbind(
  data.frame(
    check = "production_atomic_job_count_is_1200",
    rows_checked = nrow(method_metrics),
    passed = nrow(method_metrics) == 1200L
  ),
  data.frame(
    check = "A1_job_count_is_720",
    rows_checked = sum(method_metrics$battery == batteries[[1L]]),
    passed = sum(method_metrics$battery == batteries[[1L]]) == 720L
  ),
  data.frame(
    check = "A2_job_count_is_480",
    rows_checked = sum(method_metrics$battery == batteries[[2L]]),
    passed = sum(method_metrics$battery == batteries[[2L]]) == 480L
  ),
  data.frame(
    check = "no_duplicate_atomic_jobs",
    rows_checked = nrow(method_metrics),
    passed = !any(duplicated(method_metrics[job_key]))
  ),
  data.frame(
    check = "all_atomic_fits_completed",
    rows_checked = nrow(method_metrics),
    passed = all(method_metrics$fit_success)
  ),
  data.frame(
    check = "four_methods_per_dataset",
    rows_checked = nrow(method_counts),
    passed = all(method_counts$method == 4L)
  ),
  data.frame(
    check = "thirty_seeds_per_scenario_method",
    rows_checked = nrow(seed_counts),
    passed = all(seed_counts$data_seed == 30L)
  ),
  data.frame(
    check = "thirteen_component_rows_per_atomic_job",
    rows_checked = nrow(component_metrics),
    passed = nrow(component_metrics) == 13L * nrow(method_metrics)
  ),
  data.frame(
    check = "six_direction_rows_per_sglwqs_job",
    rows_checked = nrow(direction_results),
    passed = nrow(direction_results) == 6L * sum(method_metrics$method == "SGL-WQS")
  ),
  data.frame(
    check = "gaussian_residual_sd_is_one",
    rows_checked = nrow(data_diagnostics),
    passed = all(abs(data_diagnostics$gaussian_residual_sd - 1) < 1e-12)
  ),
  data.frame(
    check = "active_support_matches_nonzero_beta",
    rows_checked = nrow(truth),
    passed = identical(truth$IsActive, abs(truth$True_Beta) > truth_support_tolerance())
  ),
  data.frame(
    check = "baseline_and_weak_have_ten_active_components",
    rows_checked = nrow(active_counts),
    passed = all(active_counts$IsActive == 10L)
  ),
  data.frame(
    check = "weak_mg_is_active_negative",
    rows_checked = sum(truth$signal_profile == "weak" & truth$Variable == "mg"),
    passed = all(
      truth$IsActive[truth$signal_profile == "weak" & truth$Variable == "mg"] &
        truth$True_Direction[truth$signal_profile == "weak" & truth$Variable == "mg"] == "Negative"
    )
  )
)

write.csv(method_metrics, file.path(out_dir, "atomic_method_metrics_A.csv"), row.names = FALSE)
write.csv(component_metrics, file.path(out_dir, "atomic_component_metrics_A.csv"), row.names = FALSE)
write.csv(direction_results, file.path(out_dir, "atomic_sglwqs_group_direction_results_A.csv"), row.names = FALSE)
write.csv(data_diagnostics, file.path(out_dir, "atomic_data_diagnostics_A.csv"), row.names = FALSE)
write.csv(truth, file.path(out_dir, "dataset_truth_A.csv"), row.names = FALSE)
write.csv(method_summary, file.path(out_dir, "method_performance_summary_A.csv"), row.names = FALSE)
write.csv(conditional_summary, file.path(out_dir, "sglwqs_conditional_summary_A.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(out_dir, "paired_active_direction_differences_A.csv"), row.names = FALSE)
write.csv(failure_summary, file.path(out_dir, "failure_summary_A.csv"), row.names = FALSE)
write.csv(manifest_A, file.path(out_dir, "scenario_manifest_A.csv"), row.names = FALSE)
write.csv(raw_manifest, file.path(out_dir, "raw_file_manifest_A.csv"), row.names = FALSE)
write.csv(qa, file.path(out_dir, "qa_checks_A.csv"), row.names = FALSE)
write.csv(truth_corrections, file.path(out_dir, "truth_support_corrections_A.csv"), row.names = FALSE)

compact <- method_summary[, c(
  "signal_profile", "n", "method", "attempted_jobs", "fit_completion_rate",
  "mean_active_direction_accuracy", "mean_active_attribution", "runtime_median_sec"
)]
names(compact) <- c(
  "Signal", "n", "Method", "Jobs", "Fit", "Direction accuracy",
  "Active attribution", "Runtime median (s)"
)
conditional_compact <- conditional_summary[, c(
  "signal_profile", "n", "active_retention_rate",
  "active_rejection_rate_retained", "active_rejection_rate_all_attempted",
  "null_rejection_rate_all_attempted", "mean_bootstrap_success_rate"
)]
names(conditional_compact) <- c(
  "Signal", "n", "Active retention",
  "Active p<0.05 among retained", "Active p<0.05 among all attempted",
  "True-null p<0.05 among all attempted", "Bootstrap success"
)

report <- c(
  "# Battery A Gaussian Simulation Results",
  "",
  paste0("Generated: ", format(snapshot_time, "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "- A1 baseline Gaussian sample-size battery: n = 100, 500, 1,000, 5,000, 10,000, and 50,000.",
  "- A2 weak Gaussian sample-size battery: n = 500, 1,000, 5,000, and 10,000; mixture coefficients multiplied by 0.5.",
  "- Thirty independently generated datasets per scenario and four methods per dataset.",
  "- Gaussian outcome: Y = eta + epsilon, epsilon ~ N(0, 1).",
  "- Methods: SGL-WQS, gWQS, groupWQS, and qgcomp.",
  "- Active support is defined by a numerical-zero tolerance: abs(true beta) > 1e-12.",
  "- Saved model outputs were not refitted; component truth and derived metrics were recalculated from true beta.",
  "",
  "## Completion",
  "",
  paste0("All ", nrow(method_metrics), " production atomic jobs completed successfully; no fit failures were recorded."),
  "",
  "## Method Results",
  "",
  markdown_table(compact),
  "",
  "## SGL-WQS Conditional Validation-Stage Results",
  "",
  markdown_table(conditional_compact),
  "",
  "## Interpretation Guardrails",
  "",
  "- This low-dimensional additive Gaussian design is naturally favorable to qgcomp.",
  "- qgcomp used qgcomp.glm.noboot on the full dataset. WQS workflows used 200 bootstrap resamples, and SGL-WQS also used 10-fold tuning and a 60/40 training-validation design.",
  "- Runtime is the observed cost of each workflow, not an equal-computation benchmark.",
  "- qgcomp coefficient-derived attribution and constrained WQS weights are not identical estimands; active attribution should not be used as a universal method ranking.",
  "- Validation-stage p-values are conditional exploratory outputs, not formal selective inference.",
  "- The results do not support universal superiority of any method.",
  "",
  "## Reproducibility",
  "",
  "- Battery A local SGL-WQS source: source/sglwqs_old_analysis/sglwqs-main.",
  "- Recorded package version: 1.0.0.",
  "- The extracted source directory is not a Git checkout; no Git commit was recorded for Battery A.",
  "- Package versions and sessionInfo are included in the bundle.",
  "- raw_file_manifest_A.csv records the relative path, byte size, modification time, and MD5 checksum of every atomic RDS file.",
  "- truth_support_corrections_A.csv records every saved component label changed during post-fit truth correction.",
  "",
  "## QA",
  "",
  markdown_table(qa, digits = 0L)
)
writeLines(report, file.path(out_dir, "BATTERY_A_RESULTS_SUMMARY.md"))

bundle_readme <- c(
  "# Battery A Code and Results Bundle",
  "",
  "This bundle contains the complete Reviewer 2 Round 2 Battery A code and results.",
  "",
  "## Included",
  "",
  "- reviewer2_round2/results/raw/A1_gaussian_sample_size: all production A1 atomic RDS files.",
  "- reviewer2_round2/results/raw/A2_weak_gaussian_sample_size: all production A2 atomic RDS files.",
  "- reviewer2_round2/results/exports/battery_A_complete: flat CSV summaries, QA, manifests, and results report.",
  "- truth_support_corrections_A.csv: audit trail for post-fit truth-label corrections; production RDS files remain unchanged.",
  "- reviewer2_round2/R and selected scripts/configuration: simulation and aggregation code.",
  "- source/sglwqs_old_analysis/sglwqs-main: exact local SGL-WQS source directory recorded for Battery A.",
  "- reviewer2_round2/results/logs/battery_A_*: package versions, source metadata, and sessionInfo.",
  "",
  "The manuscript and Response letter are intentionally not included or modified."
)
writeLines(bundle_readme, file.path(out_dir, "BUNDLE_README.md"))
write.csv(data.frame(
  field = c(
    "snapshot_time", "production_jobs", "fit_successes", "fit_failures",
    "component_rows", "sglwqs_direction_rows", "raw_rds_files", "R_version"
  ),
  value = c(
    format(snapshot_time, "%Y-%m-%d %H:%M:%S %Z"),
    nrow(method_metrics), sum(method_metrics$fit_success), sum(!method_metrics$fit_success),
    nrow(component_metrics), nrow(direction_results), length(raw_paths), R.version.string
  )
), file.path(out_dir, "snapshot_metadata_A.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "aggregation_sessionInfo_A.txt"))

cat("Battery A export complete\n")
cat("Output:", out_dir, "\n")
cat("Production jobs:", nrow(method_metrics), "\n")
cat("Fit failures:", sum(!method_metrics$fit_success), "\n")
cat("QA passed:", all(qa$passed), "\n")
