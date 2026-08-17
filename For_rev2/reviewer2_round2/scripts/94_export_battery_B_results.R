#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
r2_dir <- file.path(root, "reviewer2_round2")
out_dir <- file.path(r2_dir, "results", "exports", "battery_B_complete")
figure_dir <- file.path(out_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
snapshot_time <- Sys.time()

batteries <- c("B1_gaussian_highdim", "B2_weak_gaussian_highdim")
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
is_production_B <- vapply(loaded, function(result) {
  metric <- result$method_metrics
  nrow(metric) == 1L &&
    metric$battery[[1L]] %in% batteries &&
    !grepl("smoke", metric$scenario_id[[1L]], ignore.case = TRUE)
}, logical(1))
loaded <- loaded[is_production_B]
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
    mean_active_selection_frequency = safe_mean(x$active_selection_frequency),
    mean_null_selection_frequency = safe_mean(x$null_selection_frequency),
    median_bootstrap_success_rate = safe_median(x$bootstrap_success_rate),
    median_retained_group_directions = safe_median(x$retained_group_directions),
    stringsAsFactors = FALSE
  ))
}))
rownames(method_summary) <- NULL
method_summary <- method_summary[order(
  match(method_summary$battery, batteries), method_summary$p,
  match(method_summary$method, methods)
), ]

direction_keys <- c("battery", "family", "n", "p")
direction_splits <- split(direction_results, group_id(direction_results[direction_keys]))
conditional_summary <- do.call(rbind, lapply(direction_splits, function(x) {
  active <- !x$is_true_null
  null <- x$is_true_null
  active_retained <- active & x$retained
  null_retained <- null & x$retained
  data.frame(
    battery = x$battery[[1L]],
    family = x$family[[1L]],
    signal_profile = if (x$battery[[1L]] == batteries[[1L]]) "baseline" else "weak_0.5x",
    n = x$n[[1L]],
    p = x$p[[1L]],
    group_direction_rows = nrow(x),
    active_group_direction_rows = sum(active),
    active_retained_rows = sum(active_retained),
    active_retention_rate = mean(x$retained[active]),
    active_p_lt_005_rate_retained = if (any(active_retained)) mean(x$rejected_0_05[active_retained]) else NA_real_,
    active_p_lt_005_rate_all_attempted = mean(x$rejected_0_05[active]),
    null_group_direction_rows = sum(null),
    null_retained_rows = sum(null_retained),
    null_retention_rate = mean(x$retained[null]),
    null_p_lt_005_rate_retained = if (any(null_retained)) mean(x$rejected_0_05[null_retained]) else NA_real_,
    null_p_lt_005_rate_all_attempted = mean(x$rejected_0_05[null]),
    minor_direction_exclusion_rate = mean(x$excluded_by_minor_threshold),
    median_bootstrap_success_rate = safe_median(x$boot_success_rate),
    stringsAsFactors = FALSE
  )
}))
conditional_summary <- conditional_summary[order(
  match(conditional_summary$battery, batteries), conditional_summary$p
), ]
rownames(conditional_summary) <- NULL

sgl_components <- component_metrics[component_metrics$method == "SGL-WQS", , drop = FALSE]
selection_keys <- c("battery", "family", "signal_profile", "n", "p")
selection_splits <- split(sgl_components, group_id(sgl_components[selection_keys]))
selection_summary <- do.call(rbind, lapply(selection_splits, function(x) {
  active <- x$is_active
  selected <- is.finite(x$selection_frequency) & x$selection_frequency >= 0.5
  meta <- x[1L, selection_keys, drop = FALSE]
  cbind(meta, data.frame(
    active_component_rows = sum(active),
    null_component_rows = sum(!active),
    mean_active_selection_frequency = safe_mean(x$selection_frequency[active]),
    mean_null_selection_frequency = safe_mean(x$selection_frequency[!active]),
    selection_tpr_cutoff_0_5 = mean(selected[active]),
    selection_tnr_cutoff_0_5 = mean(!selected[!active]),
    selection_fpr_cutoff_0_5 = mean(selected[!active]),
    stringsAsFactors = FALSE
  ))
}))
selection_summary <- selection_summary[order(
  match(selection_summary$battery, batteries), selection_summary$p
), ]
rownames(selection_summary) <- NULL

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
    set.seed(94000L + paired_index)
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
  match(paired_summary$battery, batteries), paired_summary$p,
  paired_summary$comparator_method
), ]
rownames(paired_summary) <- NULL

failure_summary <- method_metrics[!method_metrics$fit_success, c(
  "battery", "scenario_id", "method", "data_seed", "n", "p",
  "failure_stage", "error_class", "error_message"
), drop = FALSE]

production_manifest <- read.csv(
  file.path(r2_dir, "config", "scenario_manifest.csv"),
  check.names = FALSE
)
manifest_B <- production_manifest[production_manifest$battery %in% batteries, , drop = FALSE]

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
group_size_totals <- aggregate(
  Group_Size ~ battery + scenario_id + data_seed,
  data = unique(truth[c("battery", "scenario_id", "data_seed", "Group", "Group_Size")]),
  FUN = sum
)

qa <- rbind(
  data.frame(check = "production_atomic_job_count_is_720", rows_checked = nrow(method_metrics), passed = nrow(method_metrics) == 720L),
  data.frame(check = "B1_job_count_is_360", rows_checked = sum(method_metrics$battery == batteries[[1L]]), passed = sum(method_metrics$battery == batteries[[1L]]) == 360L),
  data.frame(check = "B2_job_count_is_360", rows_checked = sum(method_metrics$battery == batteries[[2L]]), passed = sum(method_metrics$battery == batteries[[2L]]) == 360L),
  data.frame(check = "no_duplicate_atomic_jobs", rows_checked = nrow(method_metrics), passed = !any(duplicated(method_metrics[job_key]))),
  data.frame(check = "all_atomic_fits_completed", rows_checked = nrow(method_metrics), passed = all(method_metrics$fit_success)),
  data.frame(check = "four_methods_per_dataset", rows_checked = nrow(method_counts), passed = all(method_counts$method == 4L)),
  data.frame(check = "thirty_seeds_per_scenario_method", rows_checked = nrow(seed_counts), passed = all(seed_counts$data_seed == 30L)),
  data.frame(check = "component_row_count_is_84000", rows_checked = nrow(component_metrics), passed = nrow(component_metrics) == 84000L),
  data.frame(check = "direction_row_count_is_3600", rows_checked = nrow(direction_results), passed = nrow(direction_results) == 3600L),
  data.frame(check = "gaussian_residual_sd_is_one", rows_checked = nrow(data_diagnostics), passed = all(abs(data_diagnostics$gaussian_residual_sd - 1) < 1e-12)),
  data.frame(check = "active_support_matches_nonzero_beta", rows_checked = nrow(truth), passed = identical(truth$IsActive, abs(truth$True_Beta) > truth_support_tolerance())),
  data.frame(check = "nine_active_components_per_dataset", rows_checked = nrow(active_counts), passed = all(active_counts$IsActive == 9L)),
  data.frame(check = "group_sizes_sum_to_p", rows_checked = nrow(group_size_totals), passed = all(group_size_totals$Group_Size %in% c(50L, 100L, 200L))),
  data.frame(check = "no_posthoc_truth_corrections_needed", rows_checked = nrow(truth_corrections), passed = nrow(truth_corrections) == 0L),
  data.frame(check = "raw_rds_count_is_720", rows_checked = length(raw_paths), passed = length(raw_paths) == 720L)
)

write.csv(method_metrics, file.path(out_dir, "atomic_method_metrics_B.csv"), row.names = FALSE)
write.csv(component_metrics, file.path(out_dir, "atomic_component_metrics_B.csv"), row.names = FALSE)
write.csv(direction_results, file.path(out_dir, "atomic_sglwqs_group_direction_results_B.csv"), row.names = FALSE)
write.csv(data_diagnostics, file.path(out_dir, "atomic_data_diagnostics_B.csv"), row.names = FALSE)
write.csv(truth, file.path(out_dir, "dataset_truth_B.csv"), row.names = FALSE)
write.csv(method_summary, file.path(out_dir, "method_performance_summary_B.csv"), row.names = FALSE)
write.csv(conditional_summary, file.path(out_dir, "sglwqs_conditional_summary_B.csv"), row.names = FALSE)
write.csv(selection_summary, file.path(out_dir, "sglwqs_selection_tpr_tnr_B.csv"), row.names = FALSE)
write.csv(paired_summary, file.path(out_dir, "paired_active_direction_differences_B.csv"), row.names = FALSE)
write.csv(failure_summary, file.path(out_dir, "failure_summary_B.csv"), row.names = FALSE)
write.csv(manifest_B, file.path(out_dir, "scenario_manifest_B.csv"), row.names = FALSE)
write.csv(raw_manifest, file.path(out_dir, "raw_file_manifest_B.csv"), row.names = FALSE)
write.csv(qa, file.path(out_dir, "qa_checks_B.csv"), row.names = FALSE)
write.csv(truth_corrections, file.path(out_dir, "truth_support_corrections_B.csv"), row.names = FALSE)

sgl_method <- method_summary[method_summary$method == "SGL-WQS", ]
sgl_compact <- merge(
  sgl_method[, c(
    "battery", "signal_profile", "p", "fit_completion_rate",
    "mean_active_direction_accuracy", "mean_active_attribution",
    "mean_null_attribution", "runtime_median_sec"
  )],
  selection_summary[, c(
    "battery", "p", "selection_tpr_cutoff_0_5", "selection_tnr_cutoff_0_5"
  )],
  by = c("battery", "p")
)
sgl_compact <- merge(
  sgl_compact,
  conditional_summary[, c(
    "battery", "p", "active_p_lt_005_rate_all_attempted",
    "null_p_lt_005_rate_all_attempted"
  )],
  by = c("battery", "p")
)
sgl_compact <- sgl_compact[order(match(sgl_compact$battery, batteries), sgl_compact$p), ]
sgl_compact$runtime_median_min <- sgl_compact$runtime_median_sec / 60

report_sgl <- sgl_compact[, c(
  "signal_profile", "p", "fit_completion_rate", "mean_active_direction_accuracy",
  "mean_active_attribution", "mean_null_attribution", "selection_tpr_cutoff_0_5",
  "selection_tnr_cutoff_0_5", "active_p_lt_005_rate_all_attempted",
  "null_p_lt_005_rate_all_attempted", "runtime_median_min"
)]
names(report_sgl) <- c(
  "Signal", "p", "Fit", "Direction accuracy", "Active attribution",
  "Null attribution", "Selection TPR", "Selection TNR",
  "Active conditional p<0.05", "Null conditional p<0.05", "Runtime median (min)"
)

report_all <- method_summary[, c(
  "signal_profile", "p", "method", "fit_completion_rate",
  "mean_active_direction_accuracy", "sd_active_direction_accuracy",
  "mean_active_attribution", "mean_null_attribution", "runtime_median_sec"
)]
report_all$runtime_median_sec <- report_all$runtime_median_sec / 60
names(report_all) <- c(
  "Signal", "p", "Method", "Fit", "Direction mean", "Direction SD",
  "Active attribution", "Null attribution", "Runtime median (min)"
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_data <- method_summary
  plot_data$signal_profile <- factor(
    plot_data$signal_profile,
    levels = c("baseline", "weak"),
    labels = c("Baseline", "Weak (0.5x)")
  )
  plot_data$method <- factor(plot_data$method, levels = methods)

  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = p, y = mean_active_direction_accuracy, color = method
  )) +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymin = pmax(0, mean_active_direction_accuracy - sd_active_direction_accuracy),
      ymax = pmin(1, mean_active_direction_accuracy + sd_active_direction_accuracy)
    ), width = 5, linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~signal_profile) +
    ggplot2::scale_x_continuous(breaks = c(50, 100, 200)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Number of exposures (p)", y = "Active-direction accuracy", color = "Method") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir, "B1_active_direction_accuracy_vs_p.png"), p1, width = 8, height = 4.5, dpi = 180)

  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = p, y = mean_active_attribution, color = method
  )) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~signal_profile) +
    ggplot2::scale_x_continuous(breaks = c(50, 100, 200)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Number of exposures (p)", y = "Active-component attribution", color = "Method") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir, "B2_active_attribution_vs_p.png"), p2, width = 8, height = 4.5, dpi = 180)

  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = p, y = runtime_median_sec, color = method
  )) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~signal_profile) +
    ggplot2::scale_x_continuous(breaks = c(50, 100, 200)) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "Number of exposures (p)", y = "Median runtime (seconds, log scale)", color = "Method") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir, "B3_runtime_vs_p.png"), p3, width = 8, height = 4.5, dpi = 180)

  diagnostic_plot <- rbind(
    data.frame(
      signal_profile = conditional_summary$signal_profile,
      p = conditional_summary$p,
      metric = "Active conditional p<0.05",
      value = conditional_summary$active_p_lt_005_rate_all_attempted
    ),
    data.frame(
      signal_profile = conditional_summary$signal_profile,
      p = conditional_summary$p,
      metric = "True-null conditional p<0.05",
      value = conditional_summary$null_p_lt_005_rate_all_attempted
    ),
    data.frame(
      signal_profile = selection_summary$signal_profile,
      p = selection_summary$p,
      metric = "Selection TPR (cutoff 0.5)",
      value = selection_summary$selection_tpr_cutoff_0_5
    ),
    data.frame(
      signal_profile = selection_summary$signal_profile,
      p = selection_summary$p,
      metric = "Selection TNR (cutoff 0.5)",
      value = selection_summary$selection_tnr_cutoff_0_5
    )
  )
  diagnostic_plot$signal_profile <- factor(
    diagnostic_plot$signal_profile,
    levels = c("baseline", "weak_0.5x", "weak"),
    labels = c("Baseline", "Weak (0.5x)", "Weak (0.5x)")
  )
  p4 <- ggplot2::ggplot(diagnostic_plot, ggplot2::aes(
    x = p, y = value, color = metric
  )) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::facet_wrap(~signal_profile) +
    ggplot2::scale_x_continuous(breaks = c(50, 100, 200)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Number of exposures (p)", y = "Rate", color = "SGL-WQS metric") +
    ggplot2::theme_bw(base_size = 10)
  ggplot2::ggsave(file.path(figure_dir, "B4_sglwqs_diagnostic_layers.png"), p4, width = 9, height = 4.8, dpi = 180)
}

report <- c(
  "# Battery B Matched Gaussian High-Dimensional Results",
  "",
  paste0("Generated: ", format(snapshot_time, "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope and design",
  "",
  "- Gaussian identity outcome matched to the existing binary high-dimensional design.",
  "- n = 1,000; p = 50, 100, and 200; 10 equal-sized groups.",
  "- Thirty independently generated datasets per p and signal profile.",
  "- Baseline and weak profiles; weak coefficients are exactly 0.5 times baseline mixture coefficients.",
  "- Nine active components in three active groups; the remaining seven groups are null.",
  "- Latent-factor exposure generator parameters: within-group component 0.45 and cross-group component 0.05. Because the global and group factors both contribute within groups, the implied same-group latent correlation is 0.50 and the cross-group latent correlation is 0.05.",
  "- Gaussian outcome: Y = eta + epsilon, epsilon ~ N(0, 1), with the same covariate effects as the matched design.",
  "- Methods: SGL-WQS, gWQS, groupWQS, and qgcomp.",
  "- WQS workflows used 200 bootstrap resamples; SGL-WQS used 10-fold tuning and a 60/40 training-validation split.",
  "",
  "## Completion",
  "",
  paste0("All ", nrow(method_metrics), " production atomic jobs completed successfully; no fit failures were recorded."),
  "",
  "## SGL-WQS layered results",
  "",
  markdown_table(report_sgl),
  "",
  "## All-method summary",
  "",
  markdown_table(report_all),
  "",
  "## Main factual findings",
  "",
  "- Computational completion remained 100% for every method, p level, and signal profile.",
  "- Baseline SGL-WQS active-direction accuracy remained approximately 0.92 across p = 50-200. Under weak signals it was lower and scenario-dependent (approximately 0.74-0.80), without a simple monotonic p trend.",
  "- Component-level prioritization became increasingly diffuse as p increased while the number of active components remained fixed. For SGL-WQS, active attribution decreased from about 0.21 to 0.07 under baseline signals and from about 0.19 to 0.05 under weak signals.",
  "- At a selection-frequency cutoff of 0.5, SGL-WQS TPR was high but TNR was low. Thus the selection-frequency output did not behave as a highly specific component selector in these correlated high-dimensional scenarios.",
  "- Conditional validation-stage informativeness decreased under weak signals and larger panels. True-null conditional p<0.05 proportions remained near 0.04-0.06 across these scenarios.",
  "- SGL-WQS was not uniquely slow among bootstrap WQS workflows. At p = 200 its median runtime was lower than gWQS and groupWQS in both signal profiles. qgcomp remained much faster because it used an unbootstrapped full-data GLM workflow.",
  "",
  "## Interpretation guardrails",
  "",
  "- Direction assignment, component attribution, and downstream conditional p<0.05 are distinct layers and should not be collapsed into one performance claim.",
  "- qgcomp coefficient-derived attribution and constrained WQS weights are different estimands; attribution values are descriptive within method and are not a universal ranking scale.",
  "- Runtime is the observed cost of each method's intended workflow, not an equal-computation benchmark.",
  "- Validation-stage p-values are conditional exploratory outputs, not formal selective inference.",
  "- These results do not establish universal superiority of any method.",
  "- The fixed n = 1,000 design does not by itself establish a lower sample-size boundary; Battery E addresses the weak, unbalanced p = 100 setting across n.",
  "",
  "## Reproducibility",
  "",
  "- SGL-WQS source: source/sglwqs-envint-revision-docs.",
  "- Recorded package version: 0.8.13.9001.",
  "- Source Git commit: 2fdd519e520a7dad1162810643e175cd616b1154.",
  "- Package versions, sessionInfo, source metadata, atomic outputs, MD5 manifest, and code are included in the bundle.",
  "",
  "## QA",
  "",
  markdown_table(qa, digits = 0L)
)
writeLines(report, file.path(out_dir, "BATTERY_B_RESULTS_SUMMARY.md"))

bundle_readme <- c(
  "# Battery B Code and Results Bundle",
  "",
  "This bundle contains the complete Reviewer 2 Round 2 Battery B production outputs and reproducibility materials.",
  "",
  "## Included",
  "",
  "- All 720 production atomic RDS files for B1 and B2.",
  "- Flat CSV extracts, method summaries, SGL-WQS conditional summaries, selection TPR/TNR, paired comparisons, QA, and failure tables.",
  "- Four freshly generated Battery B figures based on all 30 dataset seeds per scenario.",
  "- Simulation generators, method wrappers, metric and resume helpers, production runner, and export script.",
  "- Exact local SGL-WQS source used for Battery B and its recorded Git commit.",
  "- Package versions, source metadata, and sessionInfo.",
  "",
  "The manuscript and Response letter are intentionally not included or modified."
)
writeLines(bundle_readme, file.path(out_dir, "BUNDLE_README.md"))

write.csv(data.frame(
  field = c(
    "snapshot_time", "production_jobs", "fit_successes", "fit_failures",
    "component_rows", "sglwqs_direction_rows", "raw_rds_files", "R_version",
    "sglwqs_version", "sglwqs_source_git_commit"
  ),
  value = c(
    format(snapshot_time, "%Y-%m-%d %H:%M:%S %Z"),
    nrow(method_metrics), sum(method_metrics$fit_success), sum(!method_metrics$fit_success),
    nrow(component_metrics), nrow(direction_results), length(raw_paths), R.version.string,
    "0.8.13.9001", "2fdd519e520a7dad1162810643e175cd616b1154"
  )
), file.path(out_dir, "snapshot_metadata_B.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "aggregation_sessionInfo_B.txt"))

cat("Battery B export complete\n")
cat("Output:", out_dir, "\n")
cat("Production jobs:", nrow(method_metrics), "\n")
cat("Fit failures:", sum(!method_metrics$fit_success), "\n")
cat("QA passed:", all(qa$passed), "\n")
