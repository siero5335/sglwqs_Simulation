if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("unbalanced_factorial_generator.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("metric_helpers.R")
source_r2r2("dense_group_signal_generator.R")
source_r2r2("dense_group_signal_helpers.R")

r2r2_load_packages()
r2r2_make_dirs()

summary_dir <- r2r2_result_file("summaries", "dense_group_signal")
table_dir <- r2r2_result_file("tables", "dense_group_signal")
figure_dir <- r2r2_result_file("figures", "dense_group_signal")
invisible(lapply(c(summary_dir, table_dir, figure_dir), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

manifest <- make_dense_group_signal_manifest()
h_paths <- dense_group_signal_manifest_paths(manifest)
h_valid <- vapply(h_paths, is_complete_job_file, logical(1))
h_results <- atomic_read_results(h_paths[h_valid])

e_root <- r2r2_result_file("raw", "E_hard_setting_sample_size")
e_paths <- if (dir.exists(e_root)) {
  list.files(e_root, pattern = "^seed.*[.]rds$", recursive = TRUE, full.names = TRUE)
} else {
  character(0)
}
e_paths <- e_paths[grepl("E_(binomial|gaussian)_unbalanced_weak_heterogeneous_n(500|5000)_p100", e_paths)]
e_results <- atomic_read_results(e_paths)

h_method <- bind_result_slot(h_results, "method_metrics")
h_component <- correct_component_truth_support(bind_result_slot(h_results, "component_metrics"))
h_direction <- correct_sglwqs_direction_truth_support(
  bind_result_slot(h_results, "sglwqs_direction_results"), h_component
)
e_method <- bind_result_slot(e_results, "method_metrics")

h_summary <- summarize_method_performance(h_method)
h_component_summary <- summarize_component_strata(h_component)
h_validation <- summarize_sglwqs_scenario_validation(h_direction)
h_validation_detail <- summarize_sglwqs_validation(h_direction)
h_failures <- h_method |>
  dplyr::filter(!.data$fit_success) |>
  dplyr::select(dplyr::any_of(c(
    "scenario_id", "family", "n", "method", "data_seed", "failure_stage",
    "error_class", "error_message", "warning_summary", "runtime_sec"
  )))

e_summary <- summarize_method_performance(e_method)
h_summary$truth_profile <- "dense_group_matched"
e_summary$truth_profile <- "sparse_weak_heterogeneous"
matched_summary <- dplyr::bind_rows(e_summary, h_summary)

paired <- h_method |>
  dplyr::select(
    .data$family, .data$n, .data$method, .data$data_seed,
    dense_direction = .data$active_direction_accuracy,
    dense_active_attribution = .data$active_attribution,
    dense_null_attribution = .data$null_attribution,
    dense_runtime = .data$runtime_sec
  ) |>
  dplyr::inner_join(
    e_method |>
      dplyr::select(
        .data$family, .data$n, .data$method, .data$data_seed,
        sparse_direction = .data$active_direction_accuracy,
        sparse_active_attribution = .data$active_attribution,
        sparse_null_attribution = .data$null_attribution,
        sparse_runtime = .data$runtime_sec
      ),
    by = c("family", "n", "method", "data_seed")
  ) |>
  dplyr::mutate(
    direction_difference = .data$dense_direction - .data$sparse_direction,
    active_attribution_difference = .data$dense_active_attribution - .data$sparse_active_attribution,
    null_attribution_difference = .data$dense_null_attribution - .data$sparse_null_attribution,
    runtime_difference = .data$dense_runtime - .data$sparse_runtime
  )

paired_summary <- paired |>
  dplyr::group_by(.data$family, .data$n, .data$method) |>
  dplyr::summarise(
    paired_seeds = dplyr::n(),
    direction_difference_mean = mean(.data$direction_difference, na.rm = TRUE),
    direction_difference_sd = stats::sd(.data$direction_difference, na.rm = TRUE),
    active_attribution_difference_mean = mean(.data$active_attribution_difference, na.rm = TRUE),
    active_attribution_difference_sd = stats::sd(.data$active_attribution_difference, na.rm = TRUE),
    null_attribution_difference_mean = mean(.data$null_attribution_difference, na.rm = TRUE),
    runtime_difference_median_sec = stats::median(.data$runtime_difference, na.rm = TRUE),
    .groups = "drop"
  )

signal_spec <- make_dense_group_signal_truth()$signal_spec
completion <- data.frame(
  expected_jobs = nrow(manifest),
  valid_jobs = sum(h_valid),
  remaining_jobs = sum(!h_valid),
  fit_successes = sum(h_method$fit_success, na.rm = TRUE),
  recorded_fit_failures = sum(!h_method$fit_success, na.rm = TRUE),
  stringsAsFactors = FALSE
)

readr::write_csv(h_summary, file.path(table_dir, "01_dense_method_summary.csv"))
readr::write_csv(matched_summary, file.path(table_dir, "02_sparse_vs_dense_summary.csv"))
readr::write_csv(h_component_summary, file.path(table_dir, "03_dense_component_strata.csv"))
readr::write_csv(h_validation, file.path(table_dir, "04_dense_sglwqs_validation.csv"))
readr::write_csv(h_validation_detail, file.path(table_dir, "05_dense_sglwqs_validation_detail.csv"))
readr::write_csv(h_failures, file.path(table_dir, "06_dense_failures.csv"))
readr::write_csv(paired, file.path(table_dir, "07_sparse_vs_dense_paired_seed_results.csv"))
readr::write_csv(paired_summary, file.path(table_dir, "08_sparse_vs_dense_paired_summary.csv"))
readr::write_csv(signal_spec, file.path(table_dir, "09_dense_signal_specification.csv"))
readr::write_csv(completion, file.path(table_dir, "10_dense_completion.csv"))

save_plot <- function(plot, name, width = 10, height = 6) {
  ggplot2::ggsave(file.path(figure_dir, paste0(name, ".png")), plot,
    width = width, height = height, dpi = 180
  )
}

if (nrow(matched_summary)) {
  p_direction <- ggplot2::ggplot(
    matched_summary,
    ggplot2::aes(x = factor(.data$n), y = .data$active_direction_accuracy_mean,
      color = .data$method, group = .data$method
    )
  ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.2)) +
    ggplot2::geom_line(position = ggplot2::position_dodge(width = 0.2)) +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymin = pmax(0, .data$active_direction_accuracy_mean - .data$active_direction_accuracy_sd),
      ymax = pmin(1, .data$active_direction_accuracy_mean + .data$active_direction_accuracy_sd)
    ), width = 0.1, position = ggplot2::position_dodge(width = 0.2)) +
    ggplot2::facet_grid(.data$family ~ .data$truth_profile) +
    ggplot2::labs(x = "n", y = "Active-direction accuracy", color = "Method") +
    ggplot2::theme_bw()
  save_plot(p_direction, "01_sparse_vs_dense_direction_accuracy")

  p_attribution <- ggplot2::ggplot(
    matched_summary,
    ggplot2::aes(x = factor(.data$n), y = .data$active_attribution_mean,
      color = .data$method, group = .data$method
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_grid(.data$family ~ .data$truth_profile) +
    ggplot2::labs(x = "n", y = "Active attribution", color = "Method") +
    ggplot2::theme_bw()
  save_plot(p_attribution, "02_sparse_vs_dense_active_attribution")

  p_runtime <- ggplot2::ggplot(
    matched_summary,
    ggplot2::aes(x = factor(.data$n), y = .data$runtime_median_sec,
      color = .data$method, group = .data$method
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::scale_y_log10() +
    ggplot2::facet_grid(.data$family ~ .data$truth_profile) +
    ggplot2::labs(x = "n", y = "Median runtime (seconds, log scale)", color = "Method") +
    ggplot2::theme_bw()
  save_plot(p_runtime, "03_sparse_vs_dense_runtime")
}

if (nrow(h_validation)) {
  validation_long <- h_validation |>
    tidyr::pivot_longer(
      cols = c(
        "active_conditional_detection_all", "active_conditional_detection_retained",
        "null_rejection_all", "null_conditional_rejection_retained"
      ),
      names_to = "metric", values_to = "rate"
    )
  p_validation <- ggplot2::ggplot(
    validation_long,
    ggplot2::aes(x = factor(.data$n), y = .data$rate, color = .data$metric, group = .data$metric)
  ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~family) +
    ggplot2::labs(x = "n", y = "Conditional validation-stage proportion", color = "Metric") +
    ggplot2::theme_bw()
  save_plot(p_validation, "04_dense_sglwqs_conditional_validation")
}

report_lines <- c(
  "# Battery H: Matched Dense Group-Signal Report",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  "Battery H is a matched sensitivity analysis for Battery E. It uses the same p=100 unbalanced groups, correlation targets, covariates, data seeds, outcome families, sample sizes, train-validation splits, and method settings. The only intended DGP change is redistribution of each active group's signal from three sparse components across all components in that group while preserving the reference group-level signal variance.",
  "",
  "- Sample sizes: n=500 and n=5,000",
  "- Families: binomial-logistic and Gaussian-identity",
  "- Replicates: 30 independent seeds per scenario",
  "- Methods: SGL-WQS, gWQS, groupWQS, and qgcomp",
  "- Active components: 30/100 across G01, G05, and G10",
  "- Gaussian residual SD: 1",
  "",
  "## Completion",
  "",
  paste(capture.output(knitr::kable(completion, format = "pipe")), collapse = "\n"),
  "",
  "## Dense Signal Specification",
  "",
  paste(capture.output(knitr::kable(signal_spec, digits = 4, format = "pipe")), collapse = "\n"),
  "",
  "## Method Summary",
  "",
  if (nrow(h_summary)) paste(capture.output(knitr::kable(
    h_summary |>
      dplyr::select(
        .data$family, .data$n, .data$method, .data$attempted, .data$completed,
        .data$active_direction_accuracy_mean, .data$active_attribution_mean,
        .data$null_attribution_mean, .data$runtime_median_sec
      ),
    digits = 3, format = "pipe"
  )), collapse = "\n") else "No completed Battery H jobs are currently available.",
  "",
  "## Interpretation Guardrails",
  "",
  "- Active-direction accuracy and component-level attribution describe different performance layers.",
  "- qgcomp attribution is coefficient-derived and is not a WQS-type constrained weight.",
  "- SGL-WQS validation-stage p-values are conditional/exploratory and are not formal selective inference.",
  "- Dense and sparse truth profiles should be compared within the same family, n, method, and data seed.",
  "- Incomplete results must not be interpreted as final production summaries."
)
writeLines(report_lines, file.path(summary_dir, "BATTERY_H_DENSE_GROUP_SIGNAL_REPORT.md"))

cat(sprintf(
  "Dense group-signal aggregation complete: %d/%d valid jobs; %d recorded fit failures.\n",
  sum(h_valid), nrow(manifest), sum(!h_method$fit_success, na.rm = TRUE)
))
