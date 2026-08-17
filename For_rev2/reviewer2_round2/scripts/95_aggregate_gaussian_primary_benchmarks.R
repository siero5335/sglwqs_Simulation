if (!nzchar(Sys.getenv("R2R2_SGLWQS_SOURCE"))) {
  Sys.setenv(R2R2_SGLWQS_SOURCE = file.path(getwd(), "source", "sglwqs-envint-revision-docs"))
}
source("reviewer2_round2/R/design_helpers.R")
source_r2r2("gaussian_matched_generators.R")
source_r2r2("io_resume_helpers.R")
source_r2r2("gaussian_primary_benchmark_helpers.R")

r2r2_load_packages()
r2r2_make_dirs()

smoke_mode <- Sys.getenv("R2R2_PRIMARY_AGGREGATE_SMOKE", "false") %in% c("1", "true", "TRUE")
manifest <- if (smoke_mode) {
  make_gaussian_primary_benchmark_smoke_manifest()
} else {
  make_gaussian_primary_benchmark_manifest()
}
paths <- gaussian_primary_manifest_paths(manifest)
valid <- vapply(paths, is_complete_job_file, logical(1))
if (!any(valid)) {
  stop("No valid Gaussian primary benchmark outputs were found.", call. = FALSE)
}
results <- lapply(paths[valid], readRDS)

method_metrics <- dplyr::bind_rows(lapply(results, `[[`, "method_metrics"))
component_metrics <- dplyr::bind_rows(lapply(results, `[[`, "component_metrics"))
direction_results <- dplyr::bind_rows(lapply(results, `[[`, "sglwqs_direction_results"))

manifest_meta <- manifest |>
  dplyr::select(
    battery, scenario_id, method, data_seed,
    target_correlation = within_rho
  )
join_keys <- c("battery", "scenario_id", "method", "data_seed")
method_metrics <- dplyr::left_join(method_metrics, manifest_meta, by = join_keys)
component_metrics <- dplyr::left_join(component_metrics, manifest_meta, by = join_keys)
if (nrow(direction_results)) {
  direction_results <- dplyr::left_join(direction_results, manifest_meta, by = join_keys)
}

output_name <- if (smoke_mode) "gaussian_primary_benchmarks_smoke_aggregate" else "gaussian_primary_benchmarks"
summary_dir <- r2r2_result_file("summaries", output_name)
table_dir <- r2r2_result_file("tables", output_name)
figure_dir <- r2r2_result_file("figures", output_name)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

atomic_write_csv(method_metrics, file.path(summary_dir, "method_metrics.csv"))
atomic_write_csv(component_metrics, file.path(summary_dir, "component_metrics.csv"))
atomic_write_csv(direction_results, file.path(summary_dir, "sglwqs_direction_results.csv"))

mean_sd <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(mean = NA_real_, sd = NA_real_))
  c(mean = mean(x), sd = stats::sd(x))
}

method_summary <- method_metrics |>
  dplyr::group_by(.data$battery, .data$target_correlation, .data$method) |>
  dplyr::summarise(
    attempted = dplyr::n(),
    completed = sum(.data$fit_success),
    fit_completion_rate = mean(.data$fit_success),
    runtime_median_sec = stats::median(.data$runtime_sec, na.rm = TRUE),
    runtime_q1_sec = stats::quantile(.data$runtime_sec, 0.25, na.rm = TRUE),
    runtime_q3_sec = stats::quantile(.data$runtime_sec, 0.75, na.rm = TRUE),
    active_direction_accuracy_mean = mean(.data$active_direction_accuracy, na.rm = TRUE),
    active_direction_accuracy_sd = stats::sd(.data$active_direction_accuracy, na.rm = TRUE),
    active_attribution_mean = mean(.data$active_attribution, na.rm = TRUE),
    active_attribution_sd = stats::sd(.data$active_attribution, na.rm = TRUE),
    null_attribution_mean = mean(.data$null_attribution, na.rm = TRUE),
    null_attribution_sd = stats::sd(.data$null_attribution, na.rm = TRUE),
    bootstrap_success_mean = mean(.data$bootstrap_success_rate, na.rm = TRUE),
    .groups = "drop"
  )
atomic_write_csv(method_summary, file.path(table_dir, "01_method_completion_runtime_direction.csv"))

target_components <- component_metrics |>
  dplyr::filter(.data$group %in% c("PCBs", "Metals")) |>
  dplyr::mutate(
    target_direction_weight = dplyr::if_else(
      .data$group == "PCBs",
      .data$positive_weight,
      .data$negative_weight
    ),
    threshold_selected = .data$target_direction_weight > 0.10
  )

group_seed_metrics <- target_components |>
  dplyr::group_by(
    .data$battery, .data$target_correlation, .data$method,
    .data$data_seed, .data$group, .data$attribution_type
  ) |>
  dplyr::mutate(
    target_rank = rank(-.data$target_direction_weight, ties.method = "average"),
    normalized_rank = .data$target_rank / dplyr::n()
  ) |>
  dplyr::summarise(
    active_component_count = sum(.data$is_active),
    active_share = {
      denom <- sum(.data$target_direction_weight, na.rm = TRUE)
      if (denom > 0) sum(.data$target_direction_weight[.data$is_active], na.rm = TRUE) / denom else 0
    },
    null_share = 1 - active_share,
    pragmatic_tpr = mean(.data$threshold_selected[.data$is_active]),
    pragmatic_tnr = mean(!.data$threshold_selected[!.data$is_active]),
    active_normalized_rank = mean(.data$normalized_rank[.data$is_active]),
    .groups = "drop"
  )
atomic_write_csv(group_seed_metrics, file.path(summary_dir, "group_seed_metrics.csv"))

group_summary <- group_seed_metrics |>
  dplyr::group_by(
    .data$battery, .data$target_correlation, .data$method,
    .data$group, .data$attribution_type
  ) |>
  dplyr::summarise(
    seeds = dplyr::n(),
    active_share_mean = mean(.data$active_share, na.rm = TRUE),
    active_share_sd = stats::sd(.data$active_share, na.rm = TRUE),
    null_share_mean = mean(.data$null_share, na.rm = TRUE),
    null_share_sd = stats::sd(.data$null_share, na.rm = TRUE),
    pragmatic_tpr_mean = mean(.data$pragmatic_tpr, na.rm = TRUE),
    pragmatic_tpr_sd = stats::sd(.data$pragmatic_tpr, na.rm = TRUE),
    pragmatic_tnr_mean = mean(.data$pragmatic_tnr, na.rm = TRUE),
    pragmatic_tnr_sd = stats::sd(.data$pragmatic_tnr, na.rm = TRUE),
    active_normalized_rank_mean = mean(.data$active_normalized_rank, na.rm = TRUE),
    active_normalized_rank_sd = stats::sd(.data$active_normalized_rank, na.rm = TRUE),
    .groups = "drop"
  )
atomic_write_csv(group_summary, file.path(table_dir, "02_group_active_null_attribution_tpr_tnr.csv"))

sgl_selection <- target_components |>
  dplyr::filter(.data$method == "SGL-WQS") |>
  dplyr::mutate(stably_selected = .data$selection_frequency > 0.50) |>
  dplyr::group_by(.data$battery, .data$target_correlation, .data$data_seed) |>
  dplyr::summarise(
    selection_tpr = mean(.data$stably_selected[.data$is_active]),
    selection_tnr = mean(!.data$stably_selected[!.data$is_active]),
    active_selection_frequency = mean(.data$selection_frequency[.data$is_active], na.rm = TRUE),
    null_selection_frequency = mean(.data$selection_frequency[!.data$is_active], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(.data$battery, .data$target_correlation) |>
  dplyr::summarise(
    seeds = dplyr::n(),
    selection_tpr_mean = mean(.data$selection_tpr),
    selection_tpr_sd = stats::sd(.data$selection_tpr),
    selection_tnr_mean = mean(.data$selection_tnr),
    selection_tnr_sd = stats::sd(.data$selection_tnr),
    active_selection_frequency_mean = mean(.data$active_selection_frequency),
    null_selection_frequency_mean = mean(.data$null_selection_frequency),
    .groups = "drop"
  )
atomic_write_csv(sgl_selection, file.path(table_dir, "03_sglwqs_selection_frequency_tpr_tnr.csv"))

if (nrow(direction_results)) {
  conditional_summary <- direction_results |>
    dplyr::mutate(retained_flag = .data$retained) |>
    dplyr::group_by(
      .data$battery, .data$target_correlation, .data$group,
      .data$direction, .data$is_true_null
    ) |>
    dplyr::summarise(
      attempted = dplyr::n(),
      retained = sum(.data$retained_flag),
      retention_rate = mean(.data$retained_flag),
      conditional_p_lt_0_05_retained = ifelse(
        sum(.data$retained_flag) > 0,
        mean(dplyr::coalesce(.data$rejected_0_05[.data$retained_flag], FALSE)),
        NA_real_
      ),
      conditional_p_lt_0_05_all_attempted = mean(dplyr::coalesce(.data$rejected_0_05, FALSE)),
      estimate_mean_retained = mean(.data$estimate[.data$retained_flag], na.rm = TRUE),
      estimate_sd_retained = stats::sd(.data$estimate[.data$retained_flag], na.rm = TRUE),
      minor_direction_exclusion_rate = mean(.data$excluded_by_minor_threshold),
      bootstrap_success_mean = mean(.data$boot_success_rate, na.rm = TRUE),
      .groups = "drop"
    )
  atomic_write_csv(conditional_summary, file.path(table_dir, "04_sglwqs_conditional_retention_pvalues.csv"))
} else {
  conditional_summary <- data.frame()
}

failures <- method_metrics |>
  dplyr::filter(!.data$fit_success) |>
  dplyr::select(
    battery, scenario_id, target_correlation,
    data_seed, method, failure_stage,
    error_class, error_message, warning_summary
  )
atomic_write_csv(failures, file.path(table_dir, "05_failures.csv"))

correlation_plot_data <- group_summary |>
  dplyr::filter(grepl("gaussian_correlation_robustness", .data$battery, fixed = TRUE)) |>
  dplyr::mutate(active_share_sd = dplyr::coalesce(.data$active_share_sd, 0))
p1 <- ggplot2::ggplot(
  correlation_plot_data,
  ggplot2::aes(
    x = .data$target_correlation,
    y = .data$active_share_mean,
    color = .data$method,
    group = .data$method
  )
) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 2) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = pmax(0, .data$active_share_mean - .data$active_share_sd),
      ymax = pmin(1, .data$active_share_mean + .data$active_share_sd)
    ),
    width = 0.025,
    linewidth = 0.4
  ) +
  ggplot2::facet_wrap(~group) +
  ggplot2::scale_x_continuous(breaks = c(0.2, 0.5, 0.8, 0.95)) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(
    x = "Target within-group correlation",
    y = "Active component attribution share",
    color = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "bottom")
ggplot2::ggsave(file.path(figure_dir, "01_correlation_active_share.png"), p1, width = 8, height = 4.8, dpi = 200)
ggplot2::ggsave(file.path(figure_dir, "01_correlation_active_share.pdf"), p1, width = 8, height = 4.8)

active_plot_data <- group_seed_metrics |>
  dplyr::filter(grepl("gaussian_active_component", .data$battery, fixed = TRUE))
p2 <- ggplot2::ggplot(
  active_plot_data,
  ggplot2::aes(x = .data$method, y = .data$active_share, fill = .data$method)
) +
  ggplot2::geom_violin(alpha = 0.25, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.6) +
  ggplot2::facet_wrap(~group) +
  ggplot2::scale_y_continuous(limits = c(0, 1)) +
  ggplot2::labs(x = NULL, y = "Attribution share assigned to the true driver") +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
ggplot2::ggsave(file.path(figure_dir, "02_active_component_share.png"), p2, width = 8, height = 4.8, dpi = 200)
ggplot2::ggsave(file.path(figure_dir, "02_active_component_share.pdf"), p2, width = 8, height = 4.8)

source_commit <- git_commit_for_path(r2r2_sglwqs_source_path())
design_n <- paste(sort(unique(manifest$n)), collapse = ", ")
design_correlations <- paste(sort(unique(manifest$within_rho[grepl("correlation", manifest$battery)])), collapse = ", ")
design_seeds <- length(unique(manifest$data_seed))
report <- c(
  "# Gaussian Counterparts of the Primary Logistic Benchmarks",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Design",
  "",
  sprintf("- Correlation robustness: n = %s; r = %s; %d dataset seed(s); four methods.", design_n, design_correlations, design_seeds),
  sprintf("- Active-component prioritization: n = %s; PCB-153 and Pb as the only active components in their groups; %d dataset seed(s); four methods.", design_n, design_seeds),
  "- The published logistic X generation, covariates, coefficients, and method settings were retained. Only the outcome was replaced by y = eta + epsilon, epsilon ~ N(0, 1).",
  "- qgcomp summaries are coefficient-derived attribution, whereas WQS-type summaries are constrained index weights; they are not treated as identical estimands.",
  "",
  "## Completion",
  "",
  sprintf("- Valid atomic jobs: %d/%d", sum(valid), length(valid)),
  sprintf("- Successful fits: %d/%d", sum(method_metrics$fit_success), nrow(method_metrics)),
  sprintf("- Recorded fit failures: %d", nrow(failures)),
  sprintf("- SGL-WQS commit: `%s`", source_commit),
  "",
  "## Outputs",
  "",
  "- `01_method_completion_runtime_direction.csv`",
  "- `02_group_active_null_attribution_tpr_tnr.csv`",
  "- `03_sglwqs_selection_frequency_tpr_tnr.csv`",
  "- `04_sglwqs_conditional_retention_pvalues.csv`",
  "- `05_failures.csv`",
  "",
  "The numerical tables intentionally separate direction assignment, component prioritization, computational completion, and exploratory conditional validation-stage p-values."
)
writeLines(report, file.path(summary_dir, "GAUSSIAN_PRIMARY_BENCHMARK_REPORT.md"))
cat(paste(report, collapse = "\n"), "\n")
