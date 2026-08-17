read_all_job_results <- function() {
  paths <- list.files(r2r2_result_file("raw"), pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  paths <- paths[!grepl("checkpoint_", paths, fixed = TRUE)]
  atomic_read_results(paths)
}

bind_result_slot <- function(results, slot) {
  out <- lapply(results, `[[`, slot)
  out <- out[vapply(out, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(out)) return(data.frame())
  dplyr::bind_rows(out)
}

iqr_text <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf("%.3g [%.3g, %.3g]", q[[2]], q[[1]], q[[3]])
}

mean_sd_text <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  sprintf("%.3f (%.3f)", mean(x), stats::sd(x))
}

summarize_method_performance <- function(method_metrics) {
  if (!nrow(method_metrics)) return(data.frame())
  method_metrics |>
    dplyr::group_by(
      .data$battery, .data$scenario_id, .data$family, .data$n, .data$p,
      .data$signal_profile, .data$effect_profile, .data$group_structure, .data$method
    ) |>
    dplyr::summarise(
      attempted = dplyr::n(),
      completed = sum(.data$fit_success, na.rm = TRUE),
      fit_completion_rate = .data$completed / .data$attempted,
      failure_rate = 1 - .data$fit_completion_rate,
      runtime_median_sec = stats::median(.data$runtime_sec, na.rm = TRUE),
      runtime_iqr_sec = stats::IQR(.data$runtime_sec, na.rm = TRUE),
      runtime_median_iqr = iqr_text(.data$runtime_sec),
      active_direction_accuracy_mean = mean(.data$active_direction_accuracy, na.rm = TRUE),
      active_direction_accuracy_sd = stats::sd(.data$active_direction_accuracy, na.rm = TRUE),
      sign_assignment_accuracy_mean = mean(.data$sign_assignment_accuracy, na.rm = TRUE),
      sign_assignment_accuracy_sd = stats::sd(.data$sign_assignment_accuracy, na.rm = TRUE),
      active_attribution_mean = mean(.data$active_attribution, na.rm = TRUE),
      active_attribution_sd = stats::sd(.data$active_attribution, na.rm = TRUE),
      null_attribution_mean = mean(.data$null_attribution, na.rm = TRUE),
      null_attribution_sd = stats::sd(.data$null_attribution, na.rm = TRUE),
      active_selection_frequency_mean = mean(.data$active_selection_frequency, na.rm = TRUE),
      null_selection_frequency_mean = mean(.data$null_selection_frequency, na.rm = TRUE),
      bootstrap_success_rate_median = stats::median(.data$bootstrap_success_rate, na.rm = TRUE),
      retained_group_directions_median = stats::median(.data$retained_group_directions, na.rm = TRUE),
      between_seed_sd_active_direction = stats::sd(.data$active_direction_accuracy, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_sglwqs_validation <- function(dir_rows) {
  if (!nrow(dir_rows)) return(data.frame())
  dir_rows |>
    dplyr::group_by(
      .data$battery, .data$scenario_id, .data$family, .data$n, .data$p,
      .data$group, .data$direction, .data$is_true_null
    ) |>
    dplyr::summarise(
      attempted_group_directions = dplyr::n(),
      retained_n = sum(.data$retained, na.rm = TRUE),
      retained_rate = .data$retained_n / .data$attempted_group_directions,
      excluded_rate = mean(.data$excluded_by_minor_threshold, na.rm = TRUE),
      conditional_rejection_n = sum(.data$retained & .data$rejected_0_05, na.rm = TRUE),
      conditional_rejection_rate_retained = .data$conditional_rejection_n / .data$retained_n,
      rejection_rate_all_attempted = sum(.data$rejected_0_05, na.rm = TRUE) / .data$attempted_group_directions,
      estimate_mean = mean(.data$estimate, na.rm = TRUE),
      estimate_sd = stats::sd(.data$estimate, na.rm = TRUE),
      p_median = stats::median(.data$p_value, na.rm = TRUE),
      p_iqr = stats::IQR(.data$p_value, na.rm = TRUE),
      p_min = safe_min(.data$p_value),
      p_max = safe_max(.data$p_value),
      boot_success_rate_median = stats::median(.data$boot_success_rate, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_sglwqs_scenario_validation <- function(dir_rows) {
  if (!nrow(dir_rows)) return(data.frame())
  dir_rows |>
    dplyr::group_by(.data$battery, .data$scenario_id, .data$family, .data$n, .data$p) |>
    dplyr::summarise(
      group_direction_rows = dplyr::n(),
      retained_rate = mean(.data$retained, na.rm = TRUE),
      active_conditional_detection_retained = {
        active <- !.data$is_true_null
        mean(.data$rejected_0_05[active & .data$retained], na.rm = TRUE)
      },
      active_conditional_detection_all = {
        active <- !.data$is_true_null
        mean(.data$rejected_0_05[active], na.rm = TRUE)
      },
      null_conditional_rejection_retained = {
        null <- .data$is_true_null
        mean(.data$rejected_0_05[null & .data$retained], na.rm = TRUE)
      },
      null_rejection_all = {
        null <- .data$is_true_null
        mean(.data$rejected_0_05[null], na.rm = TRUE)
      },
      minor_direction_exclusion_rate = mean(.data$excluded_by_minor_threshold, na.rm = TRUE),
      .groups = "drop"
    )
}

summarize_component_strata <- function(component_metrics) {
  if (!nrow(component_metrics)) return(data.frame())
  component_metrics |>
    dplyr::group_by(
      .data$battery, .data$scenario_id, .data$family, .data$n, .data$p,
      .data$method, .data$group_size_tier, .data$effect_tier,
      .data$true_direction, .data$group_type
    ) |>
    dplyr::summarise(
      rows = dplyr::n(),
      active_rows = sum(.data$is_active, na.rm = TRUE),
      direction_assignment_accuracy = mean(.data$direction_correct[.data$is_active], na.rm = TRUE),
      normalized_attribution_mean = mean(.data$normalized_attribution, na.rm = TRUE),
      selection_frequency_mean = mean(.data$selection_frequency, na.rm = TRUE),
      group_omission_rate = mean(.data$is_active & (.data$normalized_attribution <= 0 | is.na(.data$normalized_attribution)), na.rm = TRUE),
      .groups = "drop"
    )
}

paired_boot_ci <- function(diff, conf = 0.95, n_boot = 2000L) {
  diff <- diff[is.finite(diff)]
  if (length(diff) < 2L) {
    return(c(estimate = mean(diff), lower = NA_real_, upper = NA_real_, p_value = NA_real_))
  }
  set.seed(61035)
  b <- replicate(n_boot, mean(sample(diff, length(diff), replace = TRUE)))
  alpha <- (1 - conf) / 2
  p <- stats::wilcox.test(diff, mu = 0, exact = FALSE)$p.value
  c(
    estimate = mean(diff),
    lower = unname(stats::quantile(b, alpha, na.rm = TRUE)),
    upper = unname(stats::quantile(b, 1 - alpha, na.rm = TRUE)),
    p_value = p
  )
}

paired_method_comparisons <- function(method_metrics,
                                      metrics = c("active_direction_accuracy", "active_attribution", "runtime_sec"),
                                      reference = "SGL-WQS") {
  if (!nrow(method_metrics)) return(data.frame())
  out <- list()
  idx <- 1L
  keys <- c("battery", "scenario_id", "family", "n", "p", "data_seed", "split_replicate")
  for (metric in metrics) {
    wide <- method_metrics |>
      dplyr::select(dplyr::all_of(c(keys, "method", metric))) |>
      tidyr::pivot_wider(names_from = .data$method, values_from = dplyr::all_of(metric))
    ref_col <- reference
    if (!ref_col %in% names(wide)) next
    methods <- setdiff(unique(method_metrics$method), reference)
    for (m in methods) {
      if (!m %in% names(wide)) next
      dif <- wide[[ref_col]] - wide[[m]]
      ci <- paired_boot_ci(dif)
      base <- wide |>
        dplyr::distinct(.data$battery, .data$scenario_id, .data$family, .data$n, .data$p)
      for (i in seq_len(nrow(base))) {
        out[[idx]] <- data.frame(
          base[i, , drop = FALSE],
          metric = metric,
          reference_method = reference,
          comparator_method = m,
          mean_paired_difference = ci[["estimate"]],
          bootstrap_ci_low = ci[["lower"]],
          bootstrap_ci_high = ci[["upper"]],
          raw_p_value = ci[["p_value"]],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  res <- dplyr::bind_rows(out)
  if (!nrow(res)) return(res)
  res |>
    dplyr::group_by(.data$battery, .data$family, .data$metric) |>
    dplyr::mutate(holm_p_value = stats::p.adjust(.data$raw_p_value, method = "holm")) |>
    dplyr::ungroup()
}

summarize_split_stability <- function(dir_rows) {
  split <- dir_rows[dir_rows$battery == "C_gaussian_split_stability", , drop = FALSE]
  if (!nrow(split)) return(data.frame())
  split |>
    dplyr::group_by(
      .data$battery, .data$scenario_id, .data$family, .data$n,
      .data$group, .data$direction, .data$is_true_null
    ) |>
    dplyr::summarise(
      splits = dplyr::n(),
      estimable_splits = sum(.data$retained & is.finite(.data$p_value), na.rm = TRUE),
      retention_rate = mean(.data$retained, na.rm = TRUE),
      conditional_rejection_rate = mean(.data$rejected_0_05[.data$retained], na.rm = TRUE),
      estimate_mean = mean(.data$estimate, na.rm = TRUE),
      estimate_sd = stats::sd(.data$estimate, na.rm = TRUE),
      positive_sign_rate = mean(.data$sign == "positive", na.rm = TRUE),
      negative_sign_rate = mean(.data$sign == "negative", na.rm = TRUE),
      p_value_median = stats::median(.data$p_value, na.rm = TRUE),
      p_value_iqr = stats::IQR(.data$p_value, na.rm = TRUE),
      p_value_min = safe_min(.data$p_value),
      p_value_max = safe_max(.data$p_value),
      p_value_cv = stats::sd(.data$p_value, na.rm = TRUE) / mean(.data$p_value, na.rm = TRUE),
      sd_log10_p = stats::sd(log10(.data$p_value[is.finite(.data$p_value) & .data$p_value > 0]), na.rm = TRUE),
      model_se_mean = mean(.data$std_error, na.rm = TRUE),
      empirical_sd = stats::sd(.data$estimate, na.rm = TRUE),
      model_se_to_empirical_sd = mean(.data$std_error, na.rm = TRUE) / stats::sd(.data$estimate, na.rm = TRUE),
      boot_success_rate_median = stats::median(.data$boot_success_rate, na.rm = TRUE),
      fit_failures = sum(!is.finite(.data$p_value) & !.data$retained, na.rm = TRUE),
      .groups = "drop"
    )
}

make_tiered_usage_evidence <- function(method_summary, sgl_summary) {
  if (!nrow(method_summary)) return(data.frame())
  sgl <- method_summary |>
    dplyr::filter(.data$method == "SGL-WQS") |>
    dplyr::mutate(
      sample_size_tier = dplyr::case_when(
        .data$n < 500 ~ "very small",
        .data$n == 500 ~ "n=500",
        .data$n <= 1000 ~ "n=1,000",
        .data$n <= 5000 ~ "n=5,000",
        TRUE ~ "large"
      ),
      dimensional_tier = dplyr::case_when(
        .data$p <= 13 ~ "low",
        .data$p == 50 ~ "moderate",
        .data$p == 100 ~ "medium-high",
        .data$p >= 200 ~ "high",
        TRUE ~ "unclassified"
      ),
      signal_tier = dplyr::case_when(
        grepl("weak", .data$signal_profile) ~ "weak",
        grepl("heterogeneous", .data$effect_profile) ~ "heterogeneous",
        TRUE ~ "baseline"
      ),
      group_balance_tier = dplyr::case_when(
        .data$group_structure == "unbalanced" ~ "unbalanced",
        .data$group_structure %in% c("balanced", "equal10") ~ "balanced/equal",
        TRUE ~ .data$group_structure
      )
    )
  needed_sgl_cols <- c(
    "battery", "scenario_id", "family", "n", "p",
    "active_conditional_detection_retained", "null_rejection_all"
  )
  if (is.data.frame(sgl_summary) && all(needed_sgl_cols %in% names(sgl_summary))) {
    sgl_val <- sgl_summary |>
      dplyr::select(
        battery, scenario_id, family, n, p,
        active_conditional_detection_retained, null_rejection_all
      )
  } else {
    sgl_val <- data.frame(
      battery = character(0),
      scenario_id = character(0),
      family = character(0),
      n = integer(0),
      p = integer(0),
      active_conditional_detection_retained = numeric(0),
      null_rejection_all = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  sgl |>
    dplyr::left_join(sgl_val, by = c("battery", "scenario_id", "family", "n", "p")) |>
    dplyr::transmute(
      family = .data$family,
      sample_size_tier = .data$sample_size_tier,
      dimensional_tier = .data$dimensional_tier,
      signal_tier = .data$signal_tier,
      group_balance_tier = .data$group_balance_tier,
      direction_accuracy = .data$active_direction_accuracy_mean,
      active_conditional_detection = .data$active_conditional_detection_retained,
      null_conditional_rejection = .data$null_rejection_all,
      component_attribution = .data$active_attribution_mean,
      split_stability = NA_character_,
      failure_rate = .data$failure_rate,
      runtime = .data$runtime_median_iqr,
      evidence_supported_caution = dplyr::case_when(
        .data$failure_rate > 0.1 ~ "computational completion requires caution",
        .data$active_direction_accuracy_mean < 0.7 ~ "direction assignment was scenario-dependent",
        .data$active_attribution_mean < 0.5 ~ "component-level attribution became less concentrated",
        TRUE ~ "no threshold recommendation assigned"
      )
    )
}

write_markdown_table <- function(x, path, caption = NULL, digits = 3) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  lines <- capture.output(knitr::kable(x, format = "markdown", digits = digits, caption = caption))
  writeLines(lines, path)
  invisible(path)
}

write_ggplot <- function(plot, name, width = 8, height = 5) {
  dir.create(r2r2_result_file("figures"), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(r2r2_result_file("figures", paste0(name, ".png")), plot,
    width = width, height = height, dpi = 150
  )
  ggplot2::ggsave(r2r2_result_file("figures", paste0(name, ".pdf")), plot,
    width = width, height = height
  )
}

generate_core_figures <- function(method_summary, sgl_summary, component_strata) {
  if (!nrow(method_summary)) return(invisible(NULL))
  sample_df <- method_summary |>
    dplyr::filter(grepl("sample_size", .data$battery))
  if (nrow(sample_df)) {
    p1 <- ggplot2::ggplot(sample_df, ggplot2::aes(x = .data$n, y = .data$active_direction_accuracy_mean, color = .data$method)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(family ~ signal_profile) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "n", y = "Active-direction accuracy")
    write_ggplot(p1, "01_active_direction_accuracy_vs_n")
  }
  if (nrow(sgl_summary)) {
    s1 <- sgl_summary |> dplyr::filter(grepl("sample_size|hard_setting", .data$battery))
    if (nrow(s1)) {
      p2 <- ggplot2::ggplot(s1, ggplot2::aes(x = .data$n, y = .data$active_conditional_detection_retained, color = .data$family)) +
        ggplot2::geom_line() +
        ggplot2::geom_point() +
        ggplot2::scale_x_log10() +
        ggplot2::labs(x = "n", y = "SGL-WQS conditional p < 0.05 among retained active directions")
      write_ggplot(p2, "02_sglwqs_conditional_detection_vs_n")
    }
  }
  high_df <- method_summary |> dplyr::filter(grepl("highdim", .data$battery))
  if (nrow(high_df)) {
    p3 <- ggplot2::ggplot(high_df, ggplot2::aes(x = .data$p, y = .data$active_direction_accuracy_mean, color = .data$method)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(family ~ signal_profile) +
      ggplot2::labs(x = "p", y = "Active-direction accuracy")
    write_ggplot(p3, "03_active_direction_accuracy_vs_p")

    p4 <- ggplot2::ggplot(high_df, ggplot2::aes(x = .data$p, y = .data$active_attribution_mean, color = .data$method)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(family ~ signal_profile) +
      ggplot2::labs(x = "p", y = "Active attribution")
    write_ggplot(p4, "04_component_attribution_vs_p")

    p5 <- ggplot2::ggplot(high_df, ggplot2::aes(x = .data$p, y = .data$runtime_median_sec, color = .data$method)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_grid(family ~ signal_profile) +
      ggplot2::labs(x = "p", y = "Median runtime, seconds")
    write_ggplot(p5, "05_runtime_vs_p")
  }
  fac <- component_strata |>
    dplyr::filter(.data$battery == "D_unbalanced_effect_factorial", .data$method == "SGL-WQS", .data$effect_tier != "null")
  if (nrow(fac)) {
    p7 <- ggplot2::ggplot(fac, ggplot2::aes(x = .data$effect_tier, y = .data$group_size_tier, fill = .data$direction_assignment_accuracy)) +
      ggplot2::geom_tile() +
      ggplot2::facet_grid(family ~ true_direction) +
      ggplot2::scale_fill_viridis_c(na.value = "grey80") +
      ggplot2::labs(x = "Effect tier", y = "Group size tier", fill = "Direction accuracy")
    write_ggplot(p7, "07_group_size_effect_tier_heatmap")
  }
  hard <- method_summary |> dplyr::filter(.data$battery == "E_hard_setting_sample_size")
  if (nrow(hard)) {
    p8 <- ggplot2::ggplot(hard, ggplot2::aes(x = .data$n, y = .data$active_direction_accuracy_mean, color = .data$method)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~family) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "n", y = "Hard-setting active-direction accuracy")
    write_ggplot(p8, "08_hard_setting_performance_vs_n")
  }
  completion <- method_summary
  if (nrow(completion)) {
    p10 <- ggplot2::ggplot(completion, ggplot2::aes(x = .data$method, y = .data$fit_completion_rate, fill = .data$family)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::facet_wrap(~battery) +
      ggplot2::coord_flip() +
      ggplot2::labs(x = NULL, y = "Fit completion rate")
    write_ggplot(p10, "10_fit_completion_panel", width = 10, height = 7)
  }
  invisible(NULL)
}

write_report_markdown <- function(method_summary,
                                  sgl_summary,
                                  split_summary,
                                  component_strata,
                                  paired,
                                  manifest_status) {
  report_path <- r2r2_file("REVIEWER2_ROUND2_SIMULATION_REPORT.md")
  completed <- sum(manifest_status$valid_complete, na.rm = TRUE)
  attempted <- nrow(manifest_status)
  failure_rows <- if (nrow(method_summary)) {
    method_summary |>
      dplyr::filter(.data$failure_rate > 0) |>
      dplyr::arrange(dplyr::desc(.data$failure_rate))
  } else {
    data.frame()
  }

  lines <- c(
    "# Reviewer #2 Round 2 Additional Simulation Report",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## 1. Existing Simulation Audit",
    "",
    "- Existing sample-size battery: `results/final_res/compare_mixture_methods_by_samplesize_30seeds.Rmd`; binary logistic outcome, 13 exposures, 3 groups, 30 seeds, q=4, n_boot=200.",
    "- Existing validation split-stability code: `analysis/validation_calibration/R/02_dgp.R` and `analysis/validation_split_stability_standalone.Rmd`; SGL-WQS-only fixed-data split design.",
    "- Existing high-dimensional checkpointed code: `analysis/compare_mixture_methods_highdim_p_checkpointed.Rmd`; repository audit found this implementation uses a Gaussian outcome branch. The new runner stores family explicitly and does not overwrite existing production results.",
    "- Local SGL-WQS source is loaded from `source/sglwqs_old_analysis/sglwqs-main` by default through `pkgload::load_all()`.",
    "",
    "## 2. Exact New DGP Definitions",
    "",
    "- Matched Gaussian sample-size DGP reuses the binary 13-exposure X, group structure, exposure coefficients, covariates, and train/validation settings. The outcome is `Y = eta_exposure + eta_covariate + epsilon`, with `epsilon ~ N(0, 1)`.",
    "- Weak-signal Gaussian scenarios multiply exposure coefficients by 0.5 and leave covariate effects unchanged.",
    "- High-dimensional and factorial scenarios use a latent-factor exposure generator with within-group target correlation 0.45 and cross-group target correlation 0.05 unless otherwise stated.",
    "- Binomial matched branches use `logit(P(Y=1)) = -3 + eta_exposure + eta_covariate`.",
    "",
    "## 3. Differences From Existing Binary Designs",
    "",
    "- Outcome family is the intended contrast. Method settings are kept common across methods: q=4, 60/40 train/validation split where relevant, n_boot=200 for production, 10-fold CV, `lambda.min`, minor-direction threshold 0.10.",
    "- qgcomp attribution is reported as coefficient-derived attribution and not interpreted as WQS-type constrained weights.",
    "- Validation-stage p-values are conditional/exploratory summaries and are not formal selective inference.",
    "",
    "## 4. Production Run Manifest",
    "",
    paste0("- Atomic run jobs requested: ", attempted),
    paste0("- Valid completed atomic jobs: ", completed),
    paste0("- Remaining/incomplete atomic jobs: ", attempted - completed),
    "",
    "## 5. Completion/Failure Summary",
    ""
  )
  if (nrow(method_summary)) {
    lines <- c(lines, capture.output(knitr::kable(
      method_summary |>
        dplyr::select(battery, family, n, p, signal_profile, effect_profile, group_structure, method, attempted, completed, fit_completion_rate, runtime_median_iqr) |>
        utils::head(80),
      format = "markdown", digits = 3
    )))
  } else {
    lines <- c(lines, "No method-level summary is available yet.")
  }
  lines <- c(lines, "", "## 6. Matched Gaussian Sample-Size Results", "")
  sample_tab <- method_summary |> dplyr::filter(grepl("sample_size", .data$battery))
  lines <- c(lines, if (nrow(sample_tab)) capture.output(knitr::kable(sample_tab, format = "markdown", digits = 3)) else "No completed sample-size results yet.")
  lines <- c(lines, "", "## 7. Matched Gaussian High-Dimensional Results", "")
  high_tab <- method_summary |> dplyr::filter(grepl("highdim", .data$battery))
  lines <- c(lines, if (nrow(high_tab)) capture.output(knitr::kable(high_tab, format = "markdown", digits = 3)) else "No completed high-dimensional results yet.")
  lines <- c(lines, "", "## 8. Binary-Gaussian Split-Stability Comparison", "")
  lines <- c(lines, if (nrow(split_summary)) capture.output(knitr::kable(split_summary, format = "markdown", digits = 3)) else "No completed split-stability results yet.")
  lines <- c(lines, "", "## 9. Unbalanced Group x Effect Magnitude Results", "")
  dtab <- method_summary |> dplyr::filter(.data$battery == "D_unbalanced_effect_factorial")
  lines <- c(lines, if (nrow(dtab)) capture.output(knitr::kable(dtab, format = "markdown", digits = 3)) else "No completed factorial results yet.")
  lines <- c(lines, "", "## 10. Weak-Signal Hard-Setting Results", "")
  etab <- method_summary |> dplyr::filter(.data$battery == "E_hard_setting_sample_size")
  lines <- c(lines, if (nrow(etab)) capture.output(knitr::kable(etab, format = "markdown", digits = 3)) else "No completed hard-setting results yet.")
  lines <- c(lines, "", "## 11. What The Results Support", "")
  lines <- c(lines, "- Use completed scenario-level summaries to describe computational completion, group-direction assignment, component-level prioritization, and validation-stage conditional detection as distinct layers.")
  lines <- c(lines, "- Treat performance as scenario-dependent, especially under weak signal, small n, and medium/high-dimensional panels.")
  lines <- c(lines, "", "## 12. What The Results Do Not Support", "")
  lines <- c(lines, "- These simulations do not establish universal reliability, formal post-selection inference, causal signal detection, or superiority over all comparators.")
  lines <- c(lines, "- qgcomp coefficient-derived attribution should not be equated with WQS-type constrained weights.")
  lines <- c(lines, "", "## 13. Factual Points Useful For The Response Letter", "")
  lines <- c(lines, "- The added Gaussian batteries directly match the binary DGP structure where applicable and isolate the outcome-family contrast.")
  lines <- c(lines, "- Failures, bootstrap success, runtime, excluded group-direction terms, and conditional validation-stage p-value denominators are preserved.")
  lines <- c(lines, "", "## 14. Candidate Tiered Usage Table", "")
  evidence <- make_tiered_usage_evidence(method_summary, sgl_summary)
  lines <- c(lines, if (nrow(evidence)) capture.output(knitr::kable(evidence, format = "markdown", digits = 3)) else "No tiered evidence table available yet.")
  lines <- c(lines, "", "## 15. Remaining Limitations", "")
  lines <- c(lines, "- Exact wording for manuscript and response letter should be drafted only after all production jobs are complete and QA tables have been reviewed.")
  if (nrow(failure_rows)) {
    lines <- c(lines, "", "## Failure Rows Requiring Review", capture.output(knitr::kable(failure_rows, format = "markdown", digits = 3)))
  }
  writeLines(lines, report_path)
  invisible(report_path)
}

drop_smoke_rows <- function(df) {
  if (!is.data.frame(df) || !nrow(df) || !"scenario_id" %in% names(df)) return(df)
  df[!grepl("^smoke_", df$scenario_id) & !grepl("^smoke_failure", df$battery %||% ""), , drop = FALSE]
}

aggregate_and_write_outputs <- function(manifest = NULL, include_smoke = FALSE) {
  r2r2_make_dirs()
  results <- read_all_job_results()
  method_metrics <- bind_result_slot(results, "method_metrics")
  component_metrics <- bind_result_slot(results, "component_metrics")
  sglwqs_dir <- bind_result_slot(results, "sglwqs_direction_results")
  data_diag <- bind_result_slot(results, "data_diagnostics")
  if (!include_smoke) {
    method_metrics <- drop_smoke_rows(method_metrics)
    component_metrics <- drop_smoke_rows(component_metrics)
    sglwqs_dir <- drop_smoke_rows(sglwqs_dir)
    data_diag <- drop_smoke_rows(data_diag)
  }
  if (!is.null(manifest)) {
    manifest_status <- summarize_manifest_status(manifest)
  } else {
    manifest_status <- data.frame()
  }

  method_summary <- summarize_method_performance(method_metrics)
  sgl_summary <- summarize_sglwqs_scenario_validation(sglwqs_dir)
  sgl_by_direction <- summarize_sglwqs_validation(sglwqs_dir)
  component_strata <- summarize_component_strata(component_metrics)
  split_summary <- summarize_split_stability(sglwqs_dir)
  paired <- paired_method_comparisons(method_metrics)
  evidence <- make_tiered_usage_evidence(method_summary, sgl_summary)

  atomic_write_csv(method_metrics, r2r2_result_file("tables", "atomic_method_metrics.csv"))
  atomic_write_csv(component_metrics, r2r2_result_file("tables", "atomic_component_metrics.csv"))
  atomic_write_csv(sglwqs_dir, r2r2_result_file("tables", "atomic_sglwqs_group_direction_results.csv"))
  atomic_write_csv(data_diag, r2r2_result_file("tables", "atomic_data_diagnostics.csv"))
  atomic_write_csv(method_summary, r2r2_result_file("tables", "01_method_performance_summary.csv"))
  atomic_write_csv(method_summary |> dplyr::filter(grepl("sample_size", .data$battery)), r2r2_result_file("tables", "01_gaussian_sample_size_performance.csv"))
  atomic_write_csv(method_summary |> dplyr::filter(grepl("highdim", .data$battery)), r2r2_result_file("tables", "03_gaussian_highdimensional_performance.csv"))
  atomic_write_csv(split_summary, r2r2_result_file("tables", "05_binary_vs_gaussian_split_stability.csv"))
  atomic_write_csv(method_summary |> dplyr::filter(.data$battery == "D_unbalanced_effect_factorial"), r2r2_result_file("tables", "06_balanced_vs_unbalanced_factorial_summary.csv"))
  atomic_write_csv(component_strata, r2r2_result_file("tables", "07_group_size_effect_tier_stratified_summary.csv"))
  atomic_write_csv(method_summary |> dplyr::filter(.data$battery == "E_hard_setting_sample_size"), r2r2_result_file("tables", "08_hard_setting_sample_size_summary.csv"))
  atomic_write_csv(method_summary |> dplyr::select(battery, scenario_id, family, n, p, method, fit_completion_rate, runtime_median_iqr), r2r2_result_file("tables", "09_fit_completion_runtime.csv"))
  atomic_write_csv(sgl_summary, r2r2_result_file("tables", "10_sglwqs_conditional_retention_rejection_summary.csv"))
  atomic_write_csv(sgl_by_direction, r2r2_result_file("tables", "10_sglwqs_group_direction_detail.csv"))
  atomic_write_csv(paired, r2r2_result_file("tables", "11_paired_method_differences.csv"))
  atomic_write_csv(evidence, r2r2_result_file("tables", "tiered_usage_evidence_table.csv"))

  write_markdown_table(method_summary, r2r2_result_file("tables", "01_method_performance_summary.md"), "Method performance summary")
  write_markdown_table(sgl_summary, r2r2_result_file("tables", "10_sglwqs_conditional_retention_rejection_summary.md"), "SGL-WQS conditional validation-stage summary")
  write_markdown_table(paired, r2r2_result_file("tables", "11_paired_method_differences.md"), "Paired method differences")
  write_markdown_table(evidence, r2r2_result_file("tables", "tiered_usage_evidence_table.md"), "Tiered usage evidence table")

  generate_core_figures(method_summary, sgl_summary, component_strata)
  write_report_markdown(method_summary, sgl_summary, split_summary, component_strata, paired, manifest_status)

  invisible(list(
    results = results,
    method_metrics = method_metrics,
    component_metrics = component_metrics,
    sglwqs_direction_results = sglwqs_dir,
    method_summary = method_summary,
    sgl_summary = sgl_summary,
    split_summary = split_summary,
    component_strata = component_strata,
    paired = paired,
    evidence = evidence,
    manifest_status = manifest_status
  ))
}
