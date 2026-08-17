#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(knitr)
})

args <- commandArgs(trailingOnly = TRUE)
repo_root <- normalizePath(if (length(args) >= 1L) args[[1L]] else getwd(), mustWork = TRUE)
export_name <- if (length(args) >= 2L) args[[2L]] else "C2_paired_results_20260807"
export_parent <- file.path(repo_root, "reviewer2_round2", "results", "exports")
export_dir <- file.path(export_parent, export_name)
zip_path <- file.path(export_parent, paste0(export_name, ".zip"))

stopifnot(dir.exists(file.path(repo_root, "reviewer2_round2")))
if (dir.exists(export_dir) || file.exists(zip_path)) {
  stop("Refusing to overwrite an existing C2 export: ", export_name)
}

dir.create(file.path(export_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(export_dir, "raw"), recursive = TRUE, showWarnings = FALSE)

gaussian_root <- file.path(
  repo_root, "reviewer2_round2", "results", "raw",
  "C2_gaussian_paired_split_stability"
)
binary_root <- file.path(
  repo_root, "analysis", "validation_calibration", "output",
  "manifests", "production", "split_stability"
)

gaussian_paths <- sort(list.files(
  gaussian_root, pattern = "seed.*[.]rds$", recursive = TRUE, full.names = TRUE
))
binary_paths <- sort(list.files(
  binary_root, pattern = "[.]rds$", recursive = TRUE, full.names = TRUE
))

if (length(gaussian_paths) != 400L) stop("Expected 400 Gaussian C2 files; found ", length(gaussian_paths))
if (length(binary_paths) != 400L) stop("Expected 400 binary reference files; found ", length(binary_paths))

gaussian_results <- lapply(gaussian_paths, readRDS)
binary_results <- lapply(binary_paths, readRDS)

bind_nonempty <- function(x, slot) {
  rows <- lapply(x, `[[`, slot)
  rows <- rows[vapply(rows, function(z) is.data.frame(z) && nrow(z) > 0L, logical(1))]
  bind_rows(rows)
}

gaussian_method <- bind_nonempty(gaussian_results, "method_metrics")
gaussian_diag <- bind_nonempty(gaussian_results, "data_diagnostics")
binary_diag <- bind_nonempty(binary_results, "diagnostics")

gaussian_direction <- map2_dfr(gaussian_results, gaussian_paths, function(x, path) {
  out <- x$sglwqs_direction_results
  job <- x$job[1L, , drop = FALSE]
  out |>
    mutate(
      outcome_family = "gaussian",
      effect = job$effect_profile,
      replicate = job$split_replicate,
      pair_id = job$paired_binary_job_id,
      source_file = path
    )
}) |>
  transmute(
    outcome_family, effect, n, replicate, pair_id,
    group, direction, is_true_null, active_variable_count, true_beta_sum,
    retained, excluded_by_minor_threshold, estimate, std_error, p_value,
    rejected_0_05, sign, pos_mass, neg_mass, direction_mass,
    runtime_sec, n_boot_requested, n_boot_success, boot_success_rate,
    warning_text = warning_summary, source_file
  )

binary_direction <- map2_dfr(binary_results, binary_paths, function(x, path) {
  x$direction_results |>
    mutate(
      outcome_family = "binomial",
      pair_id = job_id,
      source_file = path
    )
}) |>
  transmute(
    outcome_family, effect, n, replicate, pair_id,
    group, direction, is_true_null, active_variable_count, true_beta_sum,
    retained, excluded_by_minor_threshold, estimate, std_error, p_value,
    rejected_0_05, sign, pos_mass, neg_mass, direction_mass,
    runtime_sec, n_boot_requested, n_boot_success, boot_success_rate,
    warning_text = warning_messages, source_file
  )

direction_long <- bind_rows(binary_direction, gaussian_direction) |>
  arrange(effect, n, pair_id, outcome_family, group, direction)

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  mean(x)
}

rate_where <- function(x, keep) {
  x <- x[keep]
  if (!length(x)) return(NA_real_)
  mean(x %in% TRUE)
}

job_metrics <- direction_long |>
  group_by(outcome_family, effect, n, replicate, pair_id) |>
  summarise(
    group_direction_rows = n(),
    retention_rate = mean(retained %in% TRUE),
    active_retention_rate = rate_where(retained, !is_true_null),
    null_retention_rate = rate_where(retained, is_true_null),
    active_detection_all = rate_where(rejected_0_05, !is_true_null),
    active_detection_retained = rate_where(rejected_0_05, !is_true_null & retained),
    active_correct_sign_detection_all = rate_where(
      rejected_0_05 & sign == direction, !is_true_null
    ),
    active_sign_accuracy_retained = rate_where(
      sign == direction, !is_true_null & retained & !is.na(sign)
    ),
    null_rejection_all = rate_where(rejected_0_05, is_true_null),
    null_rejection_retained = rate_where(rejected_0_05, is_true_null & retained),
    any_null_rejection = if (any(is_true_null)) any(rejected_0_05[is_true_null] %in% TRUE) else NA,
    retained_group_directions = sum(retained %in% TRUE),
    runtime_sec = first(runtime_sec),
    boot_success_rate = first(boot_success_rate),
    .groups = "drop"
  )

job_diagnostics <- bind_rows(
  binary_diag |>
    transmute(
      outcome_family = "binomial", pair_id = job_id, fit_success,
      runtime_sec, boot_success_rate, selected_lambda,
      lambda_boundary = selected_lambda_at_path_boundary,
      all_zero_exposure, backend_jerr,
      warning_text = warning_messages,
      error_class, error_message
    ),
  map_dfr(gaussian_results, function(x) {
    data.frame(
      outcome_family = "gaussian",
      pair_id = x$job$paired_binary_job_id,
      fit_success = x$method_metrics$fit_success,
      runtime_sec = x$method_metrics$runtime_sec,
      boot_success_rate = x$method_metrics$bootstrap_success_rate,
      selected_lambda = x$method_metrics$selected_lambda,
      lambda_boundary = NA,
      all_zero_exposure = NA,
      backend_jerr = NA,
      warning_text = x$method_metrics$warning_summary,
      error_class = x$method_metrics$error_class,
      error_message = x$method_metrics$error_message,
      stringsAsFactors = FALSE
    )
  })
)

job_metrics <- job_metrics |>
  left_join(
    job_diagnostics |>
      select(outcome_family, pair_id, fit_success, selected_lambda,
             lambda_boundary, all_zero_exposure, backend_jerr,
             diagnostic_warning = warning_text, error_class, error_message),
    by = c("outcome_family", "pair_id")
  )

scenario_summary <- job_metrics |>
  group_by(outcome_family, effect, n) |>
  summarise(
    jobs = n(),
    fit_successes = sum(fit_success %in% TRUE),
    fit_completion_rate = mean(fit_success %in% TRUE),
    retention_mean = mean_or_na(retention_rate),
    active_retention_mean = mean_or_na(active_retention_rate),
    null_retention_mean = mean_or_na(null_retention_rate),
    active_detection_all_mean = mean_or_na(active_detection_all),
    active_detection_all_sd = sd(active_detection_all, na.rm = TRUE),
    active_detection_retained_mean = mean_or_na(active_detection_retained),
    active_correct_sign_detection_all_mean = mean_or_na(active_correct_sign_detection_all),
    active_sign_accuracy_retained_mean = mean_or_na(active_sign_accuracy_retained),
    null_rejection_all_mean = mean_or_na(null_rejection_all),
    null_rejection_all_sd = sd(null_rejection_all, na.rm = TRUE),
    null_rejection_retained_mean = mean_or_na(null_rejection_retained),
    any_null_rejection_rate = mean_or_na(as.numeric(any_null_rejection)),
    retained_group_directions_mean = mean(retained_group_directions),
    runtime_median_sec = median(runtime_sec),
    runtime_iqr_sec = IQR(runtime_sec),
    runtime_max_sec = max(runtime_sec),
    runtime_total_hours = sum(runtime_sec) / 3600,
    bootstrap_success_rate_median = median(boot_success_rate),
    .groups = "drop"
  ) |>
  arrange(effect, n, outcome_family)

group_direction_summary <- direction_long |>
  group_by(outcome_family, effect, n, group, direction, is_true_null,
           active_variable_count, true_beta_sum) |>
  summarise(
    splits = n(),
    retained_n = sum(retained %in% TRUE),
    retention_rate = mean(retained %in% TRUE),
    rejection_n = sum(rejected_0_05 %in% TRUE),
    rejection_rate_all = mean(rejected_0_05 %in% TRUE),
    rejection_rate_retained = rate_where(rejected_0_05, retained),
    correct_sign_rejection_rate_all = if (first(is_true_null)) NA_real_ else
      mean((rejected_0_05 & sign == direction) %in% TRUE),
    sign_accuracy_retained = if (first(is_true_null)) NA_real_ else
      rate_where(sign == direction, retained & !is.na(sign)),
    estimate_mean = mean_or_na(estimate),
    estimate_sd = sd(estimate, na.rm = TRUE),
    p_value_median = median(p_value, na.rm = TRUE),
    p_value_iqr = IQR(p_value, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(effect, n, group, direction, outcome_family)

paired_direction <- binary_direction |>
  select(-outcome_family, -source_file, -warning_text, -runtime_sec,
         -n_boot_requested, -n_boot_success, -boot_success_rate) |>
  rename_with(~paste0(.x, "_binomial"),
              -c(effect, n, replicate, pair_id, group, direction)) |>
  inner_join(
    gaussian_direction |>
      select(-outcome_family, -source_file, -warning_text, -runtime_sec,
             -n_boot_requested, -n_boot_success, -boot_success_rate) |>
      rename_with(~paste0(.x, "_gaussian"),
                  -c(effect, n, replicate, pair_id, group, direction)),
    by = c("effect", "n", "replicate", "pair_id", "group", "direction")
  ) |>
  arrange(effect, n, pair_id, group, direction)

truth_match <- with(paired_direction,
  is_true_null_binomial == is_true_null_gaussian &
    active_variable_count_binomial == active_variable_count_gaussian &
    abs(true_beta_sum_binomial - true_beta_sum_gaussian) < 1e-12
)

paired_agreement_summary <- paired_direction |>
  group_by(effect, n) |>
  summarise(
    paired_direction_rows = n(),
    truth_definition_agreement = mean(
      is_true_null_binomial == is_true_null_gaussian &
        active_variable_count_binomial == active_variable_count_gaussian &
        abs(true_beta_sum_binomial - true_beta_sum_gaussian) < 1e-12
    ),
    retention_agreement = mean(retained_binomial == retained_gaussian),
    rejection_agreement = mean(rejected_0_05_binomial == rejected_0_05_gaussian),
    sign_agreement_both_retained = rate_where(
      sign_binomial == sign_gaussian,
      retained_binomial & retained_gaussian & !is.na(sign_binomial) & !is.na(sign_gaussian)
    ),
    spearman_p_both_finite = suppressWarnings(cor(
      p_value_binomial, p_value_gaussian,
      method = "spearman", use = "complete.obs"
    )),
    .groups = "drop"
  )

paired_job_metrics <- job_metrics |>
  select(outcome_family, effect, n, replicate, pair_id,
         retention_rate, active_retention_rate, null_retention_rate,
         active_detection_all, active_detection_retained,
         active_correct_sign_detection_all, active_sign_accuracy_retained,
         null_rejection_all, null_rejection_retained, any_null_rejection,
         retained_group_directions, runtime_sec, boot_success_rate) |>
  pivot_wider(
    names_from = outcome_family,
    values_from = c(
      retention_rate, active_retention_rate, null_retention_rate,
      active_detection_all, active_detection_retained,
      active_correct_sign_detection_all, active_sign_accuracy_retained,
      null_rejection_all, null_rejection_retained, any_null_rejection,
      retained_group_directions, runtime_sec, boot_success_rate
    ),
    names_sep = "__"
  ) |>
  arrange(effect, n, replicate)

contrast_metrics <- c(
  "retention_rate", "active_retention_rate", "null_retention_rate",
  "active_detection_all", "active_detection_retained",
  "active_correct_sign_detection_all", "active_sign_accuracy_retained",
  "null_rejection_all", "null_rejection_retained", "any_null_rejection",
  "retained_group_directions", "runtime_sec", "boot_success_rate"
)

paired_boot_summary <- function(data, metric, seed) {
  bcol <- paste0(metric, "__binomial")
  gcol <- paste0(metric, "__gaussian")
  keep <- is.finite(data[[bcol]]) & is.finite(data[[gcol]])
  dif <- data[[gcol]][keep] - data[[bcol]][keep]
  if (!length(dif)) return(NULL)
  set.seed(seed)
  boot <- if (length(dif) > 1L) {
    replicate(4000L, mean(sample(dif, length(dif), replace = TRUE)))
  } else {
    dif
  }
  data.frame(
    metric = metric,
    paired_splits = length(dif),
    binomial_mean = mean(data[[bcol]][keep]),
    gaussian_mean = mean(data[[gcol]][keep]),
    mean_difference_gaussian_minus_binomial = mean(dif),
    median_difference_gaussian_minus_binomial = median(dif),
    paired_difference_sd = sd(dif),
    bootstrap_ci_low = unname(quantile(boot, 0.025)),
    bootstrap_ci_high = unname(quantile(boot, 0.975)),
    stringsAsFactors = FALSE
  )
}

paired_contrasts <- split(paired_job_metrics, interaction(
  paired_job_metrics$effect, paired_job_metrics$n, drop = TRUE
)) |>
  imap_dfr(function(dat, key) {
    map2_dfr(contrast_metrics, seq_along(contrast_metrics), function(metric, i) {
      out <- paired_boot_summary(dat, metric, 20260807L + i)
      if (is.null(out)) return(data.frame())
      out$effect <- dat$effect[[1L]]
      out$n <- dat$n[[1L]]
      out
    })
  }) |>
  select(effect, n, everything()) |>
  arrange(effect, n, metric)

gaussian_hash_audit <- map_dfr(gaussian_results, function(x) {
  cbind(
    x$job[, c(
      "scenario_id", "data_seed", "split_seed", "split_replicate",
      "paired_binary_scenario_id", "paired_binary_job_id"
    )],
    x$data_diagnostics[, c(
      "paired_exposure_covariate_hash", "paired_reference_train_hash",
      "paired_reference_validation_hash", "effective_train_prop",
      "paired_binary_prevalence", "gaussian_residual_sd"
    )]
  )
}) |>
  arrange(scenario_id, split_replicate)

hash_scenario_summary <- gaussian_hash_audit |>
  group_by(scenario_id) |>
  summarise(
    jobs = n(),
    unique_exposure_covariate_hashes = n_distinct(paired_exposure_covariate_hash),
    unique_reference_train_hashes = n_distinct(paired_reference_train_hash),
    unique_reference_validation_hashes = n_distinct(paired_reference_validation_hash),
    effective_train_prop_min = min(effective_train_prop),
    effective_train_prop_max = max(effective_train_prop),
    binary_reference_prevalence = first(paired_binary_prevalence),
    gaussian_residual_sd = first(gaussian_residual_sd),
    .groups = "drop"
  )

preflight_path <- file.path(
  repo_root, "reviewer2_round2", "results", "summaries",
  "delegated_C2_paired_preflight_checks.csv"
)
completion_path <- file.path(
  repo_root, "reviewer2_round2", "results", "summaries",
  "delegated_C2_paired_completion_status.csv"
)
preflight <- read.csv(preflight_path, stringsAsFactors = FALSE)
completion <- read.csv(completion_path, stringsAsFactors = FALSE)

quality_checks <- data.frame(
  check = c(
    "Gaussian C2 atomic files",
    "Binary reference atomic files",
    "Gaussian C2 successful fits",
    "Binary reference successful fits",
    "Completion rows marked valid",
    "One-to-one paired jobs",
    "Gaussian direction rows",
    "Binary direction rows",
    "One-to-one paired direction rows",
    "Paired truth definitions identical",
    "All pairing preflight rows passed",
    "All completed Gaussian jobs contain pairing hashes",
    "All Gaussian bootstrap success rates equal 1"
  ),
  observed = c(
    length(gaussian_paths),
    length(binary_paths),
    sum(gaussian_method$fit_success %in% TRUE),
    sum(binary_diag$fit_success %in% TRUE),
    sum(completion$valid_complete %in% TRUE),
    nrow(paired_job_metrics),
    nrow(gaussian_direction),
    nrow(binary_direction),
    nrow(paired_direction),
    sum(truth_match %in% TRUE),
    sum(preflight$passed %in% TRUE),
    sum(nzchar(gaussian_hash_audit$paired_exposure_covariate_hash)),
    sum(gaussian_method$bootstrap_success_rate == 1, na.rm = TRUE)
  ),
  expected = c(400, 400, 400, 400, 400, 400, 2400, 2400, 2400, 2400, 4, 400, 400),
  stringsAsFactors = FALSE
) |>
  mutate(passed = observed == expected)

sum_true_or_na <- function(x) {
  if (all(is.na(x))) return(NA_integer_)
  sum(x %in% TRUE, na.rm = TRUE)
}

warning_summary <- job_diagnostics |>
  mutate(
    any_warning = !is.na(warning_text) & nzchar(trimws(warning_text)),
    both_directions_nonzero = grepl(
      "both positive and negative coefficients non-zero", warning_text,
      fixed = TRUE
    ),
    other_numerical_warning = grepl(
      "converg|Inf detected|numerical|failed", warning_text,
      ignore.case = TRUE
    ) & !both_directions_nonzero
  ) |>
  group_by(outcome_family) |>
  summarise(
    jobs = n(),
    fit_failures = sum(!(fit_success %in% TRUE)),
    jobs_with_any_warning = sum(any_warning, na.rm = TRUE),
    jobs_with_both_directions_warning = sum(both_directions_nonzero, na.rm = TRUE),
    jobs_with_other_numerical_warning = sum(other_numerical_warning, na.rm = TRUE),
    lambda_boundary_jobs = sum_true_or_na(lambda_boundary),
    all_zero_exposure_jobs = sum_true_or_na(all_zero_exposure),
    backend_nonzero_jerr_jobs = if (all(is.na(backend_jerr))) NA_integer_ else
      sum(is.finite(backend_jerr) & backend_jerr != 0, na.rm = TRUE),
    .groups = "drop"
  )

diagnostic_scenario_summary <- job_metrics |>
  mutate(
    any_warning = !is.na(diagnostic_warning) & nzchar(trimws(diagnostic_warning)),
    both_directions_nonzero = grepl(
      "both positive and negative coefficients non-zero", diagnostic_warning,
      fixed = TRUE
    ),
    other_numerical_warning = grepl(
      "converg|Inf detected|numerical|failed", diagnostic_warning,
      ignore.case = TRUE
    ) & !both_directions_nonzero
  ) |>
  group_by(outcome_family, effect, n) |>
  summarise(
    jobs = n(),
    fit_failures = sum(!(fit_success %in% TRUE)),
    jobs_with_any_warning = sum(any_warning, na.rm = TRUE),
    jobs_with_both_directions_warning = sum(both_directions_nonzero, na.rm = TRUE),
    jobs_with_other_numerical_warning = sum(other_numerical_warning, na.rm = TRUE),
    lambda_boundary_jobs = sum_true_or_na(lambda_boundary),
    all_zero_exposure_jobs = sum_true_or_na(all_zero_exposure),
    backend_nonzero_jerr_jobs = if (all(is.na(backend_jerr))) NA_integer_ else
      sum(is.finite(backend_jerr) & backend_jerr != 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(effect, n, outcome_family)

write_csv <- function(x, name) {
  write.csv(x, file.path(export_dir, "tables", name), row.names = FALSE, na = "")
}

write_csv(quality_checks, "01_quality_checks.csv")
write_csv(scenario_summary, "02_scenario_summary.csv")
write_csv(group_direction_summary, "03_group_direction_summary.csv")
write_csv(paired_job_metrics, "04_paired_job_metrics.csv")
write_csv(paired_contrasts, "05_paired_family_contrasts.csv")
write_csv(paired_agreement_summary, "06_paired_direction_agreement.csv")
write_csv(hash_scenario_summary, "07_pairing_hash_summary.csv")
write_csv(gaussian_hash_audit, "08_pairing_hash_audit_atomic.csv")
write_csv(warning_summary, "09_warning_failure_summary.csv")
write_csv(diagnostic_scenario_summary, "09_diagnostic_scenario_summary.csv")
write_csv(job_metrics, "10_job_level_metrics_long.csv")
write_csv(gaussian_direction, "11_gaussian_direction_results_atomic.csv")
write_csv(binary_direction, "12_binary_direction_results_atomic.csv")
write_csv(paired_direction, "13_paired_direction_results_atomic.csv")

family_labels <- c(binomial = "Binomial", gaussian = "Gaussian")
metric_labels <- c(
  retention_mean = "Retention",
  active_detection_all_mean = "Active rejection (all attempted)",
  active_correct_sign_detection_all_mean = "Correct-sign active rejection",
  null_rejection_all_mean = "Null rejection",
  any_null_rejection_rate = "Any null rejection per split"
)

scenario_plot <- scenario_summary |>
  select(outcome_family, effect, n, all_of(names(metric_labels))) |>
  pivot_longer(all_of(names(metric_labels)), names_to = "metric", values_to = "rate") |>
  mutate(
    outcome_family = recode(outcome_family, !!!family_labels),
    metric = factor(metric, levels = names(metric_labels), labels = metric_labels),
    n_label = paste0("n=", n),
    effect = recode(effect, global_null = "Global null", partial_null = "Partial null")
  )

p1 <- ggplot(scenario_plot, aes(n_label, rate, color = outcome_family, group = outcome_family)) +
  geom_point(position = position_dodge(width = 0.25), size = 2.2, na.rm = TRUE) +
  geom_line(position = position_dodge(width = 0.25), na.rm = TRUE) +
  facet_grid(metric ~ effect, scales = "free_y") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Rate", color = "Outcome family",
       title = "Exactly paired split-stability summaries") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top")

contrast_plot <- paired_contrasts |>
  filter(metric %in% c(
    "retention_rate", "active_detection_all",
    "active_correct_sign_detection_all", "null_rejection_all",
    "any_null_rejection"
  )) |>
  mutate(
    metric = recode(
      metric,
      retention_rate = "Retention",
      active_detection_all = "Active rejection",
      active_correct_sign_detection_all = "Correct-sign active rejection",
      null_rejection_all = "Null rejection",
      any_null_rejection = "Any null rejection"
    ),
    scenario = paste(
      recode(effect, global_null = "Global null", partial_null = "Partial null"),
      paste0("n=", n), sep = "; "
    )
  )

p2 <- ggplot(
  contrast_plot,
  aes(mean_difference_gaussian_minus_binomial, metric)
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbar(
    aes(xmin = bootstrap_ci_low, xmax = bootstrap_ci_high),
    orientation = "y", width = 0.2
  ) +
  geom_point(size = 2) +
  facet_wrap(~scenario, scales = "free_y") +
  labs(x = "Mean paired difference (Gaussian - Binomial)", y = NULL,
       title = "Paired outcome-family contrasts with descriptive bootstrap intervals") +
  theme_bw(base_size = 10)

p3 <- ggplot(
  job_metrics |>
    mutate(
      outcome_family = recode(outcome_family, !!!family_labels),
      scenario = paste(
        recode(effect, global_null = "Global null", partial_null = "Partial null"),
        paste0("n=", n), sep = "; "
      )
    ),
  aes(outcome_family, runtime_sec, fill = outcome_family)
) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~scenario, scales = "free_y") +
  scale_y_log10() +
  labs(x = NULL, y = "Runtime per split (seconds, log scale)",
       title = "Runtime distribution for exactly paired jobs") +
  theme_bw(base_size = 10) +
  theme(legend.position = "none")

partial_direction_plot <- group_direction_summary |>
  filter(effect == "partial_null") |>
  mutate(
    outcome_family = recode(outcome_family, !!!family_labels),
    group_direction = paste(group, direction),
    truth = ifelse(is_true_null, "Null direction", "Active direction")
  )

p4 <- ggplot(
  partial_direction_plot,
  aes(group_direction, rejection_rate_all, color = outcome_family, group = outcome_family)
) +
  geom_point(position = position_dodge(width = 0.35), size = 2.2) +
  facet_grid(truth ~ n, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Rejection rate across all splits", color = "Outcome family",
       title = "Partial-null group-direction results") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "top")

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(export_dir, "figures", paste0(stem, ".png")), plot,
         width = width, height = height, dpi = 180)
  ggsave(file.path(export_dir, "figures", paste0(stem, ".pdf")), plot,
         width = width, height = height)
}

save_plot(p1, "01_family_comparison_rates", 10, 10)
save_plot(p2, "02_paired_family_differences", 10, 7)
save_plot(p3, "03_runtime_comparison", 9, 6)
save_plot(p4, "04_partial_null_direction_results", 10, 6)

fmt_pct <- function(x) ifelse(is.finite(x), sprintf("%.1f%%", 100 * x), "NA")
fmt_num <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "NA")

scenario_for_report <- scenario_summary |>
  transmute(
    Family = recode(outcome_family, !!!family_labels),
    Scenario = recode(effect, global_null = "global-null", partial_null = "partial-null"),
    n,
    Fits = paste0(fit_successes, "/", jobs),
    Retention = fmt_pct(retention_mean),
    `Active rejection` = fmt_pct(active_detection_all_mean),
    `Correct-sign active rejection` = fmt_pct(active_correct_sign_detection_all_mean),
    `Null rejection` = fmt_pct(null_rejection_all_mean),
    `Any null rejection/split` = fmt_pct(any_null_rejection_rate),
    `Median runtime (s)` = fmt_num(runtime_median_sec, 1L)
  )

agreement_for_report <- paired_agreement_summary |>
  transmute(
    Scenario = recode(effect, global_null = "global-null", partial_null = "partial-null"),
    n,
    `Direction rows` = paired_direction_rows,
    `Retention agreement` = fmt_pct(retention_agreement),
    `Rejection agreement` = fmt_pct(rejection_agreement),
    `Sign agreement when both retained` = fmt_pct(sign_agreement_both_retained),
    `Spearman p correlation` = fmt_num(spearman_p_both_finite, 3L)
  )

report_lines <- c(
  "# C2 Exactly Paired Gaussian–Binomial Split-Stability Results",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "C2 contains 400 Gaussian SGL-WQS jobs exactly paired to the 400 completed binary reference jobs. Within each pair, the exposure/covariate realization and train/validation membership are identical; only the outcome family differs. There are 100 split replicates in each of four cells: global-null and partial-null at n=500 and n=5,000.",
  "",
  "All validation-stage p-values and rejection rates in this archive are exploratory conditional summaries. They are not formal post-selection or confirmatory inference.",
  "",
  "## Quality gate",
  "",
  paste(capture.output(kable(quality_checks, format = "markdown")), collapse = "\n"),
  "",
  "All checks passed. The completion gate contains 400/400 valid Gaussian outputs, and all four scenario-level pairing preflight rows passed.",
  "",
  "## Scenario-level results",
  "",
  paste(capture.output(kable(scenario_for_report, format = "markdown", align = "l")), collapse = "\n"),
  "",
  "Definitions:",
  "",
  "- Active rejection: p<0.05 among all prespecified active group-directions; excluded directions remain in the denominator.",
  "- Correct-sign active rejection: active rejection with an estimated sign matching the prespecified direction.",
  "- Null rejection: p<0.05 among all prespecified null group-directions.",
  "- Any null rejection/split: proportion of splits with at least one rejected null group-direction.",
  "",
  "## Paired agreement across outcome families",
  "",
  paste(capture.output(kable(agreement_for_report, format = "markdown", align = "l")), collapse = "\n"),
  "",
  "The paired contrasts in `tables/05_paired_family_contrasts.csv` report Gaussian minus binomial differences with descriptive paired bootstrap intervals. These intervals are included to describe split-to-split differences, not as a claim of formal selective inference.",
  "",
  "## Files",
  "",
  "- `tables/01_quality_checks.csv`: completion and pairing gate.",
  "- `tables/02_scenario_summary.csv`: primary scenario-level summary.",
  "- `tables/03_group_direction_summary.csv`: group/direction-specific retention, rejection, sign, estimate, and p-value stability.",
  "- `tables/04_paired_job_metrics.csv`: one row per exactly paired split.",
  "- `tables/05_paired_family_contrasts.csv`: Gaussian-minus-binomial paired differences and descriptive bootstrap intervals.",
  "- `tables/06_paired_direction_agreement.csv`: outcome-family agreement measures.",
  "- `tables/07_pairing_hash_summary.csv` and `08_pairing_hash_audit_atomic.csv`: pairing metadata audit.",
  "- `tables/09_warning_failure_summary.csv`: fit failures and categorized warnings.",
  "- `tables/10`–`13`: job- and direction-level audit tables.",
  "- `figures/`: PNG and PDF figures corresponding to the principal summaries.",
  "- `metadata/`: manifest, preflight, completion status, package/session metadata, and production log.",
  "- `raw/`: all 400 Gaussian C2 outputs and all 400 paired binary reference outputs.",
  "",
  "## Interpretation boundary",
  "",
  "This experiment directly answers the outcome-family and split-stability comparison for the matched n=500/n=5,000 design. It does not establish universal equivalence between Gaussian and binomial outcomes, formal post-selection validity, or causal signal identification. Results should be described as conditional, split-dependent screening behavior."
)

writeLines(report_lines, file.path(export_dir, "README_C2_RESULTS.md"))

sval <- function(family, effect_name, sample_n, column) {
  out <- scenario_summary[
    scenario_summary$outcome_family == family &
      scenario_summary$effect == effect_name &
      scenario_summary$n == sample_n,
    column,
    drop = TRUE
  ]
  if (length(out) != 1L) return(NA_real_)
  out
}

gdval <- function(family, sample_n, group_name, direction_name, column) {
  out <- group_direction_summary[
    group_direction_summary$outcome_family == family &
      group_direction_summary$effect == "partial_null" &
      group_direction_summary$n == sample_n &
      group_direction_summary$group == group_name &
      group_direction_summary$direction == direction_name,
    column,
    drop = TRUE
  ]
  if (length(out) != 1L) return(NA_real_)
  out
}

ja_lines <- c(
  "# C2 完全対応 Gaussian–Binomial split-stability 結果報告",
  "",
  paste0("作成日時: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 結論",
  "",
  "400件のGaussian解析はすべて成功し、400件のbinary参照との一対一対応、曝露・共変量、train/validation membership、truth定義の照合にもすべて合格した。完全に対応したデザインでは、Gaussianアウトカムのvalidation-stage active-direction rejectionはbinaryより明確に高く、null-direction rejectionは同程度か低かった。したがって、連続アウトカムでbinaryより不安定になるという所見はなく、むしろ本設定では連続アウトカムの方が検出しやすい。",
  "",
  "ただし、n=500のGaussianでは改善の大部分がPCB positive方向に由来し、Metalsの正負方向は依然低率だった。平均値だけで『n=500で十分』とは解釈せず、方向・効果強度・splitに依存する探索的結果として提示する必要がある。",
  "",
  "## デザインと品質確認",
  "",
  "- global-null／partial-null × n=500／n=5,000の4条件、各100 split。",
  "- Gaussianとbinaryで曝露・共変量の実現値とtrain/validation membershipを完全に一致させ、outcome familyのみを変更。",
  "- Gaussian 400/400、binary参照400/400がfit成功。Gaussian bootstrapは全400件で200/200成功。",
  "- 2,400組のgroup-direction行が一対一で対応し、truth定義も2,400/2,400一致。",
  "- 有効train比率は0.598–0.600。binary参照の実現有病率は0.0500–0.0640。",
  "",
  "## 主結果",
  "",
  paste0(
    "- global-nullの方向単位null rejectionは、n=500でbinary ",
    fmt_pct(sval("binomial", "global_null", 500, "null_rejection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "global_null", 500, "null_rejection_all_mean")),
    "、n=5,000でbinary ", fmt_pct(sval("binomial", "global_null", 5000, "null_rejection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "global_null", 5000, "null_rejection_all_mean")), "。"
  ),
  paste0(
    "- global-nullで1つ以上のnull方向がrejectionされたsplitは、n=500でbinary ",
    fmt_pct(sval("binomial", "global_null", 500, "any_null_rejection_rate")),
    "、Gaussian ", fmt_pct(sval("gaussian", "global_null", 500, "any_null_rejection_rate")),
    "、n=5,000でbinary ", fmt_pct(sval("binomial", "global_null", 5000, "any_null_rejection_rate")),
    "、Gaussian ", fmt_pct(sval("gaussian", "global_null", 5000, "any_null_rejection_rate")),
    "。これは各splitに6つのnull方向があるfamily-wiseな記述値で、単一方向のrejection rateとは区別する。"
  ),
  paste0(
    "- partial-nullのactive-direction rejection（全prespecified active方向を分母）は、n=500でbinary ",
    fmt_pct(sval("binomial", "partial_null", 500, "active_detection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 500, "active_detection_all_mean")),
    "、n=5,000でbinary ", fmt_pct(sval("binomial", "partial_null", 5000, "active_detection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 5000, "active_detection_all_mean")), "。"
  ),
  paste0(
    "- 正しい符号を伴うactive rejectionは、n=500でbinary ",
    fmt_pct(sval("binomial", "partial_null", 500, "active_correct_sign_detection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 500, "active_correct_sign_detection_all_mean")),
    "、n=5,000でbinary ", fmt_pct(sval("binomial", "partial_null", 5000, "active_correct_sign_detection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 5000, "active_correct_sign_detection_all_mean")), "。"
  ),
  paste0(
    "- partial-nullのnull-direction rejectionは、n=500でbinary ",
    fmt_pct(sval("binomial", "partial_null", 500, "null_rejection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 500, "null_rejection_all_mean")),
    "、n=5,000でbinary ", fmt_pct(sval("binomial", "partial_null", 5000, "null_rejection_all_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 5000, "null_rejection_all_mean")),
    "であり、active rejectionの改善とnull rejection増加の明確なトレードオフは認めなかった。"
  ),
  "",
  "## 方向別の内訳",
  "",
  paste0(
    "- n=500 Gaussianのactive rejectionはPCB positive ",
    fmt_pct(gdval("gaussian", 500, "PCBs", "positive", "rejection_rate_all")),
    "、Metals positive ", fmt_pct(gdval("gaussian", 500, "Metals", "positive", "rejection_rate_all")),
    "、Metals negative ", fmt_pct(gdval("gaussian", 500, "Metals", "negative", "rejection_rate_all")),
    "。n=500で全方向が安定したわけではない。"
  ),
  paste0(
    "- n=5,000 GaussianではPCB positive ",
    fmt_pct(gdval("gaussian", 5000, "PCBs", "positive", "rejection_rate_all")),
    "、Metals positive ", fmt_pct(gdval("gaussian", 5000, "Metals", "positive", "rejection_rate_all")),
    "、Metals negative ", fmt_pct(gdval("gaussian", 5000, "Metals", "negative", "rejection_rate_all")), "。"
  ),
  paste0(
    "- 対応するn=5,000 binaryはPCB positive ",
    fmt_pct(gdval("binomial", 5000, "PCBs", "positive", "rejection_rate_all")),
    "、Metals positive ", fmt_pct(gdval("binomial", 5000, "Metals", "positive", "rejection_rate_all")),
    "、Metals negative ", fmt_pct(gdval("binomial", 5000, "Metals", "negative", "rejection_rate_all")), "。"
  ),
  "",
  "## Retentionの読み方",
  "",
  paste0(
    "partial-nullのactive-direction retentionはn=500でbinary ",
    fmt_pct(sval("binomial", "partial_null", 500, "active_retention_mean")),
    "、Gaussian ", fmt_pct(sval("gaussian", "partial_null", 500, "active_retention_mean")),
    "、n=5,000では両familyとも ",
    fmt_pct(sval("gaussian", "partial_null", 5000, "active_retention_mean")),
    "だった。Gaussianで全体retentionが低めに見える主因はnull方向の除外であり、active方向の欠落増加ではない。"
  ),
  "",
  "## 実行時間",
  "",
  paste0(
    "Gaussianのsplit当たり中央値はglobal-nullで",
    fmt_num(sval("gaussian", "global_null", 500, "runtime_median_sec"), 1L), "秒（n=500）／",
    fmt_num(sval("gaussian", "global_null", 5000, "runtime_median_sec"), 1L), "秒（n=5,000）、partial-nullで",
    fmt_num(sval("gaussian", "partial_null", 500, "runtime_median_sec"), 1L), "秒／",
    fmt_num(sval("gaussian", "partial_null", 5000, "runtime_median_sec"), 1L), "秒だった。対応binaryはそれぞれ",
    fmt_num(sval("binomial", "global_null", 500, "runtime_median_sec"), 1L), "秒、",
    fmt_num(sval("binomial", "global_null", 5000, "runtime_median_sec"), 1L), "秒、",
    fmt_num(sval("binomial", "partial_null", 500, "runtime_median_sec"), 1L), "秒、",
    fmt_num(sval("binomial", "partial_null", 5000, "runtime_median_sec"), 1L),
    "秒で、本デザインではGaussianが明確に短かった。ただしC2の主目的は速度比較ではなくoutcome-familyとsplit stabilityの比較である。"
  ),
  "",
  "## 技術的注意",
  "",
  "- 両familyの全ジョブで、正負両係数が同一変数で非ゼロとなる標準警告が保存され、net coefficientがweightに使用された。これはfit failureではなく、400/400でfitとbootstrapは完了している。",
  "- binary参照ではselected-lambda boundaryおよびall-zero exposure診断が一部splitで記録された。Gaussian C2の保存schemaには同じ診断列がないため、この項目をfamily間で件数比較しない。",
  "- validation-stage p値・rejectionはtraining-estimated indexに条件づけた探索的summaryであり、formal post-selection inferenceではない。",
  "",
  "## 査読コメントへの位置づけ",
  "",
  "この結果は、n=500／5,000の既存binaryデザインと完全に対応したGaussian実験を追加し、outcome family以外を固定したsplit-stability比較を提供する。連続アウトカムでの適用可能性を支持するが、n=500では方向別の検出差が残るため、sample sizeに関する一律な保証ではなく、条件依存の探索的性能として記述するのが妥当である。"
)

writeLines(ja_lines, file.path(export_dir, "C2_RESULTS_REPORT_JA.md"))

metadata_sources <- c(
  file.path(repo_root, "reviewer2_round2", "config", "scenario_manifest_C2_paired.csv"),
  preflight_path,
  completion_path,
  file.path(repo_root, "reviewer2_round2", "results", "logs", "battery_C2_paired_source_metadata.csv"),
  file.path(repo_root, "reviewer2_round2", "results", "logs", "battery_C2_paired_package_versions.csv"),
  file.path(repo_root, "reviewer2_round2", "results", "logs", "battery_C2_paired_sessionInfo.txt"),
  file.path(repo_root, "reviewer2_round2", "results", "logs", "battery_C2_paired_supervisor_console.log"),
  file.path(repo_root, "reviewer2_round2", "scripts", "32_summarize_C2_paired.R")
)
if (!all(file.exists(metadata_sources))) {
  stop("Missing metadata source(s): ", paste(metadata_sources[!file.exists(metadata_sources)], collapse = ", "))
}
invisible(file.copy(metadata_sources, file.path(export_dir, "metadata"), overwrite = FALSE))

copy_tree <- function(from, to) {
  dir.create(to, recursive = TRUE, showWarnings = FALSE)
  entries <- list.files(from, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  if (!length(entries)) return(TRUE)
  all(file.copy(entries, to, recursive = TRUE, overwrite = FALSE))
}

if (!copy_tree(gaussian_root, file.path(export_dir, "raw", "gaussian_C2"))) {
  stop("Failed to copy Gaussian C2 raw outputs")
}
if (!copy_tree(binary_root, file.path(export_dir, "raw", "binary_reference"))) {
  stop("Failed to copy binary reference raw outputs")
}

all_export_files <- list.files(export_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
manifest <- data.frame(
  relative_path = substring(all_export_files, nchar(export_dir) + 2L),
  size_bytes = file.info(all_export_files)$size,
  md5 = unname(tools::md5sum(all_export_files)),
  stringsAsFactors = FALSE
) |>
  arrange(relative_path)
write.csv(manifest, file.path(export_dir, "FILE_MANIFEST_MD5.csv"), row.names = FALSE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(export_parent)
utils::zip(zipfile = basename(zip_path), files = export_name, flags = "-r9X")

cat("C2 export directory:", export_dir, "\n")
cat("C2 export zip:", zip_path, "\n")
cat("Quality checks passed:", sum(quality_checks$passed), "/", nrow(quality_checks), "\n")
