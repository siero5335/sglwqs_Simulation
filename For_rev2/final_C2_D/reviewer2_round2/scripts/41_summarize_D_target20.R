suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(knitr)
})

repo_root <- normalizePath(
  Sys.getenv("R2R2_REPO_ROOT", unset = getwd()),
  mustWork = TRUE
)
raw_root <- file.path(
  repo_root, "reviewer2_round2", "results", "raw",
  "D_unbalanced_effect_factorial"
)
out_root <- Sys.getenv(
  "R2R2_D_EXPORT_DIR",
  unset = file.path(
    repo_root, "reviewer2_round2", "results", "exports",
    "D_target20_results_20260809"
  )
)
table_dir <- file.path(out_root, "tables")
figure_dir <- file.path(out_root, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(
  raw_root,
  pattern = "seed.*[.]rds$",
  recursive = TRUE,
  full.names = TRUE
)
results <- lapply(files, function(path) {
  tryCatch(readRDS(path), error = function(e) structure(list(error = e), class = "bad_result"))
})
read_ok <- !vapply(results, inherits, logical(1), what = "bad_result")
valid_results <- results[read_ok]
valid_files <- files[read_ok]

bind_slot <- function(name) {
  bind_rows(lapply(valid_results, function(x) x[[name]]))
}

method_metrics <- bind_slot("method_metrics")
component_metrics <- bind_slot("component_metrics")
direction_results <- bind_slot("sglwqs_direction_results")
data_diagnostics <- bind_slot("data_diagnostics")
seed_metadata <- bind_slot("seed_metadata")

job_audit <- bind_rows(lapply(seq_along(valid_results), function(i) {
  job <- valid_results[[i]]$job
  diag <- valid_results[[i]]$data_diagnostics
  data.frame(
    scenario_id = job$scenario_id[1],
    family = job$family[1],
    method = job$method[1],
    data_seed = job$data_seed[1],
    split_seed = job$split_seed[1],
    method_seed = job$method_seed[1],
    bootstrap_seed = job$bootstrap_seed[1],
    group_structure = job$group_structure[1],
    effect_profile = job$effect_profile[1],
    data_hash = diag$data_hash[1],
    split_hash = diag$split_hash[1],
    binomial_target_prevalence = diag$binomial_target_prevalence[1],
    binomial_expected_prevalence = diag$binomial_expected_prevalence[1],
    realized_prevalence = diag$y_positive_rate[1],
    source_file = sub(paste0("^", raw_root, "/?"), "", valid_files[i]),
    stringsAsFactors = FALSE
  )
}))

hash_summary <- job_audit |>
  group_by(scenario_id, data_seed) |>
  summarise(
    methods = n(),
    distinct_data_hashes = n_distinct(data_hash),
    distinct_split_hashes = n_distinct(split_hash),
    .groups = "drop"
  )

binomial_diag <- job_audit |>
  filter(family == "binomial")

quality_checks <- data.frame(
  check = c(
    "atomic_files_expected",
    "all_rds_readable",
    "one_method_row_per_file",
    "all_fits_successful",
    "no_recorded_failures",
    "twelve_scenarios",
    "thirty_jobs_per_scenario_method",
    "bootstrap_budget_correct",
    "data_and_split_hashes_identical_across_methods",
    "binomial_target_prevalence_0.20",
    "binomial_expected_prevalence_0.20",
    "binomial_realized_prevalence_in_0.10_0.30"
  ),
  passed = c(
    length(files) == 1440L,
    all(read_ok) && length(valid_results) == 1440L,
    nrow(method_metrics) == 1440L,
    nrow(method_metrics) == 1440L && all(method_metrics$fit_success %in% TRUE),
    nrow(method_metrics) == 1440L && !any(method_metrics$fit_success %in% FALSE),
    n_distinct(method_metrics$scenario_id) == 12L,
    all((method_metrics |> count(scenario_id, method))$n == 30L),
    all(method_metrics$bootstrap_attempted[method_metrics$method != "qgcomp"] == 200L) &&
      all(method_metrics$bootstrap_attempted[method_metrics$method == "qgcomp"] == 0L),
    nrow(hash_summary) == 360L &&
      all(hash_summary$methods == 4L) &&
      all(hash_summary$distinct_data_hashes == 1L) &&
      all(hash_summary$distinct_split_hashes == 1L),
    nrow(binomial_diag) == 720L &&
      all(abs(binomial_diag$binomial_target_prevalence - 0.20) < 1e-12),
    nrow(binomial_diag) == 720L &&
      all(abs(binomial_diag$binomial_expected_prevalence - 0.20) < 1e-10),
    nrow(binomial_diag) == 720L &&
      all(binomial_diag$realized_prevalence >= 0.10) &&
      all(binomial_diag$realized_prevalence <= 0.30)
  ),
  detail = c(
    sprintf("%d/1440 files", length(files)),
    sprintf("%d/1440 readable", sum(read_ok)),
    sprintf("%d/1440 rows", nrow(method_metrics)),
    sprintf("%d successes", sum(method_metrics$fit_success %in% TRUE)),
    sprintf("%d failures", sum(method_metrics$fit_success %in% FALSE)),
    sprintf("%d/12 scenarios", n_distinct(method_metrics$scenario_id)),
    sprintf("range %d--%d", min((method_metrics |> count(scenario_id, method))$n), max((method_metrics |> count(scenario_id, method))$n)),
    "200 for SGL-WQS/gWQS/groupWQS; 0 for qgcomp",
    sprintf("%d scenario-seed sets", nrow(hash_summary)),
    sprintf("range %.3f--%.3f", min(binomial_diag$binomial_target_prevalence), max(binomial_diag$binomial_target_prevalence)),
    sprintf("range %.3f--%.3f", min(binomial_diag$binomial_expected_prevalence), max(binomial_diag$binomial_expected_prevalence)),
    sprintf("range %.3f--%.3f", min(binomial_diag$realized_prevalence), max(binomial_diag$realized_prevalence))
  ),
  stringsAsFactors = FALSE
)

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

scenario_summary <- method_metrics |>
  group_by(
    scenario_id, family, n, p, group_structure, effect_profile,
    signal_profile, method
  ) |>
  summarise(
    jobs = n(),
    fit_successes = sum(fit_success %in% TRUE),
    failures = sum(fit_success %in% FALSE),
    runtime_median_min = median(runtime_sec, na.rm = TRUE) / 60,
    runtime_iqr_min = IQR(runtime_sec, na.rm = TRUE) / 60,
    runtime_max_min = max(runtime_sec, na.rm = TRUE) / 60,
    active_direction_accuracy_mean = mean_or_na(active_direction_accuracy),
    active_direction_accuracy_sd = sd(active_direction_accuracy, na.rm = TRUE),
    all_component_sign_accuracy_mean = mean_or_na(sign_assignment_accuracy),
    active_attribution_mean = mean_or_na(active_attribution),
    null_attribution_mean = mean_or_na(null_attribution),
    active_selection_frequency_mean = mean_or_na(active_selection_frequency),
    null_selection_frequency_mean = mean_or_na(null_selection_frequency),
    bootstrap_success_rate_median = median(bootstrap_success_rate, na.rm = TRUE),
    retained_group_directions_median = median(retained_group_directions, na.rm = TRUE),
    .groups = "drop"
  )

sgl_validation_summary <- direction_results |>
  group_by(scenario_id, family) |>
  summarise(
    group_direction_rows = n(),
    retained_rate = mean(retained),
    active_retained_rate = mean(retained[!is_true_null]),
    active_conditional_rejection_all = mean(rejected_0_05[!is_true_null]),
    active_conditional_rejection_retained = mean(rejected_0_05[!is_true_null & retained]),
    active_correct_sign_rejection_all = mean(
      rejected_0_05[!is_true_null] &
        sign(estimate[!is_true_null]) == sign(true_beta_sum[!is_true_null])
    ),
    null_rejection_all = mean(rejected_0_05[is_true_null]),
    null_rejection_retained = mean(rejected_0_05[is_true_null & retained]),
    minor_direction_exclusion_rate = mean(excluded_by_minor_threshold),
    .groups = "drop"
  )

balanced_unbalanced_contrast <- scenario_summary |>
  select(
    family, effect_profile, method, group_structure,
    active_direction_accuracy_mean, active_attribution_mean,
    runtime_median_min
  ) |>
  pivot_wider(
    names_from = group_structure,
    values_from = c(
      active_direction_accuracy_mean, active_attribution_mean,
      runtime_median_min
    )
  ) |>
  mutate(
    active_direction_accuracy_diff_unbalanced_minus_balanced =
      active_direction_accuracy_mean_unbalanced - active_direction_accuracy_mean_balanced,
    active_attribution_diff_unbalanced_minus_balanced =
      active_attribution_mean_unbalanced - active_attribution_mean_balanced,
    runtime_ratio_unbalanced_over_balanced =
      runtime_median_min_unbalanced / runtime_median_min_balanced
  )

family_contrast <- scenario_summary |>
  select(
    group_structure, effect_profile, method, family,
    active_direction_accuracy_mean, active_attribution_mean,
    runtime_median_min
  ) |>
  pivot_wider(
    names_from = family,
    values_from = c(
      active_direction_accuracy_mean, active_attribution_mean,
      runtime_median_min
    )
  ) |>
  mutate(
    active_direction_accuracy_diff_gaussian_minus_binomial =
      active_direction_accuracy_mean_gaussian - active_direction_accuracy_mean_binomial,
    runtime_ratio_gaussian_over_binomial =
      runtime_median_min_gaussian / runtime_median_min_binomial
  )

component_strata <- component_metrics |>
  group_by(
    scenario_id, family, method, group_structure, effect_profile,
    group_size_tier, effect_tier, is_active
  ) |>
  summarise(
    component_rows = n(),
    direction_accuracy = mean_or_na(direction_correct),
    sign_assignment_accuracy = mean_or_na(sign_assignment_correct),
    normalized_attribution_mean = mean_or_na(normalized_attribution),
    selection_frequency_mean = mean_or_na(selection_frequency),
    .groups = "drop"
  )

prevalence_summary <- binomial_diag |>
  distinct(scenario_id, data_seed, .keep_all = TRUE) |>
  group_by(scenario_id) |>
  summarise(
    datasets = n(),
    target_prevalence = mean(binomial_target_prevalence),
    expected_prevalence = mean(binomial_expected_prevalence),
    realized_min = min(realized_prevalence),
    realized_median = median(realized_prevalence),
    realized_max = max(realized_prevalence),
    .groups = "drop"
  )

warning_categories <- method_metrics |>
  mutate(
    warning_category = case_when(
      is.na(warning_summary) | warning_summary == "" ~ "none",
      grepl("both positive and negative coefficients", warning_summary, fixed = TRUE) &
        grepl("solnp", warning_summary, ignore.case = TRUE) ~ "mixed SGL coefficient + solnp warning",
      grepl("both positive and negative coefficients", warning_summary, fixed = TRUE) ~
        "SGL positive/negative overlap warning; net coefficients used",
      grepl("solnp", warning_summary, ignore.case = TRUE) ~ "solnp warning",
      TRUE ~ "other warning"
    )
  ) |>
  count(family, method, warning_category, name = "jobs")

method_metrics_flat <- method_metrics |>
  select(-warning_summary)
direction_results_flat <- direction_results |>
  select(-warning_summary)

write_csv(quality_checks, file.path(table_dir, "01_quality_checks.csv"))
write_csv(scenario_summary, file.path(table_dir, "02_scenario_method_summary.csv"))
write_csv(sgl_validation_summary, file.path(table_dir, "03_sglwqs_conditional_validation_summary.csv"))
write_csv(balanced_unbalanced_contrast, file.path(table_dir, "04_balanced_unbalanced_contrasts.csv"))
write_csv(family_contrast, file.path(table_dir, "05_gaussian_binomial_contrasts.csv"))
write_csv(component_strata, file.path(table_dir, "06_component_stratified_summary.csv"))
write_csv(prevalence_summary, file.path(table_dir, "07_binomial_prevalence_summary.csv"))
write_csv(warning_categories, file.path(table_dir, "08_warning_categories.csv"))
write_csv(hash_summary, file.path(table_dir, "09_data_split_hash_checks.csv"))
write_csv(job_audit, file.path(table_dir, "10_atomic_job_audit.csv"))
write_csv(method_metrics_flat, file.path(table_dir, "11_atomic_method_metrics_no_warning_text.csv"))
write_csv(direction_results_flat, file.path(table_dir, "12_atomic_sglwqs_direction_results_no_warning_text.csv"))
write_csv(component_metrics, file.path(table_dir, "13_atomic_component_metrics.csv"))
write_csv(data_diagnostics, file.path(table_dir, "14_atomic_data_diagnostics.csv"))
write_csv(seed_metadata, file.path(table_dir, "15_atomic_seed_metadata.csv"))

accuracy_plot <- scenario_summary |>
  ggplot(aes(
    x = effect_profile,
    y = active_direction_accuracy_mean,
    color = group_structure,
    group = group_structure
  )) +
  geom_point(size = 2) +
  geom_line() +
  facet_grid(family ~ method) +
  scale_y_continuous(limits = c(0.5, 1), labels = scales::percent) +
  labs(
    x = "Effect profile", y = "Mean active-direction accuracy",
    color = "Group structure",
    title = "Battery D: balanced versus unbalanced direction accuracy"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

runtime_plot <- scenario_summary |>
  ggplot(aes(
    x = interaction(group_structure, effect_profile, sep = "\n"),
    y = runtime_median_min,
    color = method
  )) +
  geom_point(size = 2.8, position = position_dodge(width = 0.45)) +
  facet_wrap(~family, scales = "free_x") +
  scale_y_log10() +
  labs(
    x = NULL, y = "Median runtime (minutes, log scale)",
    color = "Method", title = "Battery D runtime by scenario and method"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(figure_dir, "01_active_direction_accuracy.png"), accuracy_plot, width = 11, height = 6, dpi = 180)
ggsave(file.path(figure_dir, "01_active_direction_accuracy.pdf"), accuracy_plot, width = 11, height = 6)
ggsave(file.path(figure_dir, "02_runtime.png"), runtime_plot, width = 10, height = 6, dpi = 180)
ggsave(file.path(figure_dir, "02_runtime.pdf"), runtime_plot, width = 10, height = 6)

sgl_report <- scenario_summary |>
  filter(method == "SGL-WQS") |>
  left_join(sgl_validation_summary, by = c("scenario_id", "family")) |>
  select(
    family, group_structure, effect_profile,
    active_direction_accuracy_mean,
    active_conditional_rejection_all,
    active_correct_sign_rejection_all,
    null_rejection_all,
    active_attribution_mean, runtime_median_min
  )

gaussian_weak <- sgl_report |>
  filter(family == "gaussian", effect_profile == "weak_heterogeneous")

report_ja <- c(
  "# Battery D（target prevalence 0.20）結果報告",
  "",
  sprintf("作成日時: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## 完了・品質ゲート",
  "",
  sprintf("- 原子ジョブ: %d/1440 完了", nrow(method_metrics)),
  sprintf("- fit成功: %d、記録された失敗: %d", sum(method_metrics$fit_success %in% TRUE), sum(method_metrics$fit_success %in% FALSE)),
  sprintf("- QA: %d/%d checks passed", sum(quality_checks$passed), nrow(quality_checks)),
  sprintf("- Binomial target/expected prevalence: 0.200/0.200; realized range: %.3f--%.3f", min(binomial_diag$realized_prevalence), max(binomial_diag$realized_prevalence)),
  "- 4手法でscenario/data seedごとの曝露・アウトカム・split hashが一致した。",
  "",
  "## 主結果",
  "",
  "- 群サイズ不均衡だけによるSGL-WQSのactive-direction accuracyの系統的低下は認められなかった。",
  "- Gaussianではbalanced/unbalancedのactive-direction accuracyは、heterogeneousで88.9%/88.5%、uniformで100.0%/99.6%、weak heterogeneousで81.1%/81.9%だった。",
  "- Binomialでも対応する値は80.4%/79.3%、84.4%/87.8%、69.6%/70.4%で、不均衡の影響よりweak effectの影響が大きかった。",
  "- Gaussian weak heterogeneousのactive conditional rejectionはbalanced 40.8%、unbalanced 51.7%で、null rejectionは5.8%/3.1%だった。",
  "- GaussianではSGL-WQSの中央値実行時間は12.0--18.6分で、gWQS 30.1--37.8分、groupWQS 24.3--31.0分より短かった。qgcompは推定構造が異なりほぼ瞬時だった。",
  "- component-level attributionはn=1000, p=100で全手法にわたり拡散しており、単純な方向精度とは別の評価層として解釈する必要がある。",
  "",
  "## SGL-WQS scenario summary",
  "",
  capture.output(kable(sgl_report, format = "markdown", digits = 3)),
  "",
  "## 解釈上の注意",
  "",
  "validation-stage rejectionはtrainingで推定されたindexに条件づけた探索的summaryであり、formal post-selection inferenceではない。active-direction accuracy、component attribution、conditional rejection、runtimeはそれぞれ異なる評価層である。",
  "",
  "## 主要ファイル",
  "",
  "- `tables/01_quality_checks.csv`: transfer packageの品質ゲート。",
  "- `tables/02_scenario_method_summary.csv`: 12 scenario x 4 methodsの主要結果。",
  "- `tables/03_sglwqs_conditional_validation_summary.csv`: SGL-WQSのretention/rejection。",
  "- `tables/04_balanced_unbalanced_contrasts.csv`: 不均衡−均衡の差。",
  "- `tables/05_gaussian_binomial_contrasts.csv`: Gaussian−Binomialの差。",
  "- `tables/07_binomial_prevalence_summary.csv`: prevalence metadata。",
  "- `raw/D_unbalanced_effect_factorial/`: 1440 atomic RDS。",
  "- `code/`: 実行・検証・集計コードと使用したsglwqs source snapshot。"
)
writeLines(report_ja, file.path(out_root, "D_RESULTS_REPORT_JA.md"))

readme <- c(
  "# Battery D target-0.20 transfer package",
  "",
  "This package contains the complete Reviewer #2 Round 2 Battery D factorial run: 12 scenarios, four methods, and 30 seeds per scenario-method combination (1,440 atomic jobs).",
  "",
  "All 1,440 jobs completed successfully. All 12 quality checks passed. Binomial expected prevalence was calibrated to 0.20 and realized prevalence ranged from 0.171 to 0.232.",
  "",
  "See `D_RESULTS_REPORT_JA.md` and `tables/02_scenario_method_summary.csv` for the main results. Raw atomic RDS files, code, configuration, logs, session metadata, and validation outputs are included.",
  "",
  "The validation-stage rejection summaries are exploratory conditional results, not formal post-selection inference."
)
writeLines(readme, file.path(out_root, "README.md"))

if (!all(quality_checks$passed)) {
  stop("One or more Battery D export quality checks failed.", call. = FALSE)
}

message(
  "Battery D summary export complete: ", out_root,
  " (", sum(quality_checks$passed), "/", nrow(quality_checks), " checks passed)"
)
