#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(knitr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

out_dir <- file.path("analysis", "results", "highdim_overall_summary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

method_levels <- c("qgcomp", "gWQS", "groupWQS", "SGL-WQS")

make_group_definitions <- function(p, n_groups = 10L) {
  vars <- sprintf("X%03d", seq_len(p))
  group_ids <- cut(seq_len(p),
                   breaks = n_groups,
                   labels = sprintf("G%02d", seq_len(n_groups)),
                   include.lowest = TRUE)
  split(vars, group_ids)
}

make_truth <- function(p, n_groups = 10L, truth_profile = "baseline3") {
  groups <- make_group_definitions(p, n_groups)
  vars <- unlist(groups, use.names = FALSE)
  beta <- setNames(rep(0, length(vars)), vars)
  profile <- match.arg(truth_profile, c("weak3", "baseline3", "strong3",
                                        "sparse3", "dense3", "many5"))
  signal_scale <- switch(profile, weak3 = 0.60, strong3 = 1.60, 1.00)

  add_group_effect <- function(group_index, direction, n_active, magnitude = c(0.08, 0.14)) {
    if (group_index > length(groups) || n_active <= 0L) return(invisible(NULL))
    group_vars <- groups[[group_index]]
    n_use <- min(length(group_vars), n_active)
    if (n_use <= 0L) return(invisible(NULL))
    effect <- signal_scale * seq(magnitude[1], magnitude[2], length.out = n_use)
    if (identical(direction, "negative")) effect <- -rev(effect)
    beta[group_vars[seq_len(n_use)]] <<- effect
    invisible(NULL)
  }

  add_mixed_group <- function(group_index, n_neg, n_pos) {
    if (group_index > length(groups)) return(invisible(NULL))
    group_vars <- groups[[group_index]]
    n_neg <- min(length(group_vars), n_neg)
    remaining <- max(length(group_vars) - n_neg, 0L)
    n_pos <- min(remaining, n_pos)
    if (n_neg > 0L) {
      beta[group_vars[seq_len(n_neg)]] <<- signal_scale * seq(-0.16, -0.11, length.out = n_neg)
    }
    if (n_pos > 0L) {
      pos_idx <- seq.int(n_neg + 1L, n_neg + n_pos)
      beta[group_vars[pos_idx]] <<- signal_scale * seq(0.10, 0.14, length.out = n_pos)
    }
    invisible(NULL)
  }

  if (profile %in% c("weak3", "baseline3", "strong3")) {
    g1 <- groups[[1]]
    g2 <- groups[[2]]
    g3 <- groups[[3]]
    add_group_effect(1, "positive", max(3L, ceiling(length(g1) * 0.15)), c(0.08, 0.12))
    add_mixed_group(2,
                    max(2L, floor(length(g2) * 0.08)),
                    max(2L, floor(length(g2) * 0.06)))
    add_group_effect(3, "negative", max(2L, ceiling(length(g3) * 0.08)), c(0.10, 0.14))
  } else if (profile == "sparse3") {
    add_group_effect(1, "positive", 2L, c(0.10, 0.13))
    add_mixed_group(2, 1L, 1L)
    add_group_effect(3, "negative", 2L, c(0.10, 0.13))
  } else if (profile == "dense3") {
    g1 <- groups[[1]]
    g2 <- groups[[2]]
    g3 <- groups[[3]]
    add_group_effect(1, "positive", max(4L, ceiling(length(g1) * 0.30)), c(0.06, 0.11))
    add_mixed_group(2,
                    max(3L, ceiling(length(g2) * 0.15)),
                    max(3L, ceiling(length(g2) * 0.15)))
    add_group_effect(3, "negative", max(4L, ceiling(length(g3) * 0.20)), c(0.07, 0.12))
  } else if (profile == "many5") {
    for (idx in seq_len(min(5L, length(groups)))) {
      direction <- if (idx %% 2L == 0L) "negative" else "positive"
      add_group_effect(idx, direction, 2L, c(0.09, 0.13))
    }
  }

  tibble(
    Variable = vars,
    True_Beta = as.numeric(beta),
    Group = rep(names(groups), lengths(groups))
  ) %>%
    mutate(
      True_Direction = case_when(
        True_Beta > 0.01 ~ "Positive",
        True_Beta < -0.01 ~ "Negative",
        TRUE ~ "None"
      ),
      IsActive = True_Direction != "None",
      IsActiveGroup = as.logical(ave(IsActive, Group, FUN = any))
    )
}

result_meta <- function(res, i) {
  p <- res$p
  tibble(
    ResultID = i,
    Scenario = res$Scenario %||% sprintf("p%03d", as.integer(p)),
    p = p,
    n = res$n %||% 1000L,
    n_groups = res$n_groups %||% 10L,
    truth_profile = res$truth_profile %||% "baseline3",
    within_rho = res$within_rho %||% NA_real_,
    cross_rho = res$cross_rho %||% NA_real_,
    Seed = res$seed
  )
}

summarize_results <- function(res, subdir, scenario_cols) {
  sub_out <- file.path(out_dir, subdir)
  dir.create(sub_out, showWarnings = FALSE, recursive = TRUE)

  meta <- purrr::imap_dfr(res, result_meta)
  truth_df <- meta %>%
    distinct(Scenario, p, n, n_groups, truth_profile, within_rho, cross_rho) %>%
    pmap_dfr(function(Scenario, p, n, n_groups, truth_profile, within_rho, cross_rho) {
      make_truth(p = p, n_groups = n_groups, truth_profile = truth_profile) %>%
        mutate(
          Scenario = Scenario,
          p = p,
          n = n,
          n_groups = n_groups,
          truth_profile = truth_profile,
          within_rho = within_rho,
          cross_rho = cross_rho
        )
    })

  accuracy_df <- purrr::imap_dfr(res, function(x, i) {
    m <- result_meta(x, i)
    purrr::imap_dfr(x$accuracy, function(acc, method) {
      bind_cols(m, tibble(Method = method, Active_Direction_Accuracy = as.numeric(acc)))
    })
  })

  timing_df <- purrr::imap_dfr(res, function(x, i) {
    m <- result_meta(x, i)
    purrr::imap_dfr(x$timing, function(elapsed, method) {
      bind_cols(m, tibble(Method = method, Seconds = as.numeric(elapsed)))
    })
  })

  weights_df <- purrr::imap_dfr(res, function(x, i) {
    m <- result_meta(x, i)
    bind_rows(x$weights) %>%
      mutate(
        Scenario = m$Scenario,
        n = m$n,
        n_groups = m$n_groups,
        truth_profile = m$truth_profile,
        within_rho = m$within_rho,
        cross_rho = m$cross_rho
      )
  })

  effects_df <- purrr::imap_dfr(res, function(x, i) {
    m <- result_meta(x, i)
    bind_rows(x$effects) %>%
      mutate(
        Scenario = m$Scenario,
        n = m$n,
        n_groups = m$n_groups,
        truth_profile = m$truth_profile,
        within_rho = m$within_rho,
        cross_rho = m$cross_rho
      )
  })

  error_raw <- purrr::imap_dfr(res, function(x, i) {
    if (!length(x$errors)) return(tibble())
    m <- result_meta(x, i)
    purrr::imap_dfr(x$errors, function(msg, method) {
      bind_cols(m, tibble(Method = method, Error = as.character(msg)))
    })
  })

  attempts <- meta %>%
    select(all_of(c(scenario_cols, "Seed"))) %>%
    distinct() %>%
    tidyr::crossing(Method = method_levels)

  if (nrow(error_raw) == 0L) {
    error_df <- attempts[0, ] %>% mutate(Error = character())
  } else {
    error_df <- error_raw
  }

  failure_summary <- attempts %>%
    left_join(error_df %>% mutate(Failed = TRUE),
              by = c(scenario_cols, "Seed", "Method")) %>%
    mutate(Failed = replace_na(Failed, FALSE)) %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      Attempts = n(),
      Failures = sum(Failed),
      Failure_Rate = Failures / Attempts,
      .groups = "drop"
    )

  accuracy_summary <- accuracy_df %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      Mean = mean(Active_Direction_Accuracy, na.rm = TRUE),
      SD = sd(Active_Direction_Accuracy, na.rm = TRUE),
      Median = median(Active_Direction_Accuracy, na.rm = TRUE),
      N = sum(!is.na(Active_Direction_Accuracy)),
      .groups = "drop"
    )

  timing_summary <- timing_df %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      Median_Seconds = median(Seconds, na.rm = TRUE),
      IQR_Seconds = IQR(Seconds, na.rm = TRUE),
      Mean_Seconds = mean(Seconds, na.rm = TRUE),
      N = sum(!is.na(Seconds)),
      .groups = "drop"
    )

  weight_share_df <- weights_df %>%
    group_by(Method, p, Seed, Scenario, Variable) %>%
    summarize(MaxWeight = max(Weight, na.rm = TRUE), .groups = "drop") %>%
    left_join(
      truth_df %>%
        select(Scenario, p, Variable, True_Beta, Group, True_Direction, IsActive, IsActiveGroup),
      by = c("Scenario", "p", "Variable")
    ) %>%
    group_by(Method, p, Seed, Scenario) %>%
    mutate(
      TotalWeight = sum(MaxWeight, na.rm = TRUE),
      NormWeight = if_else(TotalWeight > 0, MaxWeight / TotalWeight, NA_real_)
    ) %>%
    ungroup() %>%
    left_join(meta %>% distinct(Scenario, p, n, n_groups, truth_profile, within_rho, cross_rho),
              by = c("Scenario", "p"))

  active_weight_summary <- weight_share_df %>%
    group_by(across(all_of(scenario_cols)), Method, Seed, IsActive) %>%
    summarize(WeightShare = sum(NormWeight, na.rm = TRUE), .groups = "drop") %>%
    filter(IsActive) %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      Mean = mean(WeightShare, na.rm = TRUE),
      SD = sd(WeightShare, na.rm = TRUE),
      Median = median(WeightShare, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    )

  null_group_summary <- weight_share_df %>%
    group_by(across(all_of(scenario_cols)), Method, Seed, Group, IsActiveGroup) %>%
    summarize(GroupShare = sum(NormWeight, na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(scenario_cols)), Method, Seed, IsActiveGroup) %>%
    summarize(TotalShare = sum(GroupShare, na.rm = TRUE), .groups = "drop") %>%
    filter(!IsActiveGroup) %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      Mean = mean(TotalShare, na.rm = TRUE),
      SD = sd(TotalShare, na.rm = TRUE),
      Median = median(TotalShare, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    )

  effect_summary <- effects_df %>%
    group_by(across(all_of(scenario_cols)), Method, Effect_Type) %>%
    summarize(
      Mean = mean(Estimate, na.rm = TRUE),
      SD = sd(Estimate, na.rm = TRUE),
      Median = median(Estimate, na.rm = TRUE),
      N = sum(!is.na(Estimate)),
      .groups = "drop"
    )

  paired_tests <- accuracy_df %>%
    select(all_of(c(scenario_cols, "Seed")), Method, Active_Direction_Accuracy) %>%
    pivot_wider(names_from = Method, values_from = Active_Direction_Accuracy) %>%
    group_by(across(all_of(scenario_cols))) %>%
    group_modify(function(dat, key) {
      purrr::map_dfr(setdiff(method_levels, "SGL-WQS"), function(method) {
        d <- dat[["SGL-WQS"]] - dat[[method]]
        tibble(
          Comparison = paste("SGL-WQS -", method),
          N = sum(!is.na(d)),
          Mean_Diff = mean(d, na.rm = TRUE),
          Median_Diff = median(d, na.rm = TRUE),
          T_P = tryCatch(t.test(d, mu = 0)$p.value, error = function(e) NA_real_),
          Wilcoxon_P = tryCatch(
            wilcox.test(dat[["SGL-WQS"]], dat[[method]], paired = TRUE, exact = FALSE)$p.value,
            error = function(e) NA_real_
          )
        )
      })
    }) %>%
    ungroup() %>%
    mutate(
      T_P_BH = p.adjust(T_P, method = "BH"),
      Wilcoxon_P_BH = p.adjust(Wilcoxon_P, method = "BH")
    )

  one_sample_tests <- accuracy_df %>%
    group_by(across(all_of(scenario_cols)), Method) %>%
    summarize(
      N = n(),
      Mean = mean(Active_Direction_Accuracy, na.rm = TRUE),
      SD = sd(Active_Direction_Accuracy, na.rm = TRUE),
      T_P_gt_0_5 = tryCatch(
        t.test(Active_Direction_Accuracy, mu = 0.5, alternative = "greater")$p.value,
        error = function(e) NA_real_
      ),
      Wilcoxon_P_gt_0_5 = tryCatch(
        wilcox.test(Active_Direction_Accuracy, mu = 0.5, alternative = "greater", exact = FALSE)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      T_P_BH = p.adjust(T_P_gt_0_5, method = "BH"),
      Wilcoxon_P_BH = p.adjust(Wilcoxon_P_gt_0_5, method = "BH")
    )

  selection_df <- purrr::imap_dfr(res, function(x, i) {
    m <- result_meta(x, i)
    bind_rows(x$selection) %>%
      mutate(
        Scenario = m$Scenario,
        n = m$n,
        n_groups = m$n_groups,
        truth_profile = m$truth_profile,
        within_rho = m$within_rho,
        cross_rho = m$cross_rho
      )
  })

  if (nrow(selection_df) > 0L) {
    var_sel <- selection_df %>%
      group_by(across(all_of(scenario_cols)), Seed, Variable, IsActive) %>%
      summarize(MaxSelFreq = max(SelFreq, na.rm = TRUE), .groups = "drop")

    sgl_tnr_summary <- purrr::map_dfr(c(0.20, 0.50, 0.80), function(threshold) {
      var_sel %>%
        group_by(across(all_of(scenario_cols))) %>%
        summarize(
          Threshold = threshold,
          Seeds = n_distinct(Seed),
          Active_Rows = sum(IsActive),
          Null_Rows = sum(!IsActive),
          TPR = mean(MaxSelFreq[IsActive] >= threshold, na.rm = TRUE),
          TNR = mean(MaxSelFreq[!IsActive] < threshold, na.rm = TRUE),
          FPR = mean(MaxSelFreq[!IsActive] >= threshold, na.rm = TRUE),
          Active_Mean_MaxSelFreq = mean(MaxSelFreq[IsActive], na.rm = TRUE),
          Null_Mean_MaxSelFreq = mean(MaxSelFreq[!IsActive], na.rm = TRUE),
          .groups = "drop"
        )
    })
  } else {
    sgl_tnr_summary <- tibble()
  }

  write.csv(failure_summary, file.path(sub_out, "failure_summary.csv"), row.names = FALSE)
  write.csv(accuracy_summary, file.path(sub_out, "active_direction_accuracy_summary.csv"), row.names = FALSE)
  write.csv(active_weight_summary, file.path(sub_out, "active_weight_share_summary.csv"), row.names = FALSE)
  write.csv(null_group_summary, file.path(sub_out, "null_group_weight_share_summary.csv"), row.names = FALSE)
  write.csv(timing_summary, file.path(sub_out, "runtime_summary.csv"), row.names = FALSE)
  write.csv(effect_summary, file.path(sub_out, "effect_summary.csv"), row.names = FALSE)
  write.csv(error_df, file.path(sub_out, "error_detail.csv"), row.names = FALSE)
  write.csv(paired_tests, file.path(sub_out, "active_direction_accuracy_paired_tests_sgl_vs_others.csv"), row.names = FALSE)
  write.csv(one_sample_tests, file.path(sub_out, "active_direction_accuracy_one_sample_gt_0p5.csv"), row.names = FALSE)
  write.csv(sgl_tnr_summary, file.path(sub_out, "sglwqs_selection_tnr_summary.csv"), row.names = FALSE)

  list(
    failure = failure_summary,
    accuracy = accuracy_summary,
    active_weight = active_weight_summary,
    null_group = null_group_summary,
    timing = timing_summary,
    paired = paired_tests,
    one_sample = one_sample_tests,
    tnr = sgl_tnr_summary
  )
}

read_required_rds <- function(path) {
  if (!file.exists(path)) stop("Required RDS not found: ", path)
  readRDS(path)
}

p_res <- read_required_rds(file.path(
  "analysis", "analysis_cache", "highdim_p_checkpointed",
  "highdim_p_checkpointed_n1000_groups10_boot200_seeds30.rds"
))
group_res <- read_required_rds(file.path(
  "analysis", "analysis_cache", "highdim_group_truth_checkpointed",
  "highdim_group_truth_checkpointed_boot100_seeds10.rds"
))
corr_res <- read_required_rds(file.path(
  "analysis", "analysis_cache", "highdim_correlation_checkpointed",
  "highdim_correlation_checkpointed_boot100_seeds10.rds"
))

p_sum <- summarize_results(p_res, "p_sensitivity", c("p"))
group_sum <- summarize_results(group_res, "group_signal_sensitivity",
                               c("Scenario", "p", "n", "n_groups", "truth_profile"))
corr_sum <- summarize_results(corr_res, "correlation_signal_sensitivity",
                              c("Scenario", "p", "n", "n_groups", "truth_profile", "within_rho", "cross_rho"))

fmt_num <- function(x, digits = 3) ifelse(is.na(x), NA, round(x, digits))

make_wide_metric <- function(df, id_cols, value_col = "Mean") {
  df %>%
    mutate(Method = factor(Method, levels = method_levels)) %>%
    select(all_of(id_cols), Method, all_of(value_col)) %>%
    mutate("{value_col}" := fmt_num(.data[[value_col]])) %>%
    pivot_wider(names_from = Method, values_from = all_of(value_col)) %>%
    arrange(across(all_of(id_cols)))
}

p_accuracy_wide <- make_wide_metric(p_sum$accuracy, "p")
p_weight_wide <- make_wide_metric(p_sum$active_weight, "p")
p_runtime_wide <- make_wide_metric(p_sum$timing, "p", "Median_Seconds")
group_accuracy_wide <- make_wide_metric(group_sum$accuracy, c("Scenario", "n_groups", "truth_profile"))
corr_accuracy_wide <- make_wide_metric(corr_sum$accuracy, c("Scenario", "within_rho", "truth_profile"))
p_tnr_wide <- p_sum$tnr %>%
  filter(Threshold %in% c(0.50, 0.80)) %>%
  mutate(across(c(TPR, TNR, FPR, Active_Mean_MaxSelFreq, Null_Mean_MaxSelFreq), fmt_num)) %>%
  arrange(p, Threshold)

write.csv(p_accuracy_wide, file.path(out_dir, "summary_table_p_active_direction_accuracy_wide.csv"), row.names = FALSE)
write.csv(p_weight_wide, file.path(out_dir, "summary_table_p_active_weight_share_wide.csv"), row.names = FALSE)
write.csv(p_runtime_wide, file.path(out_dir, "summary_table_p_runtime_median_seconds_wide.csv"), row.names = FALSE)
write.csv(group_accuracy_wide, file.path(out_dir, "summary_table_group_signal_accuracy_wide.csv"), row.names = FALSE)
write.csv(corr_accuracy_wide, file.path(out_dir, "summary_table_correlation_signal_accuracy_wide.csv"), row.names = FALSE)
write.csv(p_tnr_wide, file.path(out_dir, "summary_table_p_sglwqs_tnr.csv"), row.names = FALSE)

md <- c(
  "# High-Dimensional Simulation Summary",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "- p-sensitivity: n=1000, groups=10, p=50/100/200, 30 seeds.",
  "- group/signal sensitivity: p=100, n=1000, group count and signal strength perturbations, 70 scenario-seeds.",
  "- correlation/signal sensitivity: p=100, n=1000, groups=10, within-group rho=0.20/0.45/0.80 and signal strength perturbations, 70 scenario-seeds.",
  "",
  "## Metric Notes",
  "",
  "- Active direction accuracy is direction-only among truly active components. It does not require component-level GLM significance.",
  "- SGL-WQS TNR is based on variable-level max positive/negative bootstrap selection frequency below a threshold among null variables.",
  "- qgcomp and gWQS produce non-sparse weights for all components; sparse TNR is therefore reported only for SGL-WQS.",
  "",
  "## Main Interpretation",
  "",
  "- All high-dimensional simulations completed without method failures.",
  "- SGL-WQS remained computationally stable through p=200 and across group/signal/correlation perturbations.",
  "- Active direction recovery was generally above chance, but SGL-WQS did not show adjusted statistically significant superiority over comparators in the p-sensitivity analysis.",
  "- Active-variable weight concentration decreased as dimensionality increased.",
  "- SGL-WQS bootstrap selection frequencies separated active and null components weakly in high-dimensional settings; null components often had high selection frequencies.",
  "- Overall, these results support feasibility and robustness against outright failure rather than strong high-dimensional component-selection performance.",
  "",
  "## p-Sensitivity: Active Direction Accuracy",
  "",
  capture.output(knitr::kable(p_accuracy_wide, format = "pipe")),
  "",
  "## p-Sensitivity: Active Weight Share",
  "",
  capture.output(knitr::kable(p_weight_wide, format = "pipe")),
  "",
  "## p-Sensitivity: Median Runtime Seconds",
  "",
  capture.output(knitr::kable(p_runtime_wide, format = "pipe")),
  "",
  "## p-Sensitivity: SGL-WQS TPR/TNR by Selection-Frequency Threshold",
  "",
  capture.output(knitr::kable(p_tnr_wide, format = "pipe")),
  "",
  "## Group/Signal Sensitivity: Active Direction Accuracy",
  "",
  capture.output(knitr::kable(group_accuracy_wide, format = "pipe")),
  "",
  "## Correlation/Signal Sensitivity: Active Direction Accuracy",
  "",
  capture.output(knitr::kable(corr_accuracy_wide, format = "pipe")),
  "",
  "## Primary Output Files",
  "",
  "- `p_sensitivity/`: p=50/100/200 summaries.",
  "- `group_signal_sensitivity/`: group count and signal strength summaries.",
  "- `correlation_signal_sensitivity/`: within-group correlation and signal strength summaries.",
  "- `summary_table_*.csv`: compact cross-method tables used in this report."
)

writeLines(md, file.path(out_dir, "HIGH_DIMENSIONAL_SIMULATION_SUMMARY.md"))

manifest <- tibble(
  file = list.files(out_dir, recursive = TRUE, full.names = FALSE),
  description = case_when(
    grepl("HIGH_DIMENSIONAL_SIMULATION_SUMMARY", file) ~ "Human-readable markdown summary",
    grepl("failure_summary", file) ~ "Method failure counts and rates",
    grepl("active_direction_accuracy_summary", file) ~ "Active direction accuracy summary",
    grepl("active_weight_share_summary", file) ~ "Normalized weight share assigned to active components",
    grepl("null_group_weight_share_summary", file) ~ "Normalized weight share assigned to null groups",
    grepl("runtime_summary", file) ~ "Runtime summary in seconds",
    grepl("effect_summary", file) ~ "Mixture/index effect estimate summary",
    grepl("paired_tests", file) ~ "Paired tests comparing SGL-WQS active direction accuracy with comparators",
    grepl("one_sample", file) ~ "One-sample tests for active direction accuracy greater than 0.5",
    grepl("sglwqs_selection_tnr", file) ~ "SGL-WQS TPR/TNR/FPR by selection-frequency threshold",
    grepl("summary_table", file) ~ "Compact summary table used in markdown report",
    grepl("error_detail", file) ~ "Raw error messages; empty if no failures",
    TRUE ~ "High-dimensional simulation summary output"
  )
)
write.csv(manifest, file.path(out_dir, "manifest.csv"), row.names = FALSE)

cat("Wrote high-dimensional overall summary to:", out_dir, "\n")
cat("Markdown:", file.path(out_dir, "HIGH_DIMENSIONAL_SIMULATION_SUMMARY.md"), "\n")
