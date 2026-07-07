suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

cache_dir <- file.path(
  "analysis",
  "analysis_cache",
  "highdim_p_checkpointed",
  "scenarios_n1000_groups10_boot200"
)
out_dir <- file.path("analysis", "results", "highdim_wrong_direction")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_group_definitions <- function(p, n_groups = 10) {
  vars <- sprintf("X%03d", seq_len(p))
  group_ids <- cut(seq_len(p),
                   breaks = n_groups,
                   labels = sprintf("G%02d", seq_len(n_groups)),
                   include.lowest = TRUE)
  split(vars, group_ids)
}

make_truth_df <- function(p, n_groups = 10) {
  groups <- make_group_definitions(p, n_groups)
  vars <- unlist(groups, use.names = FALSE)
  beta <- setNames(rep(0, length(vars)), vars)

  g1 <- groups[[1]]
  n_g1 <- min(length(g1), max(3L, ceiling(length(g1) * 0.15)))
  beta[g1[seq_len(n_g1)]] <- seq(0.08, 0.12, length.out = n_g1)

  g2 <- groups[[2]]
  n_g2_neg <- min(length(g2), max(2L, floor(length(g2) * 0.08)))
  remaining_g2 <- max(length(g2) - n_g2_neg, 0L)
  n_g2_pos <- min(remaining_g2, max(2L, floor(length(g2) * 0.06)))
  if (n_g2_neg > 0L) {
    beta[g2[seq_len(n_g2_neg)]] <- seq(-0.16, -0.11, length.out = n_g2_neg)
  }
  if (n_g2_pos > 0L) {
    pos_idx <- seq.int(n_g2_neg + 1L, n_g2_neg + n_g2_pos)
    beta[g2[pos_idx]] <- seq(0.10, 0.14, length.out = n_g2_pos)
  }

  g3 <- groups[[3]]
  n_g3 <- min(length(g3), max(2L, ceiling(length(g3) * 0.08)))
  beta[g3[seq_len(n_g3)]] <- seq(-0.14, -0.10, length.out = n_g3)

  tibble(
    p = p,
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
      IsActiveGroup = Group %in% names(groups)[1:3]
    )
}

scenario_files <- list.files(cache_dir, pattern = "^p[0-9]+_seed[0-9]+[.]rds$", full.names = TRUE)
if (!length(scenario_files)) {
  stop("No checkpointed high-dimensional scenario files were found in: ", cache_dir)
}

scenario_results <- lapply(scenario_files, readRDS)
weights_df <- purrr::map_dfr(scenario_results, function(res) dplyr::bind_rows(res$weights))
selection_df <- purrr::map_dfr(scenario_results, function(res) dplyr::bind_rows(res$selection))

completed <- tibble(
  p = vapply(scenario_results, `[[`, numeric(1), "p"),
  Seed = vapply(scenario_results, `[[`, numeric(1), "seed")
) %>%
  count(p, name = "Completed_Scenarios") %>%
  arrange(p)

truth_df <- purrr::map_dfr(sort(unique(completed$p)), make_truth_df)

wrong_direction_weight_detail <- weights_df %>%
  filter(Direction %in% c("Positive", "Negative")) %>%
  group_by(Method, p, Seed, Variable, Direction) %>%
  summarize(Weight = sum(Weight, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Direction, values_from = Weight, values_fill = 0) %>%
  left_join(truth_df, by = c("p", "Variable")) %>%
  filter(IsActive) %>%
  mutate(
    CorrectWeight = if_else(True_Direction == "Positive", Positive, Negative),
    WrongWeight = if_else(True_Direction == "Positive", Negative, Positive),
    WrongDominant = WrongWeight > CorrectWeight,
    WrongShare = if_else(
      CorrectWeight + WrongWeight > 0,
      WrongWeight / (CorrectWeight + WrongWeight),
      NA_real_
    )
  )

wrong_direction_weight_summary <- wrong_direction_weight_detail %>%
  group_by(Method, p) %>%
  summarize(
    Completed_Scenarios = n_distinct(Seed),
    Active_Rows = n(),
    Wrong_Dominant_Rate = mean(WrongDominant, na.rm = TRUE),
    Mean_Wrong_Share = mean(WrongShare, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(p, Method)

write.csv(
  wrong_direction_weight_detail,
  file.path(out_dir, "highdim_wrong_direction_weight_detail.csv"),
  row.names = FALSE
)
write.csv(
  wrong_direction_weight_summary,
  file.path(out_dir, "highdim_wrong_direction_weight_summary.csv"),
  row.names = FALSE
)

cat("Completed scenarios:\n")
print(completed, n = Inf)
cat("\nWrong direction by positive/negative weights among active variables:\n")
print(wrong_direction_weight_summary, n = Inf)

if (nrow(selection_df)) {
  wrong_direction_selection_detail <- selection_df %>%
    filter(Direction %in% c("Positive", "Negative"), IsActive) %>%
    select(Method, p, Seed, Variable, Direction, SelFreq, True_Direction) %>%
    tidyr::pivot_wider(names_from = Direction, values_from = SelFreq, values_fill = 0) %>%
    mutate(
      CorrectSelFreq = if_else(True_Direction == "Positive", Positive, Negative),
      WrongSelFreq = if_else(True_Direction == "Positive", Negative, Positive),
      WrongDominant = WrongSelFreq > CorrectSelFreq,
      WrongSelFreq_ge_0_50 = WrongSelFreq >= 0.50,
      WrongSelFreq_ge_0_80 = WrongSelFreq >= 0.80
    )

  wrong_direction_selection_summary <- wrong_direction_selection_detail %>%
    group_by(Method, p) %>%
    summarize(
      Completed_Scenarios = n_distinct(Seed),
      Active_Rows = n(),
      Wrong_Dominant_Rate = mean(WrongDominant, na.rm = TRUE),
      Wrong_SelFreq_ge_0_50 = mean(WrongSelFreq_ge_0_50, na.rm = TRUE),
      Wrong_SelFreq_ge_0_80 = mean(WrongSelFreq_ge_0_80, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(p, Method)

  write.csv(
    wrong_direction_selection_detail,
    file.path(out_dir, "highdim_wrong_direction_sglwqs_selection_detail.csv"),
    row.names = FALSE
  )
  write.csv(
    wrong_direction_selection_summary,
    file.path(out_dir, "highdim_wrong_direction_sglwqs_selection_summary.csv"),
    row.names = FALSE
  )

  cat("\nSGL-WQS wrong direction by bootstrap selection frequency among active variables:\n")
  print(wrong_direction_selection_summary, n = Inf)
}

cat("\nWrote wrong-direction summaries to: ", out_dir, "\n", sep = "")
