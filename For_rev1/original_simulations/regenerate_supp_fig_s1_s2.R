#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  normalizePath("final_res")
}
repo_dir <- normalizePath(file.path(script_dir, ".."))
out_dir <- file.path(script_dir, "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cache_dir <- file.path(
  script_dir,
  "compare_mixture_methods_by_samplesize_30seeds_cache",
  "html"
)
cache_rdb <- list.files(
  cache_dir,
  pattern = "^run-analysis_.*[.]rdb$",
  full.names = TRUE
)
if (length(cache_rdb) == 0) {
  stop("Could not find the run-analysis knitr cache in: ", cache_dir)
}
cache_base <- sub("[.]rdb$", "", cache_rdb[[1]])
cache_env <- new.env(parent = emptyenv())
lazyLoad(cache_base, envir = cache_env)
all_results <- cache_env$all_results

if (!is.list(all_results) || length(all_results) == 0) {
  stop("Cached all_results object was not found or was empty.")
}

sample_sizes <- c(100, 500, 1000, 5000, 10000, 50000)
method_order <- c("SGL-WQS", "gWQS", "groupWQS", "qgcomp")
# Match the method-color mapping used by main Figure 3
# (scale_fill_brewer(palette = "Set1") in the correlation robustness Rmd).
method_palette <- c(
  "SGL-WQS" = "#984EA3",
  "gWQS" = "#377EB8",
  "groupWQS" = "#E41A1C",
  "qgcomp" = "#4DAF4A"
)
direction_palette <- c(
  "Positive" = "#0072B2",
  "Negative" = "#D55E00",
  "Combined" = "#6A737D"
)

pcbs <- c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_187")
metals <- c("cd", "pb", "hg", "se", "mg")
pfass <- c("PFOA", "PFNA", "PFOS")
all_exposures <- c(pcbs, metals, pfass)

true_directions <- tibble(
  Variable = all_exposures,
  True_Beta = c(
    0.04, 0.15, 0.20, 0.08, 0.05,
    -0.12, -0.18, 0.10, 0.14, -0.02,
    0, 0, 0
  ),
  Group = c(rep("PCBs", 5), rep("Metals", 5), rep("PFASs", 3))
) |>
  mutate(
    True_Direction = case_when(
      True_Beta > 0.01 ~ "Positive",
      True_Beta < -0.01 ~ "Negative",
      TRUE ~ "None"
    )
  )

theme_supp <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3, margin = margin(b = 5)),
      plot.subtitle = element_text(size = base_size, margin = margin(b = 8)),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(12, 14, 12, 14)
    )
}

save_plot <- function(plot, name, width, height) {
  pdf_file <- file.path(out_dir, paste0(name, ".pdf"))
  png_file <- file.path(out_dir, paste0(name, ".png"))

  ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::pdf,
    useDingbats = FALSE
  )
  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = ragg::agg_png,
    bg = "white"
  )

  invisible(c(pdf_file, png_file))
}

make_nhanes_dag_plot <- function() {
  dag_nodes <- tibble(
    id = c(
      "age", "sex", "race", "pir", "cotinine", "bmi", "ucr",
      "voc", "ggt"
    ),
    label = c(
      "Age", "Sex", "Race/\nethnicity", "Family income-\nto-poverty ratio",
      "Log cotinine\n(smoking)", "BMI", "Log urinary\ncreatinine",
      "Urinary VOC\nmetabolite mixture", "Log-GGT"
    ),
    x = c(1.15, 2.95, 1.15, 2.95, 1.15, 2.95, 2.05, 5.9, 8.55),
    y = c(4.0, 4.0, 3.1, 3.1, 2.15, 2.15, 1.15, 2.75, 2.75),
    type = c(
      rep("DAG-based adjustment variables", 6),
      "Urinary dilution adjustment",
      "Exposure mixture",
      "Outcome"
    )
  )

  internal_edges <- tibble(
    x = c(1.55, 2.95, 1.55, 2.95, 2.6),
    y = c(3.1, 3.1, 4.0, 4.0, 1.95),
    xend = c(2.35, 1.35, 2.55, 2.8, 2.25),
    yend = c(3.1, 2.35, 2.35, 2.35, 1.35)
  )

  ggplot() +
    annotate(
      "rect",
      xmin = 0.35, xmax = 3.75, ymin = 0.55, ymax = 4.65,
      fill = "#F7F8FA", color = "#BAC3CF", linewidth = 0.55
    ) +
    annotate(
      "text",
      x = 2.05, y = 4.92,
      label = "DAG-based adjustment variables",
      size = 4.3, fontface = "bold", color = "#202A33"
    ) +
    geom_segment(
      data = internal_edges,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "#AEB7C2",
      linewidth = 0.45,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.12, "in"))
    ) +
    geom_curve(
      aes(x = 3.72, y = 3.3, xend = 5.05, yend = 2.95),
      curvature = -0.14,
      color = "#5E6A76",
      linewidth = 0.75,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.16, "in"))
    ) +
    geom_curve(
      aes(x = 3.72, y = 3.95, xend = 8.08, yend = 3.08),
      curvature = -0.22,
      color = "#5E6A76",
      linewidth = 0.75,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.16, "in"))
    ) +
    geom_curve(
      aes(x = 2.55, y = 1.15, xend = 5.05, yend = 2.45),
      curvature = 0.18,
      color = "#B66D00",
      linetype = "dashed",
      linewidth = 0.7,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.16, "in"))
    ) +
    geom_segment(
      aes(x = 6.72, y = 2.75, xend = 7.85, yend = 2.75),
      color = "#1F78B4",
      linewidth = 1.1,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.18, "in"))
    ) +
    annotate(
      "label",
      x = 4.55, y = 3.55,
      label = "potential\ncommon causes",
      size = 3.4,
      linewidth = 0,
      fill = "white",
      color = "#4E5965"
    ) +
    annotate(
      "label",
      x = 4.35, y = 1.65,
      label = "urinary dilution\nadjustment",
      size = 3.3,
      linewidth = 0,
      fill = "white",
      color = "#8A5200"
    ) +
    annotate(
      "label",
      x = 7.25, y = 2.25,
      label = "target\ncontrast",
      size = 3.3,
      linewidth = 0,
      fill = "white",
      color = "#1F5F8B"
    ) +
    geom_label(
      data = dag_nodes,
      aes(x = x, y = y, label = label, fill = type),
      color = "#111820",
      fontface = "bold",
      size = 3.6,
      linewidth = 0.45,
      label.r = grid::unit(0.12, "lines"),
      label.padding = grid::unit(0.28, "lines")
    ) +
    scale_fill_manual(
      values = c(
        "DAG-based adjustment variables" = "#FFFFFF",
        "Urinary dilution adjustment" = "#FFF2CC",
        "Exposure mixture" = "#DDEEFF",
        "Outcome" = "#E3F2E8"
      )
    ) +
    labs(
      title = "Simplified DAG for the NHANES VOC demonstration",
      subtitle = "The primary covariate set adjusts for DAG-based common causes; urinary creatinine is included as a dilution adjustment"
    ) +
    coord_cartesian(xlim = c(0, 9.5), ylim = c(0.35, 5.4), clip = "off") +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.02, margin = margin(b = 4)),
      plot.subtitle = element_text(size = 12.5, hjust = 0.02, color = "#3D4752", margin = margin(b = 8)),
      legend.position = "none",
      plot.margin = margin(14, 18, 12, 18)
    )
}

nhanes_dag <- make_nhanes_dag_plot()
save_plot(
  nhanes_dag,
  "Supplementary_Figure_S1_NHANES_VOC_DAG",
  width = 10.2,
  height = 5.9
)

accuracy_rows <- lapply(all_results, function(res) {
  rows <- list()

  if (!is.null(res$qgcomp_accuracy)) {
    rows[[length(rows) + 1L]] <- tibble(
      n = res$n,
      Seed = res$seed,
      Method = "qgcomp",
      Accuracy = res$qgcomp_accuracy
    )
  }
  if (!is.null(res$gwqs_accuracy)) {
    rows[[length(rows) + 1L]] <- tibble(
      n = res$n,
      Seed = res$seed,
      Method = "gWQS",
      Accuracy = res$gwqs_accuracy
    )
  }
  if (!is.null(res$groupwqs_coefs) && length(res$groupwqs_coefs) >= 3) {
    coefs <- res$groupwqs_coefs
    group_dirs <- c(
      setNames(rep(ifelse(coefs[[1]] >= 0, "Positive", "Negative"), length(pcbs)), pcbs),
      setNames(rep(ifelse(coefs[[2]] >= 0, "Positive", "Negative"), length(metals)), metals),
      setNames(rep(ifelse(coefs[[3]] >= 0, "Positive", "Negative"), length(pfass)), pfass)
    )
    group_check <- tibble(
      Variable = names(group_dirs),
      Detected = as.character(group_dirs)
    ) |>
      left_join(true_directions |> select(Variable, True_Direction), by = "Variable")
    group_accuracy <- sum(
      group_check$Detected == group_check$True_Direction |
        group_check$True_Direction == "None"
    ) / nrow(group_check)
    rows[[length(rows) + 1L]] <- tibble(
      n = res$n,
      Seed = res$seed,
      Method = "groupWQS",
      Accuracy = group_accuracy
    )
  }
  if (!is.null(res$sglwqs_accuracy)) {
    rows[[length(rows) + 1L]] <- tibble(
      n = res$n,
      Seed = res$seed,
      Method = "SGL-WQS",
      Accuracy = res$sglwqs_accuracy
    )
  }

  bind_rows(rows)
})

accuracy_df <- bind_rows(accuracy_rows) |>
  mutate(
    n = factor(n, levels = sample_sizes),
    n_numeric = as.numeric(as.character(n)),
    Method = factor(Method, levels = method_order)
  )

accuracy_summary <- accuracy_df |>
  group_by(n, n_numeric, Method) |>
  summarize(
    Mean = mean(Accuracy, na.rm = TRUE),
    SD = sd(Accuracy, na.rm = TRUE),
    Successful_Seeds = dplyr::n(),
    .groups = "drop"
  ) |>
  arrange(n_numeric, Method)

readr::write_csv(
  accuracy_summary |>
    mutate(
      Mean_Percent = Mean * 100,
      SD_Percent = SD * 100
    ),
  file.path(out_dir, "Supplementary_Figure_S2_accuracy_summary.csv")
)

s1 <- ggplot(
  accuracy_summary,
  aes(x = n_numeric, y = Mean, color = Method, group = Method)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.8) +
  geom_errorbar(
    aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
    width = 0.04,
    linewidth = 0.65,
    alpha = 0.85
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.4) +
  scale_x_log10(
    breaks = sample_sizes,
    labels = comma_format()
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.05),
    breaks = seq(0, 1, 0.25),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_color_manual(values = method_palette, breaks = method_order, drop = FALSE) +
  labs(
    title = "Direction detection accuracy across sample sizes",
    subtitle = "Mean ± SD across successful fits; 30 seeds were run per scenario",
    x = "Sample size (log scale)",
    y = "Accuracy"
  ) +
  theme_supp(base_size = 14) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

save_plot(
  s1,
  "Supplementary_Figure_S2_Direction_Detection_Accuracy",
  width = 7.5,
  height = 5.8
)

weight_rows <- lapply(all_results, function(res) {
  bind_rows(
    res$qgcomp_weights,
    res$gwqs_weights,
    res$groupwqs_weights,
    res$sglwqs_weights
  )
})
all_weights <- bind_rows(weight_rows) |>
  filter(!is.na(Weight)) |>
  left_join(true_directions, by = "Variable") |>
  mutate(
    Method = factor(Method, levels = method_order),
    Variable = factor(Variable, levels = all_exposures),
    Sample_Size = factor(
      paste0("n = ", comma(n)),
      levels = paste0("n = ", comma(sample_sizes))
    ),
    Direction = factor(Direction, levels = c("Positive", "Negative", "Combined"))
  )

readr::write_csv(
  all_weights |>
    group_by(Method, n, Variable, Direction) |>
    summarize(
      Mean_Weight = mean(Weight, na.rm = TRUE),
      SD_Weight = sd(Weight, na.rm = TRUE),
      Successful_Seeds = dplyr::n_distinct(Seed),
      .groups = "drop"
    ),
  file.path(out_dir, "Supplementary_Figure_S3_weight_summary.csv")
)

method_labels <- c(
  "SGL-WQS" = "SGL-WQS",
  "gWQS" = "gWQS",
  "groupWQS" = "groupWQS",
  "qgcomp" = "qgcomp"
)
panel_labels <- c(
  "SGL-WQS" = "S3A",
  "gWQS" = "S3B",
  "groupWQS" = "S3C",
  "qgcomp" = "S3D"
)
panel_titles <- c(
  "SGL-WQS" = "S3A. SGL-WQS weights by sample size",
  "gWQS" = "S3B. gWQS weights by sample size",
  "groupWQS" = "S3C. groupWQS weights by sample size",
  "qgcomp" = "S3D. qgcomp coefficient-derived attribution by sample size"
)
panel_subtitles <- c(
  "SGL-WQS" = "Violin plots and dots show seed-level component weights; diamonds indicate means",
  "gWQS" = "Violin plots and dots show seed-level component weights; diamonds indicate means",
  "groupWQS" = "Violin plots and dots show seed-level component weights; diamonds indicate means",
  "qgcomp" = "Violin plots and dots show seed-level post-hoc attribution; diamonds indicate means"
)
file_labels <- c(
  "SGL-WQS" = "SGL-WQS",
  "gWQS" = "gWQS",
  "groupWQS" = "groupWQS",
  "qgcomp" = "qgcomp"
)

make_s2_plot <- function(method_name) {
  plot_data <- all_weights |>
    filter(Method == method_name)

  if (nrow(plot_data) == 0) {
    stop("No weight data found for method: ", method_name)
  }

  violin_data <- plot_data |>
    add_count(Sample_Size, Variable, Direction, name = "Group_N") |>
    filter(Group_N >= 2)

  ggplot(
    plot_data,
    aes(x = Variable, y = Weight, fill = Direction)
  ) +
    geom_vline(
      xintercept = c(length(pcbs) + 0.5, length(pcbs) + length(metals) + 0.5),
      linetype = "dotted",
      color = "grey70",
      linewidth = 0.4
    ) +
    geom_violin(
      data = violin_data,
      alpha = 0.32,
      width = 0.82,
      linewidth = 0.35,
      trim = FALSE,
      scale = "width",
      position = position_dodge(width = 0.78)
    ) +
    geom_point(
      size = 0.55,
      alpha = 0.38,
      stroke = 0,
      position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.78)
    ) +
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 23,
      fill = "white",
      color = "black",
      stroke = 0.35,
      size = 1.7,
      position = position_dodge(width = 0.78)
    ) +
    facet_wrap(~ Sample_Size, ncol = 2) +
    scale_fill_manual(values = direction_palette, drop = TRUE) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.25),
      expand = expansion(mult = c(0, 0.04))
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = panel_titles[[method_name]],
      subtitle = panel_subtitles[[method_name]],
      x = NULL,
      y = ifelse(method_name == "qgcomp", "Attribution", "Weight")
    ) +
    theme_supp(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      panel.grid.major.x = element_blank()
    )
}

for (method_name in method_order) {
  s2_plot <- make_s2_plot(method_name)
  save_plot(
    s2_plot,
    paste0(
      "Supplementary_Figure_",
      panel_labels[[method_name]],
      "_",
      file_labels[[method_name]],
      "_Weights_by_Sample_Size"
    ),
    width = 13.5,
    height = 9.5
  )
}

cat("Regenerated Supplementary Figure S1, S2, and S3A-D in:\n")
cat(out_dir, "\n")
cat("Method order: ", paste(method_order, collapse = ", "), "\n", sep = "")
