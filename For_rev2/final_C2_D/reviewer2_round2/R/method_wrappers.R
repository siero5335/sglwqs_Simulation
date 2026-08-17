capture_fit <- function(expr) {
  warnings_seen <- character(0)
  start <- Sys.time()
  value <- tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  list(
    ok = !inherits(value, "error"),
    value = if (inherits(value, "error")) NULL else value,
    error = if (inherits(value, "error")) value else NULL,
    warnings = unique(warnings_seen),
    elapsed = elapsed
  )
}

qa_generated_data <- function(dat, job) {
  groups <- attr(dat, "groups")
  truth <- attr(dat, "truth")
  vars <- unlist(groups, use.names = FALSE)
  assert_or_stop(sum(lengths(groups)) == length(vars), "Group size metadata is inconsistent.")
  assert_or_stop(all(vars %in% names(dat)), "Not all exposure variables are present in generated data.")
  assert_or_stop(all(truth$Variable %in% vars), "Truth table contains variables outside groups.")
  null_groups <- names(which(vapply(groups, function(g) all(truth$True_Beta[match(g, truth$Variable)] == 0), logical(1))))
  if (length(null_groups)) {
    null_truth <- truth[truth$Group %in% null_groups, , drop = FALSE]
    assert_or_stop(all(null_truth$True_Beta == 0), "A null group has nonzero true coefficients.")
  }
  if (identical(job$family, "binomial")) {
    assert_or_stop(length(unique(dat$Y)) == 2L, "Binary outcome does not contain both classes.")
  } else {
    assert_or_stop(all(is.finite(dat$Y)), "Gaussian outcome contains non-finite values.")
    assert_or_stop(identical(attr(dat, "family"), "gaussian"), "Gaussian family metadata is missing.")
    assert_or_stop(isTRUE(all.equal(attr(dat, "gaussian_residual_sd"), 1)), "Gaussian residual SD is not 1.")
    assert_or_stop(length(attr(dat, "eta_linear")) == nrow(dat), "Gaussian identity linear predictor is missing.")
  }
  expected_direction <- ifelse(truth$True_Beta > 0, "Positive",
    ifelse(truth$True_Beta < 0, "Negative", "None")
  )
  assert_or_stop(all(truth$True_Direction == expected_direction), "Truth direction coding is inconsistent with beta signs.")
  if (job$battery %in% c("D_unbalanced_effect_factorial", "E_hard_setting_sample_size")) {
    expected_sizes <- factorial_group_sizes(job$group_structure)
    assert_or_stop(identical(as.integer(lengths(groups)), as.integer(expected_sizes)), "Factorial group sizes are incorrect.")
    active_counts <- vapply(groups, function(g) sum(truth$True_Beta[match(g, truth$Variable)] != 0), integer(1))
    assert_or_stop(identical(unname(active_counts[c("G01", "G05", "G10")]), rep(3L, 3L)), "Each active factorial group must contain three active components.")
    assert_or_stop(all(active_counts[setdiff(names(groups), c("G01", "G05", "G10"))] == 0L), "A nominally null factorial group is active.")
    direction_sets <- lapply(c("G01", "G05", "G10"), function(g) {
      unique(truth$True_Direction[truth$Group == g & truth$True_Direction != "None"])
    })
    assert_or_stop(
      identical(direction_sets[[1]], "Positive") &&
        identical(direction_sets[[2]], "Negative") &&
        setequal(direction_sets[[3]], c("Positive", "Negative")),
      "Factorial active-group direction plan is incorrect."
    )
    if (identical(job$family, "binomial")) {
      target_prevalence <- attr(dat, "binomial_target_prevalence")
      expected_prevalence <- attr(dat, "binomial_expected_prevalence")
      binomial_intercept <- attr(dat, "binomial_intercept")
      assert_or_stop(
        is.finite(target_prevalence) && is.finite(expected_prevalence) && is.finite(binomial_intercept),
        "Calibrated binomial prevalence metadata is missing."
      )
      assert_or_stop(
        abs(expected_prevalence - target_prevalence) < 1e-10,
        "Expected binomial prevalence was not calibrated to its target."
      )
      assert_or_stop(
        abs(mean(dat$Y) - target_prevalence) < 0.12,
        "Realized binomial prevalence is implausibly far from its calibrated target."
      )
    }
  }
  TRUE
}

generate_job_data <- function(job) {
  if (job$battery %in% c("A1_gaussian_sample_size", "A2_weak_gaussian_sample_size")) {
    return(generate_legacy13_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      family = job$family,
      signal_profile = job$signal_profile
    ))
  }
  if (job$battery %in% c("B1_gaussian_highdim", "B2_weak_gaussian_highdim")) {
    return(generate_highdim_p_data_matched(
      p = as.integer(job$p),
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      family = job$family,
      signal_profile = job$signal_profile
    ))
  }
  if (job$battery == "C_gaussian_split_stability") {
    return(generate_validation_split_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      family = job$family,
      effect_profile = job$effect_profile,
      within_r = as.numeric(job$within_rho),
      between_r = as.numeric(job$cross_rho)
    ))
  }
  if (job$battery == "C2_gaussian_paired_split_stability") {
    return(generate_paired_gaussian_split_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      split_seed = as.integer(job$split_seed),
      effect_profile = job$effect_profile,
      within_r = as.numeric(job$within_rho),
      between_r = as.numeric(job$cross_rho)
    ))
  }
  if (job$battery %in% c("D_unbalanced_effect_factorial", "E_hard_setting_sample_size")) {
    binomial_target_prevalence <- job$binomial_target_prevalence %||%
      factorial_binomial_target_prevalence()
    if (!is.finite(binomial_target_prevalence)) {
      binomial_target_prevalence <- factorial_binomial_target_prevalence()
    }
    return(generate_factorial_data(
      n = as.integer(job$n),
      data_seed = as.integer(job$data_seed),
      family = job$family,
      group_structure = job$group_structure,
      effect_profile = job$effect_profile,
      within_rho = as.numeric(job$within_rho),
      cross_rho = as.numeric(job$cross_rho),
      binomial_target_prevalence = as.numeric(binomial_target_prevalence)
    ))
  }
  stop("Unknown battery: ", job$battery, call. = FALSE)
}

truth_group_direction_rows <- function(truth, job) {
  groups <- split(truth, truth$Group)
  do.call(rbind, lapply(names(groups), function(group) {
    g <- groups[[group]]
    data.frame(
      scenario_id = job$scenario_id,
      battery = job$battery,
      family = job$family,
      method = job$method,
      data_seed = as.integer(job$data_seed),
      split_replicate = as.integer(job$split_replicate %||% NA_integer_),
      n = as.integer(job$n),
      p = as.integer(job$p),
      group = group,
      direction = c("positive", "negative"),
      is_true_null = c(!any(g$True_Beta > 0), !any(g$True_Beta < 0)),
      active_variable_count = c(sum(g$True_Beta > 0), sum(g$True_Beta < 0)),
      true_beta_sum = c(sum(g$True_Beta[g$True_Beta > 0]), sum(g$True_Beta[g$True_Beta < 0])),
      group_size = nrow(g),
      group_size_tier = unique(g$Group_Size_Tier)[[1]],
      group_type = unique(g$Group_Type)[[1]],
      stringsAsFactors = FALSE
    )
  }))
}

component_metrics_from_directional_weights <- function(weights,
                                                       truth,
                                                       job,
                                                       attribution_type,
                                                       selection_freq = NULL) {
  vars <- truth$Variable
  if (is.null(weights) || !nrow(weights)) {
    weights <- data.frame(
      Variable = character(0), Direction = character(0), Weight = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  weights$Direction <- ifelse(tolower(weights$Direction) %in% c("positive", "pos"),
    "Positive",
    ifelse(tolower(weights$Direction) %in% c("negative", "neg"), "Negative", weights$Direction)
  )
  wide <- weights |>
    dplyr::filter(.data$Direction %in% c("Positive", "Negative")) |>
    dplyr::group_by(.data$Variable, .data$Direction) |>
    dplyr::summarise(weight = sum(.data$Weight, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .data$Direction, values_from = .data$weight, values_fill = 0)
  if (!"Positive" %in% names(wide)) wide$Positive <- numeric(nrow(wide))
  if (!"Negative" %in% names(wide)) wide$Negative <- numeric(nrow(wide))
  scored <- truth |>
    dplyr::select(
      Variable, Group, True_Beta, True_Direction, IsActive, IsActiveGroup,
      Group_Size, Group_Size_Tier, Effect_Tier, Group_Type
    ) |>
    dplyr::left_join(wide, by = "Variable") |>
    dplyr::mutate(
      Positive = tidyr::replace_na(.data$Positive, 0),
      Negative = tidyr::replace_na(.data$Negative, 0),
      raw_attribution = pmax(abs(.data$Positive), abs(.data$Negative), na.rm = TRUE),
      estimated_direction = dplyr::case_when(
        .data$Positive == 0 & .data$Negative == 0 ~ "None",
        .data$Positive >= .data$Negative ~ "Positive",
        TRUE ~ "Negative"
      )
    )
  total <- sum(scored$raw_attribution, na.rm = TRUE)
  scored$normalized_attribution <- if (is.finite(total) && total > 0) scored$raw_attribution / total else NA_real_
  if (!is.null(selection_freq) && nrow(selection_freq)) {
    sf <- selection_freq |>
      dplyr::group_by(.data$Variable) |>
      dplyr::summarise(selection_frequency = max(.data$selection_frequency, na.rm = TRUE), .groups = "drop")
    scored <- scored |> dplyr::left_join(sf, by = "Variable")
  } else {
    scored$selection_frequency <- NA_real_
  }
  scored |>
    dplyr::mutate(
      scenario_id = job$scenario_id,
      battery = job$battery,
      family = job$family,
      method = job$method,
      data_seed = as.integer(job$data_seed),
      split_replicate = as.integer(job$split_replicate %||% NA_integer_),
      n = as.integer(job$n),
      p = as.integer(job$p),
      group_structure = job$group_structure,
      effect_profile = job$effect_profile,
      signal_profile = job$signal_profile,
      attribution_type = attribution_type,
      direction_correct = .data$IsActive & .data$estimated_direction == .data$True_Direction,
      sign_assignment_correct = dplyr::case_when(
        .data$IsActive ~ .data$estimated_direction == .data$True_Direction,
        !.data$IsActive ~ .data$estimated_direction == "None",
        TRUE ~ NA
      )
    ) |>
    dplyr::rename(
      variable = Variable,
      group = Group,
      true_beta = True_Beta,
      true_direction = True_Direction,
      is_active = IsActive,
      is_active_group = IsActiveGroup,
      group_size = Group_Size,
      group_size_tier = Group_Size_Tier,
      effect_tier = Effect_Tier,
      group_type = Group_Type,
      positive_weight = Positive,
      negative_weight = Negative
    )
}

component_metrics_from_groupwqs <- function(weights, group_coefs, truth, groups, job) {
  if (is.null(weights) || !nrow(weights)) {
    weights <- data.frame(Variable = truth$Variable, Weight = 0, stringsAsFactors = FALSE)
  }
  weight_by_var <- weights |>
    dplyr::group_by(.data$Variable) |>
    dplyr::summarise(raw_attribution = max(abs(.data$Weight), na.rm = TRUE), .groups = "drop")
  group_dir <- setNames(rep("None", length(groups)), names(groups))
  if (length(group_coefs)) {
    n_coef <- min(length(group_coefs), length(groups))
    group_dir[names(groups)[seq_len(n_coef)]] <- ifelse(group_coefs[seq_len(n_coef)] >= 0, "Positive", "Negative")
  }
  scored <- truth |>
    dplyr::select(
      Variable, Group, True_Beta, True_Direction, IsActive, IsActiveGroup,
      Group_Size, Group_Size_Tier, Effect_Tier, Group_Type
    ) |>
    dplyr::left_join(weight_by_var, by = "Variable") |>
    dplyr::mutate(
      raw_attribution = tidyr::replace_na(.data$raw_attribution, 0),
      estimated_direction = unname(group_dir[.data$Group]),
      positive_weight = ifelse(.data$estimated_direction == "Positive", .data$raw_attribution, 0),
      negative_weight = ifelse(.data$estimated_direction == "Negative", .data$raw_attribution, 0)
    )
  total <- sum(scored$raw_attribution, na.rm = TRUE)
  scored$normalized_attribution <- if (is.finite(total) && total > 0) scored$raw_attribution / total else NA_real_
  scored$selection_frequency <- NA_real_
  scored |>
    dplyr::mutate(
      scenario_id = job$scenario_id,
      battery = job$battery,
      family = job$family,
      method = job$method,
      data_seed = as.integer(job$data_seed),
      split_replicate = as.integer(job$split_replicate %||% NA_integer_),
      n = as.integer(job$n),
      p = as.integer(job$p),
      group_structure = job$group_structure,
      effect_profile = job$effect_profile,
      signal_profile = job$signal_profile,
      attribution_type = "WQS-type constrained weights",
      direction_correct = .data$IsActive & .data$estimated_direction == .data$True_Direction,
      sign_assignment_correct = dplyr::case_when(
        .data$IsActive ~ .data$estimated_direction == .data$True_Direction,
        !.data$IsActive ~ .data$estimated_direction == "None",
        TRUE ~ NA
      )
    ) |>
    dplyr::rename(
      variable = Variable,
      group = Group,
      true_beta = True_Beta,
      true_direction = True_Direction,
      is_active = IsActive,
      is_active_group = IsActiveGroup,
      group_size = Group_Size,
      group_size_tier = Group_Size_Tier,
      effect_tier = Effect_Tier,
      group_type = Group_Type
    )
}

summarize_component_metrics <- function(component_metrics) {
  if (is.null(component_metrics) || !nrow(component_metrics)) {
    return(data.frame(
      active_direction_accuracy = NA_real_,
      sign_assignment_accuracy = NA_real_,
      active_attribution = NA_real_,
      null_attribution = NA_real_,
      active_selection_frequency = NA_real_,
      null_selection_frequency = NA_real_
    ))
  }
  data.frame(
    active_direction_accuracy = mean(component_metrics$direction_correct[component_metrics$is_active], na.rm = TRUE),
    sign_assignment_accuracy = mean(component_metrics$sign_assignment_correct, na.rm = TRUE),
    active_attribution = sum(component_metrics$normalized_attribution[component_metrics$is_active], na.rm = TRUE),
    null_attribution = sum(component_metrics$normalized_attribution[!component_metrics$is_active], na.rm = TRUE),
    active_selection_frequency = mean(component_metrics$selection_frequency[component_metrics$is_active], na.rm = TRUE),
    null_selection_frequency = mean(component_metrics$selection_frequency[!component_metrics$is_active], na.rm = TRUE)
  )
}

empty_component_result <- function(truth, job, attribution_type = NA_character_) {
  truth |>
    dplyr::transmute(
      scenario_id = job$scenario_id,
      battery = job$battery,
      family = job$family,
      method = job$method,
      data_seed = as.integer(job$data_seed),
      split_replicate = as.integer(job$split_replicate %||% NA_integer_),
      n = as.integer(job$n),
      p = as.integer(job$p),
      group_structure = job$group_structure,
      effect_profile = job$effect_profile,
      signal_profile = job$signal_profile,
      attribution_type = attribution_type,
      variable = .data$Variable,
      group = .data$Group,
      true_beta = .data$True_Beta,
      true_direction = .data$True_Direction,
      is_active = .data$IsActive,
      is_active_group = .data$IsActiveGroup,
      group_size = .data$Group_Size,
      group_size_tier = .data$Group_Size_Tier,
      effect_tier = .data$Effect_Tier,
      group_type = .data$Group_Type,
      positive_weight = NA_real_,
      negative_weight = NA_real_,
      raw_attribution = NA_real_,
      normalized_attribution = NA_real_,
      estimated_direction = NA_character_,
      selection_frequency = NA_real_,
      direction_correct = NA,
      sign_assignment_correct = NA
    )
}

extract_sglwqs_selection_frequency <- function(fit) {
  if (is.null(fit$boot_info)) return(NULL)
  pos <- fit$boot_info$selection_freq_pos %||% NULL
  neg <- fit$boot_info$selection_freq_neg %||% NULL
  out <- list()
  if (!is.null(pos)) {
    out[[length(out) + 1L]] <- data.frame(
      Variable = names(pos), direction = "Positive",
      selection_frequency = as.numeric(pos),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(neg)) {
    out[[length(out) + 1L]] <- data.frame(
      Variable = names(neg), direction = "Negative",
      selection_frequency = as.numeric(neg),
      stringsAsFactors = FALSE
    )
  }
  if (length(out)) dplyr::bind_rows(out) else NULL
}

extract_sglwqs_direction_results <- function(fit, truth, job, settings, elapsed, warnings_seen) {
  rows <- truth_group_direction_rows(truth, job)
  rows$retained <- FALSE
  rows$excluded_by_minor_threshold <- FALSE
  rows$estimate <- NA_real_
  rows$std_error <- NA_real_
  rows$p_value <- NA_real_
  rows$rejected_0_05 <- FALSE
  rows$sign <- NA_character_
  rows$pos_mass <- NA_real_
  rows$neg_mass <- NA_real_
  rows$direction_mass <- NA_real_
  rows$runtime_sec <- elapsed
  rows$n_boot_requested <- as.integer(settings$n_boot)
  boot_success <- fit$boot_info$boot_success %||% logical(0)
  rows$n_boot_success <- if (length(boot_success)) sum(boot_success, na.rm = TRUE) else NA_integer_
  rows$boot_success_rate <- if (length(boot_success)) rows$n_boot_success[[1]] / length(boot_success) else NA_real_
  rows$warning_summary <- paste(unique(warnings_seen), collapse = " | ")

  gr <- fit$validation_info$group_results %||% list()
  for (i in seq_len(nrow(rows))) {
    group <- rows$group[[i]]
    direction <- rows$direction[[i]]
    prefix <- if (identical(direction, "positive")) "pos" else "neg"
    rows$pos_mass[i] <- as.numeric((fit$pos_index_sum_by_group %||% list())[[group]] %||% NA_real_)
    rows$neg_mass[i] <- as.numeric((fit$neg_index_sum_by_group %||% list())[[group]] %||% NA_real_)
    rows$direction_mass[i] <- if (identical(direction, "positive")) rows$pos_mass[i] else rows$neg_mass[i]
    if (!is.null(gr[[group]])) {
      g <- gr[[group]]
      rows$excluded_by_minor_threshold[i] <- isTRUE(g[[paste0(prefix, "_excluded")]] %||% FALSE)
      rows$estimate[i] <- finite_or_na(g[[paste0(prefix, "_estimate")]] %||% NA_real_)
      rows$std_error[i] <- finite_or_na(g[[paste0(prefix, "_se")]] %||% NA_real_)
      rows$p_value[i] <- finite_or_na(g[[paste0(prefix, "_pvalue")]] %||% NA_real_)
    }
  }
  rows$retained <- is.finite(rows$p_value) & !rows$excluded_by_minor_threshold
  rows$rejected_0_05 <- rows$retained & is.finite(rows$p_value) & rows$p_value < settings$alpha
  rows$sign <- ifelse(is.finite(rows$estimate) & rows$estimate > 0, "positive",
    ifelse(is.finite(rows$estimate) & rows$estimate < 0, "negative", NA_character_)
  )
  rows
}

empty_sglwqs_direction_results <- function(truth, job, settings, elapsed = NA_real_, warnings_seen = character(0)) {
  rows <- truth_group_direction_rows(truth, job)
  rows$retained <- FALSE
  rows$excluded_by_minor_threshold <- NA
  rows$estimate <- NA_real_
  rows$std_error <- NA_real_
  rows$p_value <- NA_real_
  rows$rejected_0_05 <- FALSE
  rows$sign <- NA_character_
  rows$pos_mass <- NA_real_
  rows$neg_mass <- NA_real_
  rows$direction_mass <- NA_real_
  rows$runtime_sec <- elapsed
  rows$n_boot_requested <- as.integer(settings$n_boot)
  rows$n_boot_success <- NA_integer_
  rows$boot_success_rate <- NA_real_
  rows$warning_summary <- paste(unique(warnings_seen), collapse = " | ")
  rows
}

method_metric_row <- function(job,
                              fit_success,
                              elapsed,
                              warnings_seen,
                              error = NULL,
                              component_summary = NULL,
                              bootstrap_attempted = NA_integer_,
                              bootstrap_succeeded = NA_integer_,
                              selected_lambda = NA_real_,
                              retained_group_directions = NA_integer_,
                              failure_stage = "fit") {
  component_summary <- component_summary %||% summarize_component_metrics(NULL)
  data.frame(
    scenario_id = job$scenario_id,
    battery = job$battery,
    family = job$family,
    method = job$method,
    data_seed = as.integer(job$data_seed),
    split_replicate = as.integer(job$split_replicate %||% NA_integer_),
    n = as.integer(job$n),
    p = as.integer(job$p),
    group_structure = job$group_structure,
    effect_profile = job$effect_profile,
    signal_profile = job$signal_profile,
    fit_success = isTRUE(fit_success),
    failure_stage = if (isTRUE(fit_success)) NA_character_ else failure_stage,
    error_class = if (isTRUE(fit_success)) NA_character_ else classify_error(conditionMessage(error)),
    error_message = if (isTRUE(fit_success)) NA_character_ else conditionMessage(error),
    warning_summary = paste(unique(warnings_seen), collapse = " | "),
    runtime_sec = elapsed,
    bootstrap_attempted = as.integer(bootstrap_attempted),
    bootstrap_succeeded = as.integer(bootstrap_succeeded),
    bootstrap_success_rate = ifelse(is.finite(bootstrap_attempted) && bootstrap_attempted > 0,
      bootstrap_succeeded / bootstrap_attempted, NA_real_
    ),
    selected_lambda = selected_lambda,
    retained_group_directions = as.integer(retained_group_directions),
    active_direction_accuracy = component_summary$active_direction_accuracy,
    sign_assignment_accuracy = component_summary$sign_assignment_accuracy,
    active_attribution = component_summary$active_attribution,
    null_attribution = component_summary$null_attribution,
    active_selection_frequency = component_summary$active_selection_frequency,
    null_selection_frequency = component_summary$null_selection_frequency,
    stringsAsFactors = FALSE
  )
}

fit_qgcomp_job <- function(dat, truth, groups, job, settings) {
  vars <- unlist(groups, use.names = FALSE)
  formula <- stats::as.formula(paste("Y ~", paste(c(vars, legacy13_covariates()), collapse = " + ")))
  res <- capture_fit({
    fit <- qgcomp::qgcomp.glm.noboot(
      formula,
      expnms = vars,
      data = dat,
      family = if (identical(job$family, "binomial")) stats::binomial() else stats::gaussian(),
      q = settings$quantiles
    )
    weights <- data.frame(
      Variable = c(names(fit$pos.weights), names(fit$neg.weights)),
      Direction = c(rep("Positive", length(fit$pos.weights)), rep("Negative", length(fit$neg.weights))),
      Weight = c(fit$pos.weights, fit$neg.weights),
      stringsAsFactors = FALSE
    )
    list(fit = fit, weights = weights)
  })
  if (res$ok) {
    comp <- component_metrics_from_directional_weights(
      res$value$weights, truth, job,
      attribution_type = "coefficient-derived attribution"
    )
    metrics <- method_metric_row(job, TRUE, res$elapsed, res$warnings,
      component_summary = summarize_component_metrics(comp),
      bootstrap_attempted = 0L, bootstrap_succeeded = 0L
    )
  } else {
    comp <- empty_component_result(truth, job, "coefficient-derived attribution")
    metrics <- method_metric_row(job, FALSE, res$elapsed, res$warnings, error = res$error)
  }
  list(method_metrics = metrics, component_metrics = comp, sglwqs_direction_results = data.frame())
}

fit_gwqs_job <- function(dat, truth, groups, job, settings, split) {
  vars <- unlist(groups, use.names = FALSE)
  res <- capture_fit({
    validation_rows <- list(seq_len(nrow(dat)) %in% split$validation)
    fit <- gWQS::gwqs(
      stats::as.formula("Y ~ pwqs + nwqs + sex + age + bmi"),
      mix_name = vars,
      data = dat,
      q = settings$quantiles,
      validation = 1 - settings$train_prop,
      validation_rows = validation_rows,
      b = settings$n_boot,
      b1_pos = TRUE,
      rh = 1,
      family = job$family,
      seed = as.integer(job$method_seed),
      plan_strategy = "sequential"
    )
    fw <- fit$final_weights
    weights <- dplyr::bind_rows(
      data.frame(
        Variable = fw$mix_name,
        Direction = "Positive",
        Weight = fw$mean_weight_p,
        stringsAsFactors = FALSE
      ),
      data.frame(
        Variable = fw$mix_name,
        Direction = "Negative",
        Weight = fw$mean_weight_n,
        stringsAsFactors = FALSE
      )
    )
    list(fit = fit, weights = weights)
  })
  if (res$ok) {
    comp <- component_metrics_from_directional_weights(
      res$value$weights, truth, job,
      attribution_type = "WQS-type constrained weights"
    )
    metrics <- method_metric_row(job, TRUE, res$elapsed, res$warnings,
      component_summary = summarize_component_metrics(comp),
      bootstrap_attempted = settings$n_boot, bootstrap_succeeded = settings$n_boot
    )
  } else {
    comp <- empty_component_result(truth, job, "WQS-type constrained weights")
    metrics <- method_metric_row(job, FALSE, res$elapsed, res$warnings, error = res$error,
      bootstrap_attempted = settings$n_boot, bootstrap_succeeded = NA_integer_
    )
  }
  list(method_metrics = metrics, component_metrics = comp, sglwqs_direction_results = data.frame())
}

fit_groupwqs_job <- function(dat, truth, groups, job, settings, split) {
  covars <- legacy13_covariates()
  group_list <- unname(groups)
  res <- capture_fit({
    set.seed(as.integer(job$method_seed))
    x_gw <- groupWQS::make.X(dat, num.groups = length(group_list), groups = group_list)
    x_s <- groupWQS::make.x.s(dat, num.groups = length(group_list), groups = group_list)
    fit <- groupWQS::gwqs.fit(
      y = dat$Y[split$validation],
      y.train = dat$Y[split$train],
      x = x_gw[split$validation, , drop = FALSE],
      x.train = x_gw[split$train, , drop = FALSE],
      z = dat[split$validation, covars, drop = FALSE],
      z.train = dat[split$train, covars, drop = FALSE],
      x.s = x_s,
      B = settings$n_boot,
      n.quantiles = settings$quantiles,
      func = if (identical(job$family, "binomial")) "binary" else "continuous"
    )
    weights <- purrr::imap_dfr(fit$weights, function(w, i) {
      data.frame(
        Variable = names(w),
        Weight = as.numeric(w),
        Direction = "Combined",
        stringsAsFactors = FALSE
      )
    })
    coefs <- stats::coef(fit$fit)[grep("^GWQS", names(stats::coef(fit$fit)))]
    list(fit = fit, weights = weights, group_coefs = coefs)
  })
  if (res$ok) {
    comp <- component_metrics_from_groupwqs(res$value$weights, res$value$group_coefs, truth, groups, job)
    metrics <- method_metric_row(job, TRUE, res$elapsed, res$warnings,
      component_summary = summarize_component_metrics(comp),
      bootstrap_attempted = settings$n_boot, bootstrap_succeeded = settings$n_boot
    )
  } else {
    comp <- empty_component_result(truth, job, "WQS-type constrained weights")
    metrics <- method_metric_row(job, FALSE, res$elapsed, res$warnings, error = res$error,
      bootstrap_attempted = settings$n_boot, bootstrap_succeeded = NA_integer_
    )
  }
  list(method_metrics = metrics, component_metrics = comp, sglwqs_direction_results = data.frame())
}

sglwqs_call <- function(dat, groups, job, settings, checkpoint_dir) {
  vars <- unlist(groups, use.names = FALSE)
  args <- list(
    X = dat[, vars, drop = FALSE],
    y = dat$Y,
    covariates = dat[, legacy13_covariates(), drop = FALSE],
    groups = groups,
    family = job$family,
    n_quantiles = settings$quantiles,
    lambda = settings$lambda,
    nfolds = settings$nfolds,
    bootstrap = TRUE,
    n_boot = settings$n_boot,
    parallel = FALSE,
    keep_boot_matrices = FALSE,
    checkpoint_dir = checkpoint_dir,
    checkpoint_interval = 50L,
    cleanup_checkpoint = TRUE,
    validation = TRUE,
    train_prop = settings$train_prop,
    # This seed reproduces the same split represented by `job$split_seed`.
    # Bootstrap draws then continue deterministically from that RNG stream.
    seed = as.integer(job$split_seed),
    minor_threshold = settings$minor_threshold,
    asparse = settings$asparse,
    verbose = FALSE
  )
  tryCatch(
    do.call(sglwqs::sglwqs, args),
    error = function(e) {
      if (grepl("unused argument.*asparse", conditionMessage(e))) {
        args$asparse <- NULL
        do.call(sglwqs::sglwqs, args)
      } else {
        stop(e)
      }
    }
  )
}

fit_sglwqs_job <- function(dat, truth, groups, job, settings) {
  r2r2_load_local_sglwqs()
  checkpoint_dir <- r2r2_result_file(
    "raw", safe_filename(job$battery), safe_filename(job$scenario_id),
    "SGL-WQS", paste0("checkpoint_", as.integer(job$data_seed), "_", as.integer(job$split_replicate %||% 0L))
  )
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  res <- capture_fit({
    fit <- sglwqs_call(dat, groups, job, settings, checkpoint_dir)
    weights <- sglwqs::extract_weights(fit, direction = "both") |>
      dplyr::transmute(
        Variable = .data$variable,
        Direction = tools::toTitleCase(.data$direction),
        Weight = .data$weight
      )
    selection <- extract_sglwqs_selection_frequency(fit)
    list(fit = fit, weights = weights, selection = selection)
  })
  if (res$ok) {
    comp <- component_metrics_from_directional_weights(
      res$value$weights, truth, job,
      attribution_type = "WQS-type constrained weights",
      selection_freq = res$value$selection
    )
    dir_rows <- extract_sglwqs_direction_results(res$value$fit, truth, job, settings, res$elapsed, res$warnings)
    sel <- res$value$fit$selection_diagnostics %||% list()
    boot_success <- res$value$fit$boot_info$boot_success %||% logical(0)
    metrics <- method_metric_row(job, TRUE, res$elapsed, res$warnings,
      component_summary = summarize_component_metrics(comp),
      bootstrap_attempted = if (length(boot_success)) length(boot_success) else settings$n_boot,
      bootstrap_succeeded = if (length(boot_success)) sum(boot_success, na.rm = TRUE) else NA_integer_,
      selected_lambda = finite_or_na(sel$selected_lambda %||% NA_real_),
      retained_group_directions = sum(dir_rows$retained, na.rm = TRUE)
    )
  } else {
    comp <- empty_component_result(truth, job, "WQS-type constrained weights")
    dir_rows <- empty_sglwqs_direction_results(truth, job, settings, res$elapsed, res$warnings)
    metrics <- method_metric_row(job, FALSE, res$elapsed, res$warnings, error = res$error,
      bootstrap_attempted = settings$n_boot, bootstrap_succeeded = NA_integer_
    )
  }
  list(method_metrics = metrics, component_metrics = comp, sglwqs_direction_results = dir_rows)
}

run_method_job <- function(job, settings = r2r2_settings()) {
  r2r2_set_thread_env()
  r2r2_load_packages()
  if (identical(job$method, "SGL-WQS")) {
    r2r2_load_local_sglwqs()
  }
  dat <- generate_job_data(job)
  qa_generated_data(dat, job)
  groups <- attr(dat, "groups")
  truth <- attr(dat, "truth")
  settings_local <- settings
  paired_train_prop <- attr(dat, "paired_train_prop") %||% NA_real_
  if (is.finite(paired_train_prop)) {
    settings_local$train_prop <- paired_train_prop
  }
  split <- make_train_validation_split(
    nrow(dat),
    as.integer(job$split_seed),
    train_prop = settings_local$train_prop,
    y = dat$Y,
    family = job$family
  )
  assert_or_stop(length(intersect(split$train, split$validation)) == 0L, "Train/validation split overlap.")
  assert_or_stop(length(union(split$train, split$validation)) == nrow(dat), "Train/validation split is not exhaustive.")

  out <- switch(job$method,
    qgcomp = fit_qgcomp_job(dat, truth, groups, job, settings_local),
    gWQS = fit_gwqs_job(dat, truth, groups, job, settings_local, split),
    groupWQS = fit_groupwqs_job(dat, truth, groups, job, settings_local, split),
    `SGL-WQS` = fit_sglwqs_job(dat, truth, groups, job, settings_local),
    stop("Unknown method: ", job$method, call. = FALSE)
  )
  out$schema_version <- "reviewer2_round2_v1"
  out$job <- as.data.frame(job, stringsAsFactors = FALSE)
  out$truth <- truth
  out$data_diagnostics <- data.frame(
    data_hash = digest::digest(dat, algo = "sha256"),
    split_hash = digest::digest(split, algo = "sha256"),
    train_n = length(split$train),
    validation_n = length(split$validation),
    y_mean = mean(dat$Y),
    y_sd = stats::sd(dat$Y),
    y_min = min(dat$Y),
    y_max = max(dat$Y),
    y_positive_rate = if (identical(job$family, "binomial")) mean(dat$Y == 1) else NA_real_,
    y_positive_count = if (identical(job$family, "binomial")) sum(dat$Y == 1) else NA_integer_,
    binomial_intercept = attr(dat, "binomial_intercept") %||% NA_real_,
    binomial_target_prevalence = attr(dat, "binomial_target_prevalence") %||% NA_real_,
    binomial_expected_prevalence = attr(dat, "binomial_expected_prevalence") %||% NA_real_,
    gaussian_residual_sd = attr(dat, "gaussian_residual_sd") %||% NA_real_,
    within_rho = attr(dat, "within_rho") %||% attr(dat, "within_r") %||% NA_real_,
    cross_rho = attr(dat, "cross_rho") %||% attr(dat, "between_r") %||% NA_real_,
    effective_train_prop = settings_local$train_prop,
    paired_binary_prevalence = attr(dat, "paired_binary_prevalence") %||% NA_real_,
    paired_exposure_covariate_hash = attr(dat, "paired_exposure_covariate_hash") %||% NA_character_,
    paired_reference_train_hash = if (!is.null(attr(dat, "paired_reference_train_idx"))) {
      digest::digest(sort(attr(dat, "paired_reference_train_idx")), algo = "sha256")
    } else NA_character_,
    paired_reference_validation_hash = if (!is.null(attr(dat, "paired_reference_validation_idx"))) {
      digest::digest(sort(attr(dat, "paired_reference_validation_idx")), algo = "sha256")
    } else NA_character_,
    stringsAsFactors = FALSE
  )
  out$seed_metadata <- data.frame(
    data_seed = as.integer(job$data_seed),
    split_seed = as.integer(job$split_seed),
    method_seed = as.integer(job$method_seed),
    bootstrap_seed = as.integer(job$bootstrap_seed %||% NA_integer_),
    stringsAsFactors = FALSE
  )
  out
}

run_method_job_resumable <- function(job, settings = r2r2_settings(), force = FALSE) {
  path <- job_output_path(job)
  if (!force && is_complete_job_file(path)) {
    return(readRDS(path))
  }
  result <- tryCatch(
    run_method_job(job, settings),
    error = function(e) {
      dat <- tryCatch(generate_job_data(job), error = function(e2) NULL)
      truth <- if (is.null(dat)) {
        data.frame()
      } else {
        attr(dat, "truth")
      }
      comp <- if (nrow(truth)) {
        empty_component_result(truth, job, NA_character_)
      } else {
        data.frame()
      }
      settings_local <- settings
      dir_rows <- if (identical(job$method, "SGL-WQS") && nrow(truth)) {
        empty_sglwqs_direction_results(truth, job, settings_local)
      } else {
        data.frame()
      }
      list(
        schema_version = "reviewer2_round2_v1",
        job = as.data.frame(job, stringsAsFactors = FALSE),
        method_metrics = method_metric_row(
          job, FALSE, NA_real_, character(0), error = e,
          failure_stage = if (grepl("Unknown method", conditionMessage(e), fixed = TRUE)) "dispatch" else "setup_or_generation"
        ),
        component_metrics = comp,
        sglwqs_direction_results = dir_rows,
        truth = truth,
        data_diagnostics = data.frame(),
        seed_metadata = data.frame(
          data_seed = as.integer(job$data_seed),
          split_seed = as.integer(job$split_seed),
          method_seed = as.integer(job$method_seed),
          bootstrap_seed = as.integer(job$bootstrap_seed %||% NA_integer_)
        )
      )
    }
  )
  atomic_save_rds(result, path)
  result
}
