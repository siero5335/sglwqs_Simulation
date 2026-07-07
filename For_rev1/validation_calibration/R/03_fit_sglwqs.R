load_sglwqs_source <- function() {
  activate_local_lib()
  source_path <- sglwqs_source_path()
  if (!is.na(source_path) && nzchar(source_path)) {
    if (!requireNamespace("pkgload", quietly = TRUE)) {
      stop("Package 'pkgload' is required when SGLWQS_SOURCE is set.", call. = FALSE)
    }
    pkgload::load_all(source_path, quiet = TRUE, export_all = FALSE)
    return(invisible(TRUE))
  }
  if (!requireNamespace("sglwqs", quietly = TRUE)) {
    stop("Install the 'sglwqs' package or set SGLWQS_SOURCE to a source tree.", call. = FALSE)
  }
  invisible(TRUE)
}

classify_backend_message <- function(x) {
  x <- paste(x, collapse = " | ")
  if (!nzchar(x)) {
    return(NA_character_)
  }
  if (grepl("lambda|path|boundary", x, ignore.case = TRUE)) {
    return("lambda_path")
  }
  if (grepl("conver|jerr|maxit|irls", x, ignore.case = TRUE)) {
    return("backend_convergence")
  }
  if (grepl("fold|cross-validation|cv", x, ignore.case = TRUE)) {
    return("cross_validation")
  }
  if (grepl("rank|singular|alias", x, ignore.case = TRUE)) {
    return("rank_deficiency")
  }
  "other"
}

run_sglwqs_job <- function(job, profile) {
  load_sglwqs_source()
  ensure_vc_dirs()

  n_boot <- if (identical(job$job_type, "split_stability")) {
    as.integer(profile$split_n_boot)
  } else {
    as.integer(profile$n_boot)
  }

  dat <- generate_validation_data(
    n = as.integer(job$n),
    seed = as.integer(job$data_seed),
    effect = job$effect,
    outcome_family = job$outcome_family,
    within_r = as.numeric(job$within_r),
    between_r = as.numeric(job$between_r)
  )

  checkpoint_dir <- vc_file(
    "output", "checkpoints", profile$name, job$job_type, job$job_id
  )
  dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

  warnings_seen <- character(0)
  start_time <- Sys.time()
  elapsed <- NA_real_

  fit <- tryCatch(
    withCallingHandlers(
      sglwqs::sglwqs(
        X = dat[, validation_exposures(), drop = FALSE],
        y = dat$Y,
        covariates = dat[, validation_covariates(), drop = FALSE],
        groups = validation_groups(),
        family = job$outcome_family,
        n_quantiles = as.integer(profile$n_quantiles),
        validation = TRUE,
        train_prop = as.numeric(profile$train_prop),
        bootstrap = isTRUE(profile$bootstrap),
        n_boot = n_boot,
        parallel = isTRUE(profile$parallel_inside_fit),
        nfolds = as.integer(profile$nfolds),
        lambda = profile$lambda,
        minor_threshold = as.numeric(profile$minor_threshold),
        asparse = as.numeric(profile$asparse),
        nlambda = as.integer(profile$nlambda),
        seed = as.integer(job$fit_seed),
        checkpoint_dir = checkpoint_dir,
        checkpoint_interval = as.integer(profile$checkpoint_interval),
        cleanup_checkpoint = isTRUE(profile$cleanup_checkpoint),
        verbose = FALSE
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (inherits(fit, "error")) {
    result <- make_failure_result(job, profile, fit, warnings_seen, elapsed, n_boot)
  } else {
    result <- extract_sglwqs_result(fit, job, profile, warnings_seen, elapsed, n_boot)
  }

  out_path <- vc_file(
    "output", "manifests", profile$name, job$job_type,
    paste0(job$job_id, ".rds")
  )
  atomic_save_rds(result, out_path)
  result
}

make_failure_result <- function(job, profile, err, warnings_seen, elapsed, n_boot) {
  metrics <- empty_direction_rows(job)
  metrics$profile <- profile$name
  metrics$fit_success <- FALSE
  metrics$error_message <- conditionMessage(err)
  metrics$error_class <- classify_backend_message(conditionMessage(err))
  metrics$warning_messages <- paste(unique(warnings_seen), collapse = " | ")
  metrics$runtime_sec <- elapsed
  metrics$n_boot_requested <- n_boot

  list(
    job = as.data.frame(job, stringsAsFactors = FALSE),
    direction_results = metrics,
    diagnostics = data.frame(
      job_id = job$job_id,
      profile = profile$name,
      job_type = job$job_type,
      scenario_id = job$scenario_id,
      scenario_family = job$scenario_family,
      n = as.integer(job$n),
      outcome_family = job$outcome_family,
      fit_success = FALSE,
      runtime_sec = elapsed,
      n_boot_requested = n_boot,
      n_boot_success = NA_integer_,
      boot_success_rate = NA_real_,
      error_class = classify_backend_message(conditionMessage(err)),
      error_message = conditionMessage(err),
      warning_messages = paste(unique(warnings_seen), collapse = " | "),
      stringsAsFactors = FALSE
    )
  )
}

empty_direction_rows <- function(job) {
  td <- truth_group_direction(job$effect)
  data.frame(
    profile = NA_character_,
    job_id = job$job_id,
    job_type = job$job_type,
    scenario_id = job$scenario_id,
    scenario_family = job$scenario_family,
    scenario_label = job$scenario_label,
    replicate = as.integer(job$replicate),
    n = as.integer(job$n),
    outcome_family = job$outcome_family,
    effect = job$effect,
    within_r = as.numeric(job$within_r),
    between_r = as.numeric(job$between_r),
    group = td$group,
    direction = td$direction,
    is_true_null = td$is_true_null,
    active_variable_count = td$active_variable_count,
    true_beta_sum = td$true_beta_sum,
    retained = FALSE,
    excluded_by_minor_threshold = FALSE,
    estimate = NA_real_,
    std_error = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    p_value = NA_real_,
    rejected_0_05 = FALSE,
    sign = NA_character_,
    pos_mass = NA_real_,
    neg_mass = NA_real_,
    direction_mass = NA_real_,
    fit_success = FALSE,
    runtime_sec = NA_real_,
    n_boot_requested = NA_integer_,
    n_boot_success = NA_integer_,
    boot_success_rate = NA_real_,
    lambda_path_source = NA_character_,
    lambda_path_length = NA_integer_,
    lambda_path_min = NA_real_,
    lambda_path_max = NA_real_,
    selected_lambda = NA_real_,
    selected_lambda_at_path_boundary = NA,
    all_zero_exposure = NA,
    backend_jerr = NA_integer_,
    error_class = NA_character_,
    error_message = NA_character_,
    warning_messages = NA_character_,
    stringsAsFactors = FALSE
  )
}

extract_sglwqs_result <- function(fit, job, profile, warnings_seen, elapsed, n_boot) {
  rows <- empty_direction_rows(job)
  rows$profile <- profile$name
  rows$fit_success <- TRUE
  rows$runtime_sec <- elapsed
  rows$n_boot_requested <- n_boot

  boot_success <- fit$boot_info$boot_success %||% logical(0)
  n_boot_success <- if (length(boot_success)) sum(isTRUE(boot_success) | boot_success, na.rm = TRUE) else NA_integer_
  boot_success_rate <- if (length(boot_success)) n_boot_success / length(boot_success) else NA_real_

  rows$n_boot_success <- n_boot_success
  rows$boot_success_rate <- boot_success_rate
  rows$warning_messages <- paste(unique(warnings_seen), collapse = " | ")
  rows$error_class <- classify_backend_message(warnings_seen)

  sel <- fit$selection_diagnostics %||% list()
  rows$lambda_path_source <- as.character(sel$lambda_path_source %||% NA_character_)
  rows$lambda_path_length <- as.integer(sel$lambda_path_length %||% NA_integer_)
  rows$lambda_path_min <- as.numeric(sel$lambda_path_min %||% NA_real_)
  rows$lambda_path_max <- as.numeric(sel$lambda_path_max %||% NA_real_)
  rows$selected_lambda <- as.numeric(sel$selected_lambda %||% NA_real_)
  rows$selected_lambda_at_path_boundary <- isTRUE(sel$selected_lambda_at_path_boundary)
  rows$all_zero_exposure <- isTRUE(sel$all_zero_exposure)
  rows$backend_jerr <- as.integer(sel$backend_jerr %||% NA_integer_)

  excluded <- fit$excluded_directions %||% list()
  val <- fit$validation_info %||% list()
  gr <- val$group_results %||% list()

  for (i in seq_len(nrow(rows))) {
    group <- rows$group[i]
    direction <- rows$direction[i]
    rows$pos_mass[i] <- as.numeric((fit$pos_index_sum_by_group %||% list())[[group]] %||% NA_real_)
    rows$neg_mass[i] <- as.numeric((fit$neg_index_sum_by_group %||% list())[[group]] %||% NA_real_)
    rows$direction_mass[i] <- if (identical(direction, "positive")) rows$pos_mass[i] else rows$neg_mass[i]
    rows$excluded_by_minor_threshold[i] <- identical(tolower(excluded[[group]] %||% ""), direction)
    rows$retained[i] <- !rows$excluded_by_minor_threshold[i]

    if (!is.null(gr[[group]])) {
      g <- gr[[group]]
      prefix <- if (identical(direction, "positive")) "pos" else "neg"
      rows$estimate[i] <- as.numeric(g[[paste0(prefix, "_estimate")]] %||% NA_real_)
      rows$std_error[i] <- as.numeric(g[[paste0(prefix, "_se")]] %||% NA_real_)
      rows$p_value[i] <- as.numeric(g[[paste0(prefix, "_pvalue")]] %||% NA_real_)
    }
  }

  rows$retained <- rows$retained & is.finite(rows$p_value)
  rows$conf_low <- rows$estimate - 1.96 * rows$std_error
  rows$conf_high <- rows$estimate + 1.96 * rows$std_error
  rows$rejected_0_05 <- is.finite(rows$p_value) & rows$p_value < as.numeric(profile$alpha)
  rows$sign <- ifelse(is.finite(rows$estimate) & rows$estimate > 0, "positive",
    ifelse(is.finite(rows$estimate) & rows$estimate < 0, "negative", NA_character_)
  )

  diagnostics <- data.frame(
    job_id = job$job_id,
    profile = profile$name,
    job_type = job$job_type,
    scenario_id = job$scenario_id,
    scenario_family = job$scenario_family,
    n = as.integer(job$n),
    outcome_family = job$outcome_family,
    fit_success = TRUE,
    runtime_sec = elapsed,
    n_boot_requested = n_boot,
    n_boot_success = n_boot_success,
    boot_success_rate = boot_success_rate,
    lambda_path_source = as.character(sel$lambda_path_source %||% NA_character_),
    lambda_path_length = as.integer(sel$lambda_path_length %||% NA_integer_),
    lambda_path_min = as.numeric(sel$lambda_path_min %||% NA_real_),
    lambda_path_max = as.numeric(sel$lambda_path_max %||% NA_real_),
    selected_lambda = as.numeric(sel$selected_lambda %||% NA_real_),
    selected_lambda_at_path_boundary = isTRUE(sel$selected_lambda_at_path_boundary),
    all_zero_exposure = isTRUE(sel$all_zero_exposure),
    backend_jerr = as.integer(sel$backend_jerr %||% NA_integer_),
    error_class = classify_backend_message(warnings_seen),
    error_message = NA_character_,
    warning_messages = paste(unique(warnings_seen), collapse = " | "),
    stringsAsFactors = FALSE
  )

  list(
    job = as.data.frame(job, stringsAsFactors = FALSE),
    direction_results = rows,
    diagnostics = diagnostics
  )
}
