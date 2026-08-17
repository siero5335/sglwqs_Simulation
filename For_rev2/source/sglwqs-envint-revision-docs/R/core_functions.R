#' Internal Function: Fit Single SGL Model
#'
#' Core fitting function used internally by sglwqs.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Processed covariate matrix (can be NULL).
#' @param var_names Exposure variable names.
#' @param cov_names Covariate names.
#' @param groups Group list.
#' @param group_by_compound Logical.
#' @param group_structure Character.
#' @param penalize_covariates Logical.
#' @param family Character.
#' @param lambda Character or numeric.
#' @param lambda_path Optional numeric lambda sequence passed to cv.sparsegl.
#' @param nfolds Integer.
#' @param obs_weights Optional numeric observation weights for selection loss.
#' @param ... Additional arguments to cv.sparsegl.
#'
#' @return A list with fitted model results.
#'
#' @keywords internal
fit_sgl_core <- function(X_quantile, y, cov_matrix, var_names, cov_names,
                          groups, group_by_compound, group_structure,
                          penalize_covariates, family, lambda, nfolds,
                          lambda_path = NULL,
                          obs_weights = NULL, ...) {
  
  n <- nrow(X_quantile)
  p <- ncol(X_quantile)
  q <- if (is.null(cov_matrix)) 0 else ncol(cov_matrix)
  obs_weights <- .validate_obs_weights(obs_weights, n)
  lambda_path <- .validate_lambda_path(lambda_path)
  sparsegl_family <- .coerce_sparsegl_family(family, obs_weights)
  
  # Create group mapping (for compound groups)
  if (!is.null(groups) && group_by_compound) {
    var_to_group <- rep(NA_integer_, p)
    for (i in seq_along(groups)) {
      var_idx <- which(var_names %in% groups[[i]])
      var_to_group[var_idx] <- i
    }
    unassigned <- which(is.na(var_to_group))
    if (length(unassigned) > 0) {
      next_group <- length(groups) + 1
      for (idx in unassigned) {
        var_to_group[idx] <- next_group
        next_group <- next_group + 1
      }
    }
    n_compound_groups <- max(var_to_group)
  } else {
    var_to_group <- NULL
    n_compound_groups <- NULL
  }
  
  # Build design matrix
  if (q > 0) {
    design_matrix <- cbind(X_quantile, -X_quantile, cov_matrix)
    
    if (!is.null(groups) && group_by_compound) {
      # group_structure = "direction": treat positive/negative as separate groups
      group_pos <- var_to_group
      group_neg <- var_to_group + n_compound_groups
      group_cov <- seq(2 * n_compound_groups + 1, 2 * n_compound_groups + q)
      group <- c(group_pos, group_neg, group_cov)
      
      n_total_groups <- 2 * n_compound_groups + q
      if (penalize_covariates) {
        pf_sparse <- rep(1, 2 * p + q)
        pf_group <- rep(1, n_total_groups)
        lower_bnd <- c(rep(0, 2 * n_compound_groups), rep(-Inf, q))
      } else {
        pf_sparse <- c(rep(1, 2 * p), rep(0, q))
        pf_group <- c(rep(1, 2 * n_compound_groups), rep(0, q))
        lower_bnd <- c(rep(0, 2 * n_compound_groups), rep(-Inf, q))
      }
    } else {
      # No user groups: G=1 direction-split (manuscript model)
      # All positive copies in group 1, all negative copies in group 2
      group <- c(rep(1L, p), rep(2L, p), seq(3L, 2L + q))
      if (penalize_covariates) {
        pf_sparse <- rep(1, 2 * p + q)
        pf_group <- rep(1, 2 + q)
        lower_bnd <- c(rep(0, 2), rep(-Inf, q))
      } else {
        pf_sparse <- c(rep(1, 2 * p), rep(0, q))
        pf_group <- c(rep(1, 2), rep(0, q))
        lower_bnd <- c(rep(0, 2), rep(-Inf, q))
      }
    }
  } else {
    design_matrix <- cbind(X_quantile, -X_quantile)

    if (!is.null(groups) && group_by_compound) {
      # group_structure = "direction": treat positive/negative as separate groups
      group_pos <- var_to_group
      group_neg <- var_to_group + n_compound_groups
      group <- c(group_pos, group_neg)
      pf_sparse <- rep(1, 2 * p)
      pf_group <- rep(1, 2 * n_compound_groups)
      lower_bnd <- rep(0, 2 * n_compound_groups)
    } else {
      # No user groups: G=1 direction-split (manuscript model)
      group <- c(rep(1L, p), rep(2L, p))
      pf_sparse <- rep(1, 2 * p)
      pf_group <- rep(1, 2)
      lower_bnd <- rep(0, 2)
    }
  }
  
  backend_warnings <- character(0)

  # Fit Sparse Group Lasso (with error handling)
  fit <- tryCatch({
    withCallingHandlers(
      sparsegl::cv.sparsegl(
        x = design_matrix,
        y = y,
        group = group,
        family = sparsegl_family,
        pf_sparse = pf_sparse,
        pf_group = pf_group,
        lower_bnd = lower_bnd,
        nfolds = nfolds,
        lambda = lambda_path,
        weights = obs_weights,
        ...
      ),
      warning = function(w) {
        backend_warnings <<- c(backend_warnings, conditionMessage(w))
      }
    )
  }, error = function(e) {
    # Make error messages more informative
    error_msg <- e$message
    
    if (grepl("singular|rank|collinear", error_msg, ignore.case = TRUE)) {
      stop("SGL fitting failed due to multicollinearity or rank deficiency.\n",
           "  Possible causes:\n",
           "  - Too many exposures relative to sample size (p >> n)\n",
           "  - Highly correlated variables in the design matrix\n",
           "  Try: Reduce the number of exposures or increase sample size.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else if (grepl("convergence|iterate", error_msg, ignore.case = TRUE)) {
      stop("SGL fitting failed to converge.\n",
           "  Possible causes:\n",
           "  - Data may have unusual scaling or outliers\n",
           "  - Lambda sequence may be inappropriate\n",
           "  Try: Standardize your data or adjust lambda settings.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else if (grepl("fold|cv|cross", error_msg, ignore.case = TRUE)) {
      stop("Cross-validation failed.\n",
           "  Possible causes:\n",
           "  - Sample size too small for ", nfolds, "-fold CV\n",
           "  - Some folds may have too few observations\n",
           "  Try: Reduce nfolds (e.g., nfolds = 5) or increase sample size.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else {
      stop("SGL fitting failed.\n",
           "  Sample size: ", n, ", Number of exposures: ", p, "\n",
           "  This may occur when:\n",
           "  - Sample size is too small relative to number of variables\n",
           "  - Data contains NA, Inf, or unusual values\n",
           "  Original error: ", error_msg, call. = FALSE)
    }
  })
  
  # Extract coefficients
  all_coef <- coef(fit, s = lambda)
  
  intercept <- all_coef[1]
  pos_coef <- as.numeric(all_coef[2:(p + 1)])
  neg_coef <- as.numeric(all_coef[(p + 2):(2 * p + 1)])
  
  names(pos_coef) <- var_names
  names(neg_coef) <- var_names
  
  # Covariate coefficients
  if (q > 0) {
    cov_coef <- as.numeric(all_coef[(2 * p + 2):length(all_coef)])
    names(cov_coef) <- cov_names
  } else {
    cov_coef <- NULL
  }

  selection_diagnostics <- .selection_diagnostics_from_fit(
    fit = fit,
    requested_lambda = lambda,
    lambda_path_source = if (is.null(lambda_path)) "sparsegl_auto" else "user_explicit",
    pos_coef = pos_coef,
    neg_coef = neg_coef,
    backend_warnings = backend_warnings
  )
  
  return(list(
    fit = fit,
    coefficients = all_coef,
    pos_coef = pos_coef,
    neg_coef = neg_coef,
    cov_coef = cov_coef,
    intercept = intercept,
    design_matrix = design_matrix,
    group = group,
    selection_diagnostics = selection_diagnostics
  ))
}


.validate_lambda_path <- function(lambda_path) {
  if (is.null(lambda_path)) {
    return(NULL)
  }
  if (!is.numeric(lambda_path) || length(lambda_path) < 2L ||
      any(!is.finite(lambda_path)) || any(lambda_path <= 0)) {
    stop(
      "`lambda_path` must be a numeric vector of at least two positive finite values.",
      call. = FALSE
    )
  }
  as.numeric(lambda_path)
}


#' @keywords internal
.selection_diagnostics_from_fit <- function(fit, requested_lambda,
                                            lambda_path_source = "sparsegl_auto",
                                            pos_coef, neg_coef,
                                            backend_warnings = character(0)) {
  lam <- as.numeric(fit$lambda %||% numeric(0))
  finite_lam <- lam[is.finite(lam)]
  selected <- if (is.numeric(requested_lambda)) {
    as.numeric(requested_lambda)[1]
  } else {
    as.numeric(fit[[requested_lambda]] %||% NA_real_)[1]
  }
  path_min <- if (length(finite_lam) > 0) min(finite_lam) else NA_real_
  path_max <- if (length(finite_lam) > 0) max(finite_lam) else NA_real_
  selected_on_edge <- is.finite(selected) && length(finite_lam) > 0 &&
    (isTRUE(all.equal(selected, path_min)) ||
       isTRUE(all.equal(selected, path_max)))
  nonzero_pos <- sum(abs(pos_coef) > 0, na.rm = TRUE)
  nonzero_neg <- sum(abs(neg_coef) > 0, na.rm = TRUE)

  list(
    requested_lambda = requested_lambda,
    selected_lambda = selected,
    lambda_path_source = lambda_path_source,
    lambda_path_length = length(finite_lam),
    lambda_path_min = path_min,
    lambda_path_max = path_max,
    selected_lambda_at_path_boundary = isTRUE(selected_on_edge),
    nonzero_positive_coef = nonzero_pos,
    nonzero_negative_coef = nonzero_neg,
    all_zero_exposure = (nonzero_pos + nonzero_neg) == 0L,
    backend_jerr = fit$sparsegl.fit$jerr %||% NA_integer_,
    backend_warnings = unique(backend_warnings),
    backend_error = NA_character_
  )
}


#' @keywords internal
.future_lapply <- function(...) {
  future.apply::future_lapply(...)
}


#' @keywords internal
.resolve_svrep_type <- function(survey_design, svrep_type = "auto") {
  if (!identical(svrep_type, "auto")) {
    return(svrep_type)
  }
  "bootstrap"
}


#' @keywords internal
.preserve_rng_seed <- function(seed = NULL) {
  if (is.null(seed)) {
    return(function() invisible(NULL))
  }
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(seed)
  function() {
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    invisible(NULL)
  }
}


#' @keywords internal
.sample_values <- function(x, size = length(x), replace = FALSE) {
  if (length(x) == 0L || size == 0L) {
    return(x[integer(0)])
  }
  x[sample.int(length(x), size = size, replace = replace)]
}


#' @keywords internal
.prepare_naive_boot_weights <- function(n, n_boot, family = "gaussian",
                                        y = NULL, stratified = TRUE,
                                        obs_weights = NULL, seed = NULL,
                                        survey_design_ignored = FALSE) {
  restore_seed <- .preserve_rng_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  index_matrix <- matrix(0L, nrow = n, ncol = n_boot)
  for (b in seq_len(n_boot)) {
    if (stratified && identical(family, "binomial") && !is.null(y)) {
      y_levels <- unique(y[!is.na(y)])
      idx <- integer(0)
      for (level in y_levels) {
        level_idx <- which(y == level)
        idx <- c(idx, .sample_values(level_idx, length(level_idx), replace = TRUE))
      }
      idx <- .sample_values(idx)
    } else {
      idx <- sample(n, replace = TRUE)
    }
    index_matrix[, b] <- idx
  }

  list(
    weight_matrix = NULL,
    index_matrix = index_matrix,
    obs_weights = obs_weights,
    n_boot_actual = n_boot,
    method_used = "naive",
    svrep_type_used = NA_character_,
    rep_design = NULL,
    bootstrap_design = if (isTRUE(survey_design_ignored)) {
      "naive_row_resampling_ignoring_psu_strata"
    } else {
      "naive_row_resampling"
    },
    survey_design_ignored = isTRUE(survey_design_ignored)
  )
}


#' @keywords internal
.prepare_svrep_boot_weights <- function(survey_design, svrep_type = "auto",
                                        svrep_args = list(), n_boot = 100,
                                        seed = NULL) {
  if (is.null(survey_design)) {
    stop("`boot_method = \"svrep\"` requires `survey_design`.", call. = FALSE)
  }
  .validate_survey_design(survey_design)
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop(
      "`boot_method = \"svrep\"` requires the `survey` package. ",
      "Install it with install.packages('survey').",
      call. = FALSE
    )
  }
  if (!is.list(svrep_args)) {
    stop("`svrep_args` must be a list.", call. = FALSE)
  }

  restore_seed <- .preserve_rng_seed(seed)
  on.exit(restore_seed(), add = TRUE)

  stochastic_types <- c("bootstrap", "subbootstrap", "mrbbootstrap")
  if (inherits(survey_design, "svyrep.design")) {
    if (!identical(svrep_type, "auto")) {
      message(
        "Provided `survey_design` is already a replicate design; ignoring ",
        "`svrep_type = \"", svrep_type, "\"`."
      )
    }
    rep_design <- survey_design
    svrep_type_used <- rep_design$type %||% "prebuilt"
  } else {
    svrep_type_used <- .resolve_svrep_type(survey_design, svrep_type)
    call_args <- list(
      design = survey_design,
      type = svrep_type_used,
      compress = FALSE
    )
    if (svrep_type_used %in% stochastic_types) {
      call_args$replicates <- n_boot
    }
    user_args <- svrep_args
    user_args$design <- NULL
    user_args$type <- NULL
    call_args <- utils::modifyList(call_args, user_args)
    rep_design <- tryCatch(
      do.call(survey::as.svrepdesign, call_args),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("Must use JK1 or bootstrap for an unstratified design", msg, fixed = TRUE)) {
          stop(
            "Survey replicate-weight bootstrap could not be constructed for ",
            "this unstratified design. Use `svrep_type = \"JK1\"` or ",
            "`svrep_type = \"bootstrap\"`, or set `boot_method = \"naive\"` ",
            "if ordinary weighted row bootstrap is intended. Original survey ",
            "error: ", msg,
            call. = FALSE
          )
        }
        stop(msg, call. = FALSE)
      }
    )
  }

  weight_matrix <- stats::weights(rep_design, type = "analysis")
  if (!is.matrix(weight_matrix)) {
    weight_matrix <- as.matrix(weight_matrix)
  }
  storage.mode(weight_matrix) <- "double"
  col_means <- colMeans(weight_matrix, na.rm = TRUE)
  col_means[!is.finite(col_means) | col_means <= 0] <- 1
  weight_matrix <- sweep(weight_matrix, 2, col_means, "/")

  n_rep <- ncol(weight_matrix)
  n_boot_actual <- if (svrep_type_used %in% stochastic_types) {
    min(n_rep, n_boot)
  } else {
    n_rep
  }
  variance_scale <- rep_design$scale %||% 1
  if (n_boot_actual < n_rep) {
    variance_scale <- variance_scale * n_rep / n_boot_actual
  }
  weight_matrix <- weight_matrix[, seq_len(n_boot_actual), drop = FALSE]
  variance_rscales <- rep_design$rscales %||% rep(1, n_rep)
  variance_rscales <- as.numeric(variance_rscales)
  if (length(variance_rscales) == 1L) {
    variance_rscales <- rep(variance_rscales, n_rep)
  }
  variance_rscales <- variance_rscales[seq_len(n_boot_actual)]
  variance_df <- tryCatch(
    as.numeric(survey::degf(rep_design)),
    error = function(e) NA_real_
  )

  list(
    weight_matrix = weight_matrix,
    index_matrix = NULL,
    obs_weights = NULL,
    n_boot_actual = n_boot_actual,
    method_used = "svrep",
    svrep_type_used = svrep_type_used,
    rep_design = rep_design,
    variance_scale = variance_scale,
    variance_rscales = variance_rscales,
    variance_df = variance_df,
    weights_normalized = TRUE,
    bootstrap_design = "survey_replicate_weights",
    survey_design_ignored = FALSE
  )
}


#' @keywords internal
.prepare_boot_weights <- function(n, n_boot, boot_method = c("auto", "naive", "svrep"),
                                  survey_design = NULL,
                                  svrep_type = c("auto", "JK1", "JKn", "BRR", "Fay",
                                                 "bootstrap", "subbootstrap", "mrbbootstrap"),
                                  svrep_args = list(),
                                  obs_weights = NULL,
                                  family = "gaussian",
                                  y = NULL,
                                  stratified = TRUE,
                                  seed = NULL) {
  boot_method <- match.arg(boot_method)
  svrep_type <- match.arg(svrep_type)
  if (identical(boot_method, "auto")) {
    boot_method <- if (!is.null(survey_design)) "svrep" else "naive"
  }
  if (identical(boot_method, "svrep")) {
    return(.prepare_svrep_boot_weights(
      survey_design = survey_design,
      svrep_type = svrep_type,
      svrep_args = svrep_args,
      n_boot = n_boot,
      seed = seed
    ))
  }
  .prepare_naive_boot_weights(
    n = n,
    n_boot = n_boot,
    family = family,
    y = y,
    stratified = stratified,
    obs_weights = obs_weights,
    seed = seed,
    survey_design_ignored = !is.null(survey_design)
  )
}


#' @keywords internal
.bootstrap_se <- function(values, center = NULL, boot_success = NULL,
                          boot_prep = NULL, na.rm = FALSE) {
  is_vector <- is.null(dim(values))
  if (is_vector) {
    values <- matrix(values, ncol = 1L)
  }
  if (is.null(boot_success)) {
    boot_success <- rep(TRUE, nrow(values))
  }
  block <- values[boot_success, , drop = FALSE]

  if (is.null(center)) {
    center <- colMeans(block, na.rm = na.rm)
  }
  if (length(center) == 1L && ncol(block) > 1L) {
    center <- rep(center, ncol(block))
  }

  if (is.null(boot_prep) || !identical(boot_prep$method_used, "svrep")) {
    out <- apply(block, 2, stats::sd, na.rm = na.rm)
    return(if (is_vector) out[[1]] else out)
  }

  scale <- boot_prep$variance_scale %||% 1
  rscales <- boot_prep$variance_rscales %||% rep(1, nrow(values))
  rscales <- as.numeric(rscales)[boot_success]
  if (length(rscales) == 1L && nrow(block) > 1L) {
    rscales <- rep(rscales, nrow(block))
  }

  out <- vapply(seq_len(ncol(block)), function(j) {
    x <- block[, j]
    valid <- is.finite(x) & is.finite(rscales)
    if (!na.rm && any(!valid)) {
      return(NA_real_)
    }
    x <- x[valid]
    rs <- rscales[valid]
    ctr <- center[[j]]
    if (!is.finite(ctr)) {
      ctr <- mean(x, na.rm = TRUE)
    }
    sqrt(scale * sum(rs * (x - ctr)^2, na.rm = TRUE))
  }, numeric(1))

  if (is_vector) out[[1]] else out
}


#' @keywords internal
.svrep_ci_multiplier <- function(x, conf_level = 0.95) {
  df <- x$variance_df %||% x$survey_df %||% NA_real_
  alpha <- 1 - conf_level
  if (is.finite(df) && df > 0) {
    stats::qt(1 - alpha / 2, df = df)
  } else {
    stats::qnorm(1 - alpha / 2)
  }
}


#' @keywords internal
.classify_boot_error <- function(msg) {
  msg <- msg %||% ""
  if (!nzchar(msg)) {
    return("unknown")
  }
  if (grepl("boundary|sparsegl_irls|converge|convergence|iterate", msg, ignore.case = TRUE)) {
    return("backend_convergence")
  }
  if (grepl("lambda", msg, ignore.case = TRUE)) {
    return("lambda_path")
  }
  if (grepl("fold|cv|cross", msg, ignore.case = TRUE)) {
    return("cross_validation")
  }
  if (grepl("singular|rank|collinear", msg, ignore.case = TRUE)) {
    return("rank_deficiency")
  }
  "other"
}

#' Internal Function: Bootstrap Aggregation for Weights
#'
#' Performs bootstrap aggregation to stabilize weight estimates.
#' Supports parallel processing via the future framework.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Processed covariate matrix.
#' @param var_names Exposure variable names.
#' @param cov_names Covariate names.
#' @param groups Group list.
#' @param group_by_compound Logical.
#' @param group_structure Character.
#' @param penalize_covariates Logical.
#' @param family Character.
#' @param lambda Character or numeric.
#' @param lambda_path Optional numeric lambda sequence passed to cv.sparsegl.
#' @param nfolds Integer.
#' @param n_boot Number of bootstrap iterations.
#' @param seed Random seed.
#' @param verbose Logical.
#' @param parallel Logical. Whether to use parallel processing (requires future.apply).
#' @param ... Additional arguments.
#'
#' @return A list with aggregated results.
#'
#' @keywords internal
bootstrap_sgl <- function(X_quantile, y, cov_matrix, var_names, cov_names,
                           groups, group_by_compound, group_structure,
                           penalize_covariates, family, lambda, nfolds,
                           lambda_path = NULL,
                           n_boot, seed, verbose, parallel = FALSE,
                           checkpoint_dir = NULL, checkpoint_interval = 50,
                           cleanup_checkpoint = TRUE, 
                           stratified = TRUE, obs_weights = NULL,
                           boot_method = c("auto", "naive", "svrep"),
                           survey_design = NULL,
                           svrep_type = c("auto", "JK1", "JKn", "BRR", "Fay",
                                          "bootstrap", "subbootstrap", "mrbbootstrap"),
                           svrep_args = list(), ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(X_quantile)
  p <- ncol(X_quantile)
  q_cov <- length(cov_names)
  obs_weights <- .validate_obs_weights(obs_weights, n)
  lambda_path <- .validate_lambda_path(lambda_path)
  boot_prep <- .prepare_boot_weights(
    n = n,
    n_boot = n_boot,
    boot_method = boot_method,
    survey_design = survey_design,
    svrep_type = svrep_type,
    svrep_args = svrep_args,
    obs_weights = obs_weights,
    family = family,
    y = y,
    stratified = stratified,
    seed = seed
  )
  n_boot_requested <- n_boot
  n_boot <- boot_prep$n_boot_actual
  if (verbose && !identical(boot_prep$method_used, "naive") &&
      !identical(n_boot, n_boot_requested)) {
    message(sprintf(
      "Using %d survey replicates as determined by svrep_type = \"%s\" (requested n_boot = %d).",
      n_boot, boot_prep$svrep_type_used, n_boot_requested
    ))
  }
  
  # Prepare checkpoint functionality
  use_checkpoint <- !is.null(checkpoint_dir)
  checkpoint_file <- NULL
  start_boot <- 1
  boot_pos_coef <- NULL
  boot_neg_coef <- NULL
  boot_cov_coef <- NULL
  boot_pos_index_sum <- NULL
  boot_neg_index_sum <- NULL
  boot_pos_index_sum_by_group <- NULL
  boot_neg_index_sum_by_group <- NULL
  boot_success <- NULL
  boot_error_msg <- NULL
  boot_error_class <- NULL
  completed_boots <- integer(0)
  batch_errors <- character(0)
  
  if (use_checkpoint) {
    if (!dir.exists(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE)
    }
    checkpoint_file <- file.path(checkpoint_dir, "bootstrap_checkpoint.rds")
    
    # Load existing checkpoint if available
    if (file.exists(checkpoint_file)) {
      if (verbose) message("Loading checkpoint from ", checkpoint_file)
      checkpoint_data <- readRDS(checkpoint_file)
      
      # Validate checkpoint data
      if (!is.null(checkpoint_data$boot_pos_coef) && 
          ncol(checkpoint_data$boot_pos_coef) == p &&
          checkpoint_data$n_boot == n_boot &&
          identical(checkpoint_data$boot_method %||% "naive", boot_prep$method_used) &&
          identical(checkpoint_data$svrep_type_used %||% NA_character_, boot_prep$svrep_type_used) &&
          identical(checkpoint_data$svrep_args %||% list(), svrep_args)) {
        boot_pos_coef <- checkpoint_data$boot_pos_coef
        boot_neg_coef <- checkpoint_data$boot_neg_coef
        boot_cov_coef <- checkpoint_data$boot_cov_coef
        boot_pos_index_sum <- checkpoint_data$boot_pos_index_sum
        boot_neg_index_sum <- checkpoint_data$boot_neg_index_sum
        boot_pos_index_sum_by_group <- checkpoint_data$boot_pos_index_sum_by_group
        boot_neg_index_sum_by_group <- checkpoint_data$boot_neg_index_sum_by_group
        boot_success <- checkpoint_data$boot_success
        boot_error_msg <- checkpoint_data$boot_error_msg
        boot_error_class <- checkpoint_data$boot_error_class %||% NULL
        # Older checkpoints may have marked failed batch indices as completed.
        # Trust the success vector so failed iterations are retried.
        completed_boots <- which(!is.na(boot_success) & boot_success)
        start_boot <- if (length(completed_boots) > 0) max(completed_boots) + 1 else 1
        if (verbose) message("Resuming from iteration ", start_boot, 
                            " (", length(completed_boots), " completed)")
      } else {
        if (verbose) message("Checkpoint incompatible, starting fresh")
      }
    }
  }
  
  # Initialize result matrices (if not resuming from checkpoint)
  if (is.null(boot_pos_coef)) {
    boot_pos_coef <- matrix(0, nrow = n_boot, ncol = p)
    boot_neg_coef <- matrix(0, nrow = n_boot, ncol = p)
    colnames(boot_pos_coef) <- var_names
    colnames(boot_neg_coef) <- var_names
  }
  if (q_cov > 0 && is.null(boot_cov_coef)) {
    boot_cov_coef <- matrix(NA_real_, nrow = n_boot, ncol = q_cov)
    colnames(boot_cov_coef) <- cov_names
  } else if (q_cov == 0) {
    boot_cov_coef <- NULL
  }
  if (is.null(boot_pos_index_sum)) {
    boot_pos_index_sum <- rep(NA_real_, n_boot)
  }
  if (is.null(boot_neg_index_sum)) {
    boot_neg_index_sum <- rep(NA_real_, n_boot)
  }
  if (!is.null(groups)) {
    group_names <- names(groups)
    if (is.null(boot_pos_index_sum_by_group)) {
      boot_pos_index_sum_by_group <- stats::setNames(
        lapply(group_names, function(...) rep(NA_real_, n_boot)),
        group_names
      )
    }
    if (is.null(boot_neg_index_sum_by_group)) {
      boot_neg_index_sum_by_group <- stats::setNames(
        lapply(group_names, function(...) rep(NA_real_, n_boot)),
        group_names
      )
    }
  } else {
    boot_pos_index_sum_by_group <- NULL
    boot_neg_index_sum_by_group <- NULL
  }

  # Bootstrap success tracking vector (restore from checkpoint or create new)
  if (is.null(boot_success)) {
    boot_success <- logical(n_boot)
  }
  if (is.null(boot_error_msg)) {
    boot_error_msg <- rep(NA_character_, n_boot)
  }
  if (is.null(boot_error_class)) {
    boot_error_class <- rep(NA_character_, n_boot)
  }
  
  # Remaining bootstrap indices
  remaining_boots <- setdiff(seq_len(n_boot), completed_boots)
  
  if (length(remaining_boots) == 0) {
    if (verbose) message("All bootstrap iterations already completed")
  } else {
    
    # Single bootstrap execution function
    run_single_boot <- function(b, X_quantile, y, cov_matrix, var_names, cov_names,
                                 groups, group_by_compound, group_structure,
                                 penalize_covariates, family, lambda, nfolds,
                                 lambda_path = NULL,
                                 stratified = TRUE, obs_weights = NULL,
                                 boot_prep, ...) {
      
      n <- nrow(X_quantile)
      
      if (identical(boot_prep$method_used, "svrep")) {
        X_boot <- X_quantile
        y_boot <- y
        cov_boot <- cov_matrix
        weights_boot <- boot_prep$weight_matrix[, b]
      } else {
        boot_idx <- boot_prep$index_matrix[, b]
        X_boot <- X_quantile[boot_idx, , drop = FALSE]
        y_boot <- y[boot_idx]
        cov_boot <- if (!is.null(cov_matrix)) cov_matrix[boot_idx, , drop = FALSE] else NULL
        weights_boot <- if (!is.null(boot_prep$obs_weights)) boot_prep$obs_weights[boot_idx] else NULL
      }
      
      # Fitting (with error handling)
      tryCatch({
        fit_result <- fit_sgl_core(
          X_quantile = X_boot,
          y = y_boot,
          cov_matrix = cov_boot,
          var_names = var_names,
          cov_names = cov_names,
          groups = groups,
          group_by_compound = group_by_compound,
          group_structure = group_structure,
          penalize_covariates = penalize_covariates,
          family = family,
          lambda = lambda,
          nfolds = nfolds,
          lambda_path = lambda_path,
          obs_weights = weights_boot,
          ...
        )

        weight_result <- calculate_weights(
          pos_coef = fit_result$pos_coef,
          neg_coef = fit_result$neg_coef,
          var_names = var_names,
          groups = groups,
          handle_collinearity = "net"
        )
        
        return(list(
          success = TRUE,
          pos_coef = fit_result$pos_coef,
          neg_coef = fit_result$neg_coef,
          cov_coef = fit_result$cov_coef,
          pos_index_sum = weight_result$pos_index_sum,
          neg_index_sum = weight_result$neg_index_sum,
          pos_index_sum_by_group = weight_result$pos_index_sum_by_group,
          neg_index_sum_by_group = weight_result$neg_index_sum_by_group
        ))
        
      }, error = function(e) {
        return(list(
          success = FALSE,
          pos_coef = rep(0, length(var_names)),
          neg_coef = rep(0, length(var_names)),
          cov_coef = if (!is.null(cov_names)) {
            stats::setNames(rep(NA_real_, length(cov_names)), cov_names)
          } else {
            NULL
          },
          pos_index_sum = NA_real_,
          neg_index_sum = NA_real_,
          pos_index_sum_by_group = if (!is.null(groups)) {
            stats::setNames(as.list(rep(NA_real_, length(groups))), names(groups))
          } else {
            NULL
          },
          neg_index_sum_by_group = if (!is.null(groups)) {
            stats::setNames(as.list(rep(NA_real_, length(groups))), names(groups))
          } else {
            NULL
          },
          error_msg = e$message
        ))
      })
    }
    
    # Checkpoint save function
    append_completed_boots <- function(indices) {
      if (length(indices) == 0) {
        return(invisible(NULL))
      }
      # indices are only added once so unique is unnecessary (avoids O(n^2))
      completed_boots <<- c(completed_boots, indices)
      invisible(NULL)
    }

    store_boot_result <- function(b, result) {
      if (isTRUE(result$success)) {
        boot_pos_coef[b, ] <<- result$pos_coef
        boot_neg_coef[b, ] <<- result$neg_coef
        if (!is.null(boot_cov_coef) && !is.null(result$cov_coef)) {
          boot_cov_coef[b, ] <<- result$cov_coef
        }
        boot_pos_index_sum[b] <<- result$pos_index_sum
        boot_neg_index_sum[b] <<- result$neg_index_sum
        if (!is.null(boot_pos_index_sum_by_group) &&
            !is.null(result$pos_index_sum_by_group)) {
          for (grp in names(boot_pos_index_sum_by_group)) {
            boot_pos_index_sum_by_group[[grp]][b] <<-
              result$pos_index_sum_by_group[[grp]]
            boot_neg_index_sum_by_group[[grp]][b] <<-
              result$neg_index_sum_by_group[[grp]]
          }
        }
        boot_success[b] <<- TRUE
        boot_error_msg[b] <<- NA_character_
        boot_error_class[b] <<- NA_character_
        return(TRUE)
      }

      boot_error_msg[b] <<- result$error_msg %||% "Unknown bootstrap failure"
      boot_error_class[b] <<- .classify_boot_error(boot_error_msg[b])
      FALSE
    }

    run_sequential_boot_batch <- function(batch_indices) {
      batch_completed <- integer(0)
      for (b in batch_indices) {
        result <- run_single_boot(
          b, X_quantile, y, cov_matrix, var_names, cov_names,
          groups, group_by_compound, group_structure,
          penalize_covariates, family, lambda, nfolds,
          lambda_path = lambda_path,
          stratified = stratified,
          obs_weights = obs_weights,
          boot_prep = boot_prep,
          ...
        )
        if (store_boot_result(b, result)) {
          batch_completed <- c(batch_completed, b)
        }
      }
      append_completed_boots(batch_completed)
      batch_completed
    }
    
    save_checkpoint <- function() {
      if (use_checkpoint) {
        checkpoint_data <- list(
          boot_pos_coef = boot_pos_coef,
          boot_neg_coef = boot_neg_coef,
          boot_cov_coef = boot_cov_coef,
          boot_pos_index_sum = boot_pos_index_sum,
          boot_neg_index_sum = boot_neg_index_sum,
          boot_pos_index_sum_by_group = boot_pos_index_sum_by_group,
          boot_neg_index_sum_by_group = boot_neg_index_sum_by_group,
          boot_success = boot_success,
          boot_error_msg = boot_error_msg,
          boot_error_class = boot_error_class,
          completed_boots = completed_boots,
          n_boot = n_boot,
          n_boot_requested = n_boot_requested,
          boot_method = boot_prep$method_used,
          svrep_type_used = boot_prep$svrep_type_used,
          svrep_args = svrep_args,
          bootstrap_design = boot_prep$bootstrap_design %||% NA_character_,
          survey_design_ignored = boot_prep$survey_design_ignored %||% FALSE,
          variance_scale = boot_prep$variance_scale,
          variance_rscales = boot_prep$variance_rscales,
          var_names = var_names,
          timestamp = Sys.time()
        )
        saveRDS(checkpoint_data, checkpoint_file)
        if (verbose) {
          if (length(completed_boots) > 0) {
            message("  Checkpoint saved at iteration ", max(completed_boots))
          } else {
            message("  Checkpoint saved with no successful iterations yet")
          }
        }
      }
    }
    
    # Execute parallel processing (checkpoint per batch)
    if (parallel) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        warning("Package 'future.apply' is not installed. Falling back to sequential processing.\n",
                "Install with: install.packages('future.apply')")
        parallel <- FALSE
      } else {
        if (verbose) {
          message("Running bootstrap aggregation with ", length(remaining_boots), 
                  " remaining iterations (parallel)...")
          if (use_checkpoint) {
            message("Checkpointing every ", checkpoint_interval, " iterations to ", checkpoint_dir)
          }
        }
        
        # Checkpoint via batch processing
        batch_start <- 1
        parallel_degraded <- FALSE
        
        while (batch_start <= length(remaining_boots)) {
          batch_end <- min(batch_start + checkpoint_interval - 1, length(remaining_boots))
          batch_indices <- remaining_boots[batch_start:batch_end]
          
          if (parallel_degraded) {
            batch_completed <- run_sequential_boot_batch(batch_indices)
            if (verbose) {
              message("  Parallel backend already degraded; ran batch sequentially (",
                      length(batch_completed), "/", length(batch_indices),
                      " recovered).")
            }
          } else {
            # Execute batch (catch batch-level errors too)
            batch_result <- tryCatch({
              batch_results <- .future_lapply(
                batch_indices,
                function(b) {
                  run_single_boot(b, X_quantile, y, cov_matrix, var_names, cov_names,
                                  groups, group_by_compound, group_structure,
                                  penalize_covariates, family, lambda, nfolds,
                                  lambda_path = lambda_path,
                                  stratified = stratified,
                                  obs_weights = obs_weights,
                                  boot_prep = boot_prep,
                                  ...)
                },
                future.seed = TRUE
              )
              list(success = TRUE, results = batch_results)
            }, error = function(e) {
              list(success = FALSE, error_msg = e$message)
            })

            if (batch_result$success) {
              # Store results
              batch_completed <- integer(0)
              for (i in seq_along(batch_indices)) {
                b <- batch_indices[i]
                if (store_boot_result(b, batch_result$results[[i]])) {
                  batch_completed <- c(batch_completed, b)
                }
              }
              append_completed_boots(batch_completed)

              # Free memory
              rm(batch_result)
              gc(verbose = FALSE)
            } else {
              # Entire batch failed
              batch_errors <- c(batch_errors, batch_result$error_msg)
              parallel_degraded <- TRUE
              if (verbose) {
                message("  Warning: Batch ", batch_start, "-", batch_end, " failed: ",
                        substr(batch_result$error_msg, 1, 100))
                message("  Retrying failed parallel batch sequentially and degrading ",
                        "remaining batches to sequential.")
              }
              batch_completed <- run_sequential_boot_batch(batch_indices)
              if (verbose) {
                message("  Sequential retry recovered ", length(batch_completed), "/",
                        length(batch_indices), " bootstrap iteration(s).")
              }
            }
          }
          
          # Save checkpoint
          save_checkpoint()

          batch_start <- batch_end + 1
        }
        
        if (length(batch_errors) > 0 && verbose) {
          message("Note: ", length(batch_errors), " batch(es) failed during parallel processing.")
          message("Unique errors: ", paste(unique(batch_errors)[1:min(3, length(unique(batch_errors)))], collapse = "; "))
        }
      }
    }
    
    # Sequential processing
    if (!parallel) {
      if (verbose) {
        message("Running bootstrap aggregation with ", length(remaining_boots), 
                " remaining iterations...")
        if (use_checkpoint) {
          message("Checkpointing every ", checkpoint_interval, " iterations to ", checkpoint_dir)
        }
        pb <- txtProgressBar(min = 0, max = length(remaining_boots), style = 3)
      }
      
      for (i in seq_along(remaining_boots)) {
        b <- remaining_boots[i]
        run_sequential_boot_batch(b)
        
        # Save checkpoint
        if (use_checkpoint && (i %% checkpoint_interval == 0 || i == length(remaining_boots))) {
          save_checkpoint()
        }
        
        if (verbose) setTxtProgressBar(pb, i)
      }
      
      if (verbose) close(pb)
    }
  }
  
  # Aggregate results (success determined by boot_success vector — all-zero coefficients included as valid)
  successful_boots <- sum(boot_success)

  if (verbose) {
    message("Completed ", successful_boots, "/", n_boot, " bootstrap iterations.")
  }

  if (successful_boots < n_boot * 0.5) {
    warning(
      sprintf(
        "Successful bootstrap iterations: %d / %d (%.1f%%).",
        successful_boots, n_boot, 100 * successful_boots / n_boot
      ),
      call. = FALSE
    )
  }
  if (identical(boot_prep$method_used, "svrep") && successful_boots < n_boot) {
    warning(
      sprintf(
        "Survey replicate bootstrap had %d failed replicate(s); variance estimates use only %d/%d successful replicates.",
        n_boot - successful_boots, successful_boots, n_boot
      ),
      call. = FALSE
    )
  }

  # Mean coefficients (successful bootstraps only — including all-zero results)
  if (successful_boots > 0) {
    mean_pos_coef <- colMeans(boot_pos_coef[boot_success, , drop = FALSE])
    mean_neg_coef <- colMeans(boot_neg_coef[boot_success, , drop = FALSE])

    svrep_center <- NULL
    if (identical(boot_prep$method_used, "svrep")) {
      svrep_center <- tryCatch({
        center_fit <- fit_sgl_core(
          X_quantile = X_quantile,
          y = y,
          cov_matrix = cov_matrix,
          var_names = var_names,
          cov_names = cov_names,
          groups = groups,
          group_by_compound = group_by_compound,
          group_structure = group_structure,
          penalize_covariates = penalize_covariates,
          family = family,
          lambda = lambda,
          nfolds = nfolds,
          lambda_path = lambda_path,
          obs_weights = obs_weights,
          ...
        )
        center_weights <- calculate_weights(
          pos_coef = center_fit$pos_coef,
          neg_coef = center_fit$neg_coef,
          var_names = var_names,
          groups = groups,
          handle_collinearity = "net"
        )
        list(
          pos_coef = center_fit$pos_coef,
          neg_coef = center_fit$neg_coef,
          cov_coef = center_fit$cov_coef,
          pos_index_sum = center_weights$pos_index_sum,
          neg_index_sum = center_weights$neg_index_sum,
          pos_index_sum_by_group = center_weights$pos_index_sum_by_group,
          neg_index_sum_by_group = center_weights$neg_index_sum_by_group
        )
      }, error = function(e) NULL)
    }

    if (!is.null(svrep_center)) {
      mean_pos_coef <- svrep_center$pos_coef
      mean_neg_coef <- svrep_center$neg_coef
    }

    # Standard errors
    se_pos_coef <- .bootstrap_se(
      boot_pos_coef,
      center = if (!is.null(svrep_center)) svrep_center$pos_coef else mean_pos_coef,
      boot_success = boot_success,
      boot_prep = boot_prep
    )
    se_neg_coef <- .bootstrap_se(
      boot_neg_coef,
      center = if (!is.null(svrep_center)) svrep_center$neg_coef else mean_neg_coef,
      boot_success = boot_success,
      boot_prep = boot_prep
    )
    names(se_pos_coef) <- var_names
    names(se_neg_coef) <- var_names

    # Selection frequency (proportion of non-zero across all successful bootstraps)
    selection_freq_pos <- colMeans(boot_pos_coef[boot_success, , drop = FALSE] > 0)
    selection_freq_neg <- colMeans(boot_neg_coef[boot_success, , drop = FALSE] > 0)

    index_sum_pos_block <- boot_pos_index_sum[boot_success]
    index_sum_neg_block <- boot_neg_index_sum[boot_success]
    mean_index_sum_pos <- mean(index_sum_pos_block, na.rm = TRUE)
    mean_index_sum_neg <- mean(index_sum_neg_block, na.rm = TRUE)
    if (!is.null(svrep_center)) {
      if (!is.null(svrep_center$pos_index_sum)) {
        mean_index_sum_pos <- svrep_center$pos_index_sum
      }
      if (!is.null(svrep_center$neg_index_sum)) {
        mean_index_sum_neg <- svrep_center$neg_index_sum
      }
    }
    se_index_sum_pos <- .bootstrap_se(
      boot_pos_index_sum,
      center = if (!is.null(svrep_center)) svrep_center$pos_index_sum else mean_index_sum_pos,
      boot_success = boot_success,
      boot_prep = boot_prep,
      na.rm = TRUE
    )
    se_index_sum_neg <- .bootstrap_se(
      boot_neg_index_sum,
      center = if (!is.null(svrep_center)) svrep_center$neg_index_sum else mean_index_sum_neg,
      boot_success = boot_success,
      boot_prep = boot_prep,
      na.rm = TRUE
    )

    if (!is.null(boot_pos_index_sum_by_group)) {
      mean_index_sum_by_group_pos <- lapply(boot_pos_index_sum_by_group, function(x) {
        mean(x[boot_success], na.rm = TRUE)
      })
      mean_index_sum_by_group_neg <- lapply(boot_neg_index_sum_by_group, function(x) {
        mean(x[boot_success], na.rm = TRUE)
      })
      if (!is.null(svrep_center) && !is.null(svrep_center$pos_index_sum_by_group)) {
        mean_index_sum_by_group_pos <- svrep_center$pos_index_sum_by_group
      }
      if (!is.null(svrep_center) && !is.null(svrep_center$neg_index_sum_by_group)) {
        mean_index_sum_by_group_neg <- svrep_center$neg_index_sum_by_group
      }
      se_index_sum_by_group_pos <- Map(function(x, ctr) {
        .bootstrap_se(
          x,
          center = ctr,
          boot_success = boot_success,
          boot_prep = boot_prep,
          na.rm = TRUE
        )
      }, boot_pos_index_sum_by_group, mean_index_sum_by_group_pos)
      se_index_sum_by_group_neg <- Map(function(x, ctr) {
        .bootstrap_se(
          x,
          center = ctr,
          boot_success = boot_success,
          boot_prep = boot_prep,
          na.rm = TRUE
        )
      }, boot_neg_index_sum_by_group, mean_index_sum_by_group_neg)
    } else {
      mean_index_sum_by_group_pos <- NULL
      mean_index_sum_by_group_neg <- NULL
      se_index_sum_by_group_pos <- NULL
      se_index_sum_by_group_neg <- NULL
    }

    if (!is.null(boot_cov_coef)) {
      cov_block <- boot_cov_coef[boot_success, , drop = FALSE]
      mean_cov_coef <- colMeans(cov_block, na.rm = TRUE)
      if (!is.null(svrep_center) && !is.null(svrep_center$cov_coef) &&
          length(svrep_center$cov_coef) == length(mean_cov_coef)) {
        mean_cov_coef <- svrep_center$cov_coef
      }
      se_cov_coef <- .bootstrap_se(
        boot_cov_coef,
        center = mean_cov_coef,
        boot_success = boot_success,
        boot_prep = boot_prep,
        na.rm = TRUE
      )
      names(se_cov_coef) <- cov_names
      if (identical(boot_prep$method_used, "svrep")) {
        svrep_crit <- .svrep_ci_multiplier(boot_prep, conf_level = 0.95)
        ci_lower_cov <- mean_cov_coef - svrep_crit * se_cov_coef
        ci_upper_cov <- mean_cov_coef + svrep_crit * se_cov_coef
      } else {
        ci_lower_cov <- apply(cov_block, 2, quantile, probs = 0.025, na.rm = TRUE)
        ci_upper_cov <- apply(cov_block, 2, quantile, probs = 0.975, na.rm = TRUE)
      }
    } else {
      mean_cov_coef <- NULL
      se_cov_coef <- NULL
      ci_lower_cov <- NULL
      ci_upper_cov <- NULL
    }
  } else {
    representative_errors <- unique(stats::na.omit(c(boot_error_msg, batch_errors)))
    if (length(representative_errors) > 0) {
      stop(
        "All bootstrap iterations failed. Representative error(s): ",
        paste(utils::head(representative_errors, 3), collapse = "; "),
        call. = FALSE
      )
    }
    stop("All bootstrap iterations failed.", call. = FALSE)
  }
  
  # Remove checkpoint file
  if (use_checkpoint && cleanup_checkpoint && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    if (verbose) message("Checkpoint file removed")
  }
  
  return(list(
    mean_pos_coef = mean_pos_coef,
    mean_neg_coef = mean_neg_coef,
    se_pos_coef = se_pos_coef,
    se_neg_coef = se_neg_coef,
    selection_freq_pos = selection_freq_pos,
    selection_freq_neg = selection_freq_neg,
    mean_index_sum_pos = mean_index_sum_pos,
    mean_index_sum_neg = mean_index_sum_neg,
    se_index_sum_pos = se_index_sum_pos,
    se_index_sum_neg = se_index_sum_neg,
    mean_index_sum_by_group_pos = mean_index_sum_by_group_pos,
    mean_index_sum_by_group_neg = mean_index_sum_by_group_neg,
    se_index_sum_by_group_pos = se_index_sum_by_group_pos,
    se_index_sum_by_group_neg = se_index_sum_by_group_neg,
    boot_pos_coef = boot_pos_coef,
    boot_neg_coef = boot_neg_coef,
    mean_cov_coef = mean_cov_coef,
    se_cov_coef = se_cov_coef,
    ci_lower_cov = ci_lower_cov,
    ci_upper_cov = ci_upper_cov,
    boot_cov_coef = boot_cov_coef,
    boot_pos_index_sum = boot_pos_index_sum,
    boot_neg_index_sum = boot_neg_index_sum,
    boot_pos_index_sum_by_group = boot_pos_index_sum_by_group,
    boot_neg_index_sum_by_group = boot_neg_index_sum_by_group,
    boot_success = boot_success,
    boot_error_msg = boot_error_msg,
    boot_error_class = boot_error_class,
    boot_error_counts = sort(table(stats::na.omit(boot_error_class)), decreasing = TRUE),
    parallel_batch_errors = unique(batch_errors),
    n_parallel_batch_failures = length(batch_errors),
    n_successful = successful_boots,
    n_failed = n_boot - successful_boots,
    method = boot_prep$method_used,
    svrep_type_used = boot_prep$svrep_type_used,
    svrep_args = svrep_args,
    bootstrap_design = boot_prep$bootstrap_design %||% NA_character_,
    survey_design_ignored = boot_prep$survey_design_ignored %||% FALSE,
    variance_scale = boot_prep$variance_scale,
    variance_rscales = boot_prep$variance_rscales,
    variance_df = boot_prep$variance_df,
    svrep_center_pos_coef = if (!is.null(svrep_center)) svrep_center$pos_coef else NULL,
    svrep_center_neg_coef = if (!is.null(svrep_center)) svrep_center$neg_coef else NULL,
    svrep_center_cov_coef = if (!is.null(svrep_center)) svrep_center$cov_coef else NULL,
    svrep_center_index_sum_pos = if (!is.null(svrep_center)) svrep_center$pos_index_sum else NULL,
    svrep_center_index_sum_neg = if (!is.null(svrep_center)) svrep_center$neg_index_sum else NULL,
    svrep_center_index_sum_by_group_pos = if (!is.null(svrep_center)) svrep_center$pos_index_sum_by_group else NULL,
    svrep_center_index_sum_by_group_neg = if (!is.null(svrep_center)) svrep_center$neg_index_sum_by_group else NULL,
    weights_normalized = boot_prep$weights_normalized %||% FALSE,
    n_boot_actual = n_boot,
    n_boot_requested = n_boot_requested
  ))
}


#' Internal Function: Refit WQS Model for Inference
#'
#' Fits a downstream inference model using fixed WQS indices.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Covariate matrix.
#' @param pos_weights Fixed positive weights.
#' @param neg_weights Fixed negative weights.
#' @param family Character.
#' @param groups Group list.
#' @param group_inference Logical. Whether to include group-specific indices.
#' @param engine Character. One of \code{"glm"} or \code{"svyglm"}.
#' @param survey_design Optional pre-constructed survey design object.
#' @param analysis_id Optional analysis-row identifier used to verify alignment
#'   between \code{glm_data} and \code{survey_design}.
#'
#' @return A list with refit results and extracted inference summaries.
#'
#' @keywords internal
refit_model <- function(X_quantile, y, cov_matrix,
                        pos_weights, neg_weights, family, groups,
                        group_inference = TRUE,
                        engine = c("glm", "svyglm"),
                        survey_design = NULL,
                        analysis_id = NULL,
                        obs_weights = NULL,
                        excluded_directions = list()) {
  engine <- match.arg(engine)
  var_names <- names(pos_weights)
  n_obs <- length(y)
  
  if (!is.null(groups) && group_inference) {
    original_grp_names <- names(groups)
    safe_grp_names <- make.names(original_grp_names, unique = TRUE)
    grp_name_to_safe <- stats::setNames(safe_grp_names, original_grp_names)
    
    group_indices <- data.frame(matrix(nrow = n_obs, ncol = 0))
    group_names_pos <- character(0)
    group_names_neg <- character(0)
    
    for (grp_name in original_grp_names) {
      grp_vars <- groups[[grp_name]]
      grp_idx <- which(var_names %in% grp_vars)
      if (length(grp_idx) == 0) next

      grp_pos_w <- pos_weights[grp_idx]
      grp_neg_w <- neg_weights[grp_idx]
      if (sum(grp_pos_w) > 0) grp_pos_w <- grp_pos_w / sum(grp_pos_w)
      if (sum(grp_neg_w) > 0) grp_neg_w <- grp_neg_w / sum(grp_neg_w)

      X_grp <- X_quantile[, grp_idx, drop = FALSE]
      safe_name <- grp_name_to_safe[[grp_name]]

      pos_excluded <- identical(excluded_directions[[grp_name]], "positive")
      neg_excluded <- identical(excluded_directions[[grp_name]], "negative")

      if (!pos_excluded) {
        pos_col_name <- paste0(safe_name, "_pos")
        group_indices[[pos_col_name]] <- as.numeric(X_grp %*% grp_pos_w)
        group_names_pos <- c(group_names_pos, pos_col_name)
      }

      if (!neg_excluded) {
        neg_col_name <- paste0(safe_name, "_neg")
        group_indices[[neg_col_name]] <- as.numeric(X_grp %*% grp_neg_w)
        group_names_neg <- c(group_names_neg, neg_col_name)
      }
    }
    
    wqs_indices <- group_indices
    wqs_terms <- paste(c(group_names_pos, group_names_neg), collapse = " + ")
  } else {
    wqs_indices <- create_wqs_index(X_quantile, pos_weights, neg_weights)
    wqs_terms <- "wqs_pos + wqs_neg"
  }
  
  glm_data <- wqs_indices
  glm_data$y <- y
  
  if (!is.null(cov_matrix)) {
    cov_df <- as.data.frame(cov_matrix)
    safe_cov_names <- make.names(names(cov_df), unique = TRUE)
    names(cov_df) <- safe_cov_names
    glm_data <- cbind(glm_data, cov_df)
    cov_terms <- paste(safe_cov_names, collapse = " + ")
    formula_str <- paste("y ~", wqs_terms, "+", cov_terms)
  } else {
    formula_str <- paste("y ~", wqs_terms)
  }
  
  fit_formula <- stats::as.formula(formula_str)
  glm_family <- .glm_family_from_string(family, engine = engine)
  
  if (engine == "glm") {
    if (!is.null(obs_weights)) {
      refit_fit <- stats::glm(fit_formula, data = glm_data, family = glm_family, weights = obs_weights)
    } else {
      refit_fit <- stats::glm(fit_formula, data = glm_data, family = glm_family)
    }
  } else {
    if (!requireNamespace("survey", quietly = TRUE)) {
      stop(
        "Package 'survey' is required for `refit_engine = \"svyglm\"`. ",
        "Install it with install.packages('survey').",
        call. = FALSE
      )
    }
    if (is.null(survey_design)) {
      stop("`survey_design` must be provided when `engine = \"svyglm\"`.", call. = FALSE)
    }
    .check_survey_alignment(
      survey_design = survey_design,
      glm_data = glm_data,
      analysis_id = analysis_id,
      context = "refit"
    )
    design_refit <- .inject_design_variables(survey_design, glm_data)
    refit_fit <- survey::svyglm(fit_formula, design = design_refit, family = glm_family)
  }
  
  refit_summary <- summary(refit_fit)
  coef_table <- refit_summary$coefficients
  p_col <- grep("^Pr\\(>\\|[tz]\\|\\)$", colnames(coef_table), value = TRUE)
  if (length(p_col) == 0) {
    p_col <- colnames(coef_table)[ncol(coef_table)]
    warning("Could not identify p-value column; using last column: ", p_col)
  } else {
    p_col <- p_col[1]
  }
  
  group_results <- list()
  if (!is.null(groups) && group_inference) {
    for (grp_name in original_grp_names) {
      safe_name <- grp_name_to_safe[[grp_name]]
      pos_term <- paste0(safe_name, "_pos")
      neg_term <- paste0(safe_name, "_neg")
      pos_excl <- identical(excluded_directions[[grp_name]], "positive")
      neg_excl <- identical(excluded_directions[[grp_name]], "negative")
      pos_exists <- !pos_excl && pos_term %in% rownames(coef_table)
      neg_exists <- !neg_excl && neg_term %in% rownames(coef_table)

      # Always insert — dropped/excluded directions get NA (not omitted)
      group_results[[grp_name]] <- list(
        pos_estimate = if (pos_exists) coef_table[pos_term, "Estimate"] else NA_real_,
        pos_se       = if (pos_exists) coef_table[pos_term, "Std. Error"] else NA_real_,
        pos_pvalue   = if (pos_exists) coef_table[pos_term, p_col] else NA_real_,
        neg_estimate = if (neg_exists) coef_table[neg_term, "Estimate"] else NA_real_,
        neg_se       = if (neg_exists) coef_table[neg_term, "Std. Error"] else NA_real_,
        neg_pvalue   = if (neg_exists) coef_table[neg_term, p_col] else NA_real_,
        pos_excluded = pos_excl,
        neg_excluded = neg_excl
      )
    }
  }
  
  if (!is.null(groups) && group_inference) {
    wqs_pos_estimate <- NA_real_
    wqs_pos_se <- NA_real_
    wqs_pos_pvalue <- NA_real_
    wqs_neg_estimate <- NA_real_
    wqs_neg_se <- NA_real_
    wqs_neg_pvalue <- NA_real_
  } else {
    wqs_pos_exists <- "wqs_pos" %in% rownames(coef_table)
    wqs_neg_exists <- "wqs_neg" %in% rownames(coef_table)
    wqs_pos_estimate <- if (wqs_pos_exists) coef_table["wqs_pos", "Estimate"] else NA_real_
    wqs_pos_se <- if (wqs_pos_exists) coef_table["wqs_pos", "Std. Error"] else NA_real_
    wqs_pos_pvalue <- if (wqs_pos_exists) coef_table["wqs_pos", p_col] else NA_real_
    wqs_neg_estimate <- if (wqs_neg_exists) coef_table["wqs_neg", "Estimate"] else NA_real_
    wqs_neg_se <- if (wqs_neg_exists) coef_table["wqs_neg", "Std. Error"] else NA_real_
    wqs_neg_pvalue <- if (wqs_neg_exists) coef_table["wqs_neg", p_col] else NA_real_
  }
  
  list(
    refit_fit = refit_fit,
    refit_summary = refit_summary,
    coef_table = coef_table,
    wqs_pos_estimate = wqs_pos_estimate,
    wqs_pos_se = wqs_pos_se,
    wqs_pos_pvalue = wqs_pos_pvalue,
    wqs_neg_estimate = wqs_neg_estimate,
    wqs_neg_se = wqs_neg_se,
    wqs_neg_pvalue = wqs_neg_pvalue,
    group_results = group_results,
    group_inference = group_inference && !is.null(groups),
    wqs_indices = wqs_indices,
    p_col = p_col,
    formula = formula_str,
    engine = engine,
    curve = "linear",
    n_obs = n_obs
  )
}


#' Internal Function: Validation Step for Inference
#'
#' Performs the validation-stage GLM using fixed training-estimated weights.
#' Returned p-values are exploratory conditional summaries and do not propagate
#' uncertainty from index estimation or selection.
#'
#' @param X_quantile_val Quantile-transformed validation exposure matrix.
#' @param y_val Validation outcome vector.
#' @param cov_matrix_val Validation covariate matrix.
#' @param pos_weights Fixed positive weights.
#' @param neg_weights Fixed negative weights.
#' @param family Character.
#' @param groups Group list.
#' @param group_inference Logical. Whether to include group-specific indices in GLM.
#'
#' @return A list with validation GLM results.
#'
#' @keywords internal
validation_glm <- function(X_quantile_val, y_val, cov_matrix_val,
                           pos_weights, neg_weights, family, groups,
                           group_inference = TRUE,
                           obs_weights = NULL,
                           excluded_directions = list(),
                           engine = c("glm", "svyglm"),
                           survey_design = NULL,
                           analysis_id = NULL) {
  engine <- match.arg(engine)
  val_result <- refit_model(
    X_quantile = X_quantile_val,
    y = y_val,
    cov_matrix = cov_matrix_val,
    pos_weights = pos_weights,
    neg_weights = neg_weights,
    family = family,
    groups = groups,
    group_inference = group_inference,
    engine = engine,
    survey_design = survey_design,
    analysis_id = analysis_id,
    obs_weights = obs_weights,
    excluded_directions = excluded_directions
  )
  
  val_result$glm_fit <- val_result$refit_fit
  val_result$glm_summary <- val_result$refit_summary
  val_result$n_train <- NA
  val_result$n_val <- length(y_val)
  val_result
}


#' @keywords internal
.validate_analysis_id <- function(analysis_id, n) {
  if (is.null(analysis_id)) {
    return(NULL)
  }
  if (length(analysis_id) != n) {
    stop("`analysis_id` must have length equal to the number of observations.", call. = FALSE)
  }
  analysis_id <- as.character(analysis_id)
  if (any(is.na(analysis_id)) || any(!nzchar(analysis_id))) {
    stop("`analysis_id` must not contain missing or empty values.", call. = FALSE)
  }
  analysis_id
}


#' @keywords internal
.validate_obs_weights <- function(obs_weights, n) {
  if (is.null(obs_weights)) {
    return(NULL)
  }
  if (!is.numeric(obs_weights)) {
    stop("`obs_weights` must be a numeric vector.", call. = FALSE)
  }
  if (length(obs_weights) != n) {
    stop("`obs_weights` must have length equal to the number of observations.", call. = FALSE)
  }
  if (any(!is.finite(obs_weights))) {
    stop("`obs_weights` must contain only finite values.", call. = FALSE)
  }
  if (any(obs_weights < 0)) {
    stop("`obs_weights` must be non-negative.", call. = FALSE)
  }
  if (sum(obs_weights) <= 0) {
    stop("`obs_weights` must have positive total weight.", call. = FALSE)
  }
  obs_weights
}


#' @keywords internal
.validate_survey_design <- function(survey_design) {
  if (is.null(survey_design)) {
    return(invisible(FALSE))
  }
  if (!inherits(survey_design, c("survey.design", "survey.design2", "svyrep.design"))) {
    stop(
      "`survey_design` must be a survey::svydesign() or survey::svrepdesign() object.",
      call. = FALSE
    )
  }
  if (is.null(survey_design$variables) || is.null(nrow(survey_design$variables))) {
    stop("`survey_design` must contain a valid `variables` data frame.", call. = FALSE)
  }
  invisible(TRUE)
}


#' @keywords internal
.has_automatic_rownames <- function(x) {
  rn <- rownames(x)
  is.null(rn) || identical(rn, as.character(seq_len(nrow(x))))
}


#' @keywords internal
.survey_design_ids_for_alignment <- function(survey_design, analysis_id = NULL,
                                             context = "refit") {
  vars <- survey_design$variables
  if (!is.null(vars[[".analysis_id"]])) {
    return(as.character(vars[[".analysis_id"]]))
  }

  if (!is.null(analysis_id)) {
    analysis_id <- as.character(analysis_id)
    matching_cols <- names(vars)[vapply(vars, function(col) {
      identical(as.character(col), analysis_id)
    }, logical(1))]
    if (length(matching_cols) > 0L) {
      return(as.character(vars[[matching_cols[[1L]]]]))
    }

    if (!.has_automatic_rownames(vars)) {
      return(as.character(rownames(vars)))
    }

    stop(
      "`analysis_id` was supplied, but `survey_design` has no `.analysis_id` ",
      "column, no variable matching `analysis_id`, and only automatic row names ",
      "for ", context, " alignment. Add `.analysis_id` to the survey design data ",
      "or give the design stable row names before running survey-aware analysis.",
      call. = FALSE
    )
  }

  if (.has_automatic_rownames(vars)) {
    stop(
      "`survey_design` alignment cannot be verified from automatic row names. ",
      "Supply `analysis_id` and include matching IDs in `survey_design` ",
      "via a `.analysis_id` column, a matching ID column, or stable row names.",
      call. = FALSE
    )
  }

  as.character(rownames(vars))
}


#' @keywords internal
.check_survey_alignment <- function(survey_design, glm_data, analysis_id = NULL,
                                    context = "refit") {
  if (is.null(survey_design)) {
    stop("`survey_design` must be provided when survey-aware refit is requested.", call. = FALSE)
  }
  .validate_survey_design(survey_design)
  if (!identical(nrow(survey_design$variables), nrow(glm_data))) {
    stop(
      "`survey_design` must have the same number of observations as the ", context, " data. ",
      "For NHANES analyses, build `survey_design` from the exact same final analysis dataset ",
      "used for `sglwqs()` after all filtering and row ordering are finalized.",
      call. = FALSE
    )
  }

  if (!is.null(analysis_id)) {
    analysis_id <- as.character(analysis_id)
    if (!.has_automatic_rownames(glm_data) &&
        !identical(analysis_id, as.character(rownames(glm_data)))) {
      stop(
        "`analysis_id` does not match the row names of the ", context, " data. ",
        "When stable row names are present, they must identify the same rows ",
        "in the same order as `analysis_id`.",
        call. = FALSE
      )
    }
    design_ids <- .survey_design_ids_for_alignment(
      survey_design,
      analysis_id = analysis_id,
      context = context
    )
    if (!identical(as.character(design_ids), analysis_id)) {
      stop(
        "`analysis_id` does not match the observation ordering in `survey_design`. ",
        "Build `survey_design` from the exact same final analysis dataset used for analysis.",
        call. = FALSE
      )
    }
    return(invisible(TRUE))
  }

  if (.has_automatic_rownames(glm_data)) {
    stop(
      "`survey_design` alignment cannot be verified because the ", context,
      " data have automatic row names. Supply `analysis_id` or stable row names ",
      "on both the analysis data and `survey_design`.",
      call. = FALSE
    )
  }
  design_row_names <- .survey_design_ids_for_alignment(
    survey_design,
    analysis_id = NULL,
    context = context
  )
  data_row_names <- rownames(glm_data)
  if (!identical(design_row_names, data_row_names)) {
    stop(
      "Row names of `survey_design$variables` do not match the ", context, " data. ",
      "For NHANES analyses, `survey_design` must be built from the exact same final analysis dataset ",
      "in the same row order as the analysis data.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' @keywords internal
.inject_design_variables <- function(survey_design, glm_data) {
  design_refit <- survey_design
  design_refit$variables <- as.data.frame(design_refit$variables)
  for (nm in names(glm_data)) {
    design_refit$variables[[nm]] <- glm_data[[nm]]
  }
  design_refit
}


#' @keywords internal
.subset_survey_design <- function(survey_design, row_idx, context = "analysis") {
  if (is.null(survey_design)) {
    return(NULL)
  }
  .validate_survey_design(survey_design)
  n_design <- nrow(survey_design$variables)
  if (any(is.na(row_idx)) || any(row_idx < 1L) || any(row_idx > n_design)) {
    stop("Cannot subset `survey_design` for ", context, ": row indices are out of range.",
         call. = FALSE)
  }
  keep <- rep(FALSE, n_design)
  keep[row_idx] <- TRUE
  subset(survey_design, keep)
}


#' @keywords internal
.extract_analysis_weights <- function(survey_design) {
  if (is.null(survey_design)) {
    return(NULL)
  }
  weights <- tryCatch(
    stats::weights(survey_design, type = "sampling"),
    error = function(e) {
      tryCatch(stats::weights(survey_design), error = function(e2) NULL)
    }
  )
  if (is.null(weights)) {
    return(NULL)
  }
  as.numeric(weights)
}


#' @keywords internal
.survey_design_degf <- function(design) {
  if (is.null(design) || !requireNamespace("survey", quietly = TRUE)) {
    return(NA_real_)
  }
  tryCatch(as.numeric(survey::degf(design)), error = function(e) NA_real_)
}


#' @keywords internal
.coerce_sparsegl_family <- function(family, obs_weights = NULL) {
  if (is.null(obs_weights)) {
    return(family)
  }
  .glm_family_from_string(family)
}


#' @keywords internal
.glm_family_from_string <- function(family, engine = "glm") {
  if (inherits(family, "family")) {
    return(family)
  }
  if (identical(family, "gaussian")) {
    return(stats::gaussian())
  }
  if (identical(family, "binomial")) {
    if (identical(engine, "svyglm")) {
      return(stats::quasibinomial())
    }
    return(stats::binomial())
  }
  stop("Unsupported family: ", family, call. = FALSE)
}
