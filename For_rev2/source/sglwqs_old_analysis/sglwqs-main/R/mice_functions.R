#' Fit sglwqs with Multiple Imputation (mice)
#'
#' Fits sglwqs models to multiply imputed datasets created by the mice package
#' and pools the results using Rubin's rules.
#'
#' @param mids_obj A mids object from mice::mice().
#' @param formula A formula specifying the model: y ~ x1 + x2 + ... | cov1 + cov2 + ...
#'   Variables before | are exposures, variables after | are covariates.
#'   Alternatively, use exposure_vars and covariate_vars arguments.
#' @param exposure_vars Character vector of exposure variable names.
#' @param outcome_var Character. Name of the outcome variable.
#' @param covariate_vars Character vector of covariate variable names (optional).
#' @param X Optional alias for \code{exposure_vars}. In the MI interface this should
#'   identify exposure columns inside \code{mids_obj} by name, index, logical selector,
#'   or a named matrix/data.frame.
#' @param y Optional alias for \code{outcome_var}. Must identify exactly one outcome
#'   column inside \code{mids_obj}.
#' @param covariates Optional alias for \code{covariate_vars}.
#' @param groups Named list specifying variable groups for Group Lasso penalty.
#' @param n_quantiles Integer. Number of quantiles (default: 4).
#' @param family Character. "gaussian" or "binomial" (default: "gaussian").
#' @param lambda Character or numeric. Lambda selection method (default: "lambda.min").
#' @param nfolds Integer. Number of CV folds (default: 10).
#' @param penalize_covariates Logical. Whether covariates should also be
#'   regularized (default: FALSE). Set to FALSE when adjustment covariates should
#'   always remain in the training model.
#' @param group_by_compound Logical. Passed to \code{sglwqs()}. Defaults to the
#'   same behavior as \code{sglwqs()}.
#' @param group_structure Character. Currently only \code{"direction"} is supported.
#' @param bootstrap Logical. Use bootstrap aggregation (default: FALSE).
#' @param n_boot Integer. Number of bootstrap iterations (default: 100).
#' @param parallel Logical. Use parallel processing for bootstrap (default: FALSE).
#' @param keep_boot_matrices Logical. Passed to \code{sglwqs()}; keep per-bootstrap
#'   coefficient matrices inside each imputation fit (default: FALSE).
#' @param checkpoint_dir Character. Directory to save intermediate bootstrap results 
#'   (default: NULL). Each imputation's checkpoints are saved in a subdirectory.
#' @param checkpoint_interval Integer. Bootstrap iterations between checkpoints (default: 50).
#' @param cleanup_checkpoint Logical. Delete checkpoints after completion (default: TRUE).
#' @param future_globals_max_size Numeric. Maximum size (in bytes) of global variables 
#'   for parallel workers (default: NULL, auto-detect based on data size).
#' @param stratified_bootstrap Logical. Whether to use stratified bootstrap for binomial family 
#'   (default: TRUE). When TRUE and family = "binomial", bootstrap samples maintain the original 
#'   case/control ratio, which improves convergence stability for imbalanced binary outcomes.
#' @param validation Logical. Use train/validation split (default: FALSE).
#' @param train_prop Numeric. Proportion for training (default: 0.6).
#' @param refit Character. One of \code{"none"}, \code{"full"}, or
#'   \code{"validation"}. Passed to \code{sglwqs()}. Default: \code{"none"}.
#' @param seed Integer. Random seed (default: NULL).
#' @param verbose Logical. Show progress (default: TRUE).
#' @param minor_threshold Numeric (0-1). Passed to \code{sglwqs()}.
#'   Group-direction indices whose coefficient mass ratio (minor/major) falls
#'   below this threshold are excluded from the validation GLM to prevent
#'   suppressor effects. Default: 0.10.
#' @param ... Additional arguments passed to sglwqs().
#'
#' @return An object of class \code{"sglwqs_mids"} containing:
#' \itemize{
#'   \item \code{fits}: List of \code{sglwqs} objects for each imputed dataset.
#'   \item \code{pooled}: Pooled results across imputations. Important
#'     sub-components include \code{pooled$inference}, \code{pooled$covariates},
#'     \code{pooled$bootstrap}, and \code{pooled$bootstrap_inference}.
#'   \item \code{m}: Number of imputations.
#'   \item \code{diagnostics}: Fact-only summaries returned by
#'     \code{compute_diagnostics()}.
#'   \item \code{call}: The function call.
#' }
#'
#' @details
#' This function fits sglwqs to each imputed dataset and then pools the results
#' using Rubin's rules for valid inference with multiply imputed data.
#'
#' \strong{Pooling Strategy:}
#' \itemize{
#'   \item Weights: Averaged across imputations
#'   \item Selection frequency: Proportion selected across all imputations
#'   \item Downstream GLM estimates (including covariates): Pooled term-wise using Rubin's rules
#'   \item Bootstrap-only summaries: Combined across imputations as
#'     \code{pooled$bootstrap_inference} when \code{bootstrap = TRUE} and no
#'     downstream GLM is requested
#' }
#'
#' WQS composite-index coefficients and covariate coefficients are Rubin-pooled
#' whenever each imputation fit includes a downstream GLM, either via
#' \code{validation = TRUE} or \code{refit = "full"}. When both information
#' sources are present, \code{refit_info} is preferred over
#' \code{validation_info} because it uses the full analysis sample.
#'
#' Access pooled downstream inference via:
#' \itemize{
#'   \item \code{fit$pooled$inference$wqs_pos}
#'   \item \code{fit$pooled$inference$group_results}
#'   \item \code{fit$pooled$covariates}
#' }
#'
#' When \code{bootstrap = TRUE} but no downstream GLM is requested,
#' \code{fit$pooled$bootstrap_inference} stores Rubin-pooled bootstrap summaries
#' for WQS index-sum estimates and covariates.
#'
#' For backward compatibility, \code{fit$pooled$validation} remains populated
#' when validation-based inference is used.
#'
#' @examples
#' \dontrun{
#' library(mice)
#' library(sglwqs)
#'
#' # Create example data with missing values
#' data <- data.frame(
#'   y = rnorm(100),
#'   x1 = rnorm(100), x2 = rnorm(100), x3 = rnorm(100),
#'   cov1 = rnorm(100)
#' )
#' data[sample(100, 10), "x1"] <- NA
#' data[sample(100, 10), "x2"] <- NA
#'
#' # Multiple imputation
#' imp <- mice(data, m = 5, printFlag = FALSE)
#'
#' # Fit sglwqs with multiple imputation
#' fit <- sglwqs_mice(
#'   imp,
#'   X = c("x1", "x2", "x3"),
#'   y = "y",
#'   covariates = "cov1",
#'   validation = TRUE,
#'   penalize_covariates = FALSE
#' )
#'
#' summary(fit)
#' }
#'
#' @export
sglwqs_mice <- function(mids_obj,
                         formula = NULL,
                         exposure_vars = NULL,
                         outcome_var = NULL,
                         covariate_vars = NULL,
                         X = NULL,
                         y = NULL,
                         covariates = NULL,
                         groups = NULL,
                         n_quantiles = 4,
                         family = c("gaussian", "binomial"),
                         lambda = "lambda.min",
                         nfolds = 10,
                         penalize_covariates = FALSE,
                         group_by_compound = NULL,
                         group_structure = "direction",
                         bootstrap = FALSE,
                         n_boot = 100,
                         parallel = FALSE,
                         keep_boot_matrices = FALSE,
                         checkpoint_dir = NULL,
                         checkpoint_interval = 50,
                         cleanup_checkpoint = TRUE,
                         future_globals_max_size = NULL,
                         stratified_bootstrap = TRUE,
                         validation = FALSE,
                         train_prop = 0.6,
                         refit = c("none", "full", "validation"),
                         seed = NULL,
                         verbose = TRUE,
                         minor_threshold = 0.10,
                         ...) {
  refit_missing <- missing(refit)

  if (!requireNamespace("mice", quietly = TRUE)) {
    stop("Package 'mice' is required. Install with: install.packages('mice')")
  }

  if (!inherits(mids_obj, "mids")) {
    stop("mids_obj must be a 'mids' object from mice::mice()")
  }

  family <- match.arg(family)
  refit <- match.arg(refit)
  if (validation && !refit_missing && !identical(refit, "validation")) {
    stop(
      "`validation = TRUE` is only compatible with `refit = \"validation\"`.",
      call. = FALSE
    )
  }
  if (validation) {
    refit <- "validation"
  }
  if (!identical(group_structure, "direction")) {
    warning("group_structure = '", group_structure, "' is no longer supported; using 'direction'.",
            call. = FALSE)
    group_structure <- "direction"
  }

  m <- mids_obj$m
  template_data <- mice::complete(mids_obj, action = 1)
  data_names <- names(template_data)

  if (verbose) {
    message("Fitting sglwqs to ", m, " imputed datasets...")
  }

  if (!is.null(formula)) {
    formula_parts <- parse_sglwqs_formula(formula)
    if (is.null(exposure_vars)) exposure_vars <- formula_parts$exposures
    if (is.null(outcome_var)) outcome_var <- formula_parts$outcome
    if (is.null(covariate_vars)) covariate_vars <- formula_parts$covariates
  }

  if (!is.null(X)) {
    X_alias <- .resolve_mids_var_spec(X, data_names, "X")
    if (is.null(exposure_vars)) {
      exposure_vars <- X_alias
    } else if (!identical(.resolve_mids_var_spec(exposure_vars, data_names, "exposure_vars"), X_alias)) {
      stop("`X` and `exposure_vars` refer to different variables.")
    }
  }

  if (!is.null(y)) {
    y_alias <- .resolve_single_mids_var(y, data_names, "y")
    if (is.null(outcome_var)) {
      outcome_var <- y_alias
    } else if (!identical(.resolve_single_mids_var(outcome_var, data_names, "outcome_var"), y_alias)) {
      stop("`y` and `outcome_var` refer to different variables.")
    }
  }

  if (!is.null(covariates)) {
    cov_alias <- .resolve_mids_var_spec(covariates, data_names, "covariates")
    if (is.null(covariate_vars)) {
      covariate_vars <- cov_alias
    } else if (!identical(.resolve_mids_var_spec(covariate_vars, data_names, "covariate_vars"), cov_alias)) {
      stop("`covariates` and `covariate_vars` refer to different variables.")
    }
  }

  if (!is.null(exposure_vars)) {
    exposure_vars <- .resolve_mids_var_spec(exposure_vars, data_names, "exposure_vars")
  }
  if (!is.null(outcome_var)) {
    outcome_var <- .resolve_single_mids_var(outcome_var, data_names, "outcome_var")
  }
  if (!is.null(covariate_vars)) {
    covariate_vars <- .resolve_mids_var_spec(covariate_vars, data_names, "covariate_vars")
  }

  if (is.null(exposure_vars) || is.null(outcome_var)) {
    stop("Specify exposures/outcome via `exposure_vars` + `outcome_var`, `X` + `y`, or `formula`.")
  }

  # ----- Resolve group spec once (same rules as sglwqs()) -----
  if (is.null(group_by_compound)) {
    group_by_compound <- !is.null(groups)
  }
  if (!isTRUE(group_by_compound) && !is.null(groups)) {
    warning("`groups` is ignored when `group_by_compound = FALSE`.", call. = FALSE)
    groups <- NULL
  } else {
    groups <- .validate_and_complete_groups(exposure_vars, groups)
  }
  # From here on, `groups` is the canonical, validated version (or NULL).

  if (!is.null(seed)) set.seed(seed)
  dots <- list(...)

  .advise_on_settings(
    bootstrap = bootstrap,
    validation = validation,
    refit = if (validation) "validation" else refit,
    is_mids = TRUE,
    verbose = verbose
  )

  fits <- vector("list", m)

  for (i in seq_len(m)) {
    if (verbose) {
      message("  Imputation ", i, "/", m, "...")
    }

    complete_data <- mice::complete(mids_obj, action = i)

    X_data <- complete_data[, exposure_vars, drop = FALSE]
    y_data <- complete_data[[outcome_var]]

    if (!is.null(covariate_vars)) {
      covariate_data <- complete_data[, covariate_vars, drop = FALSE]
    } else {
      covariate_data <- NULL
    }

    imp_checkpoint_dir <- NULL
    if (!is.null(checkpoint_dir)) {
      imp_checkpoint_dir <- file.path(checkpoint_dir, paste0("imp_", i))
    }

    fits[[i]] <- tryCatch({
      sglwqs(
        X = X_data,
        y = y_data,
        covariates = covariate_data,
        groups = groups,
        n_quantiles = n_quantiles,
        family = family,
        lambda = lambda,
        nfolds = nfolds,
        penalize_covariates = penalize_covariates,
        group_by_compound = group_by_compound,
        group_structure = group_structure,
        bootstrap = bootstrap,
        n_boot = n_boot,
        parallel = parallel,
        keep_boot_matrices = keep_boot_matrices,
        checkpoint_dir = imp_checkpoint_dir,
        checkpoint_interval = checkpoint_interval,
        cleanup_checkpoint = cleanup_checkpoint,
        future_globals_max_size = future_globals_max_size,
        stratified_bootstrap = stratified_bootstrap,
        validation = validation,
        train_prop = train_prop,
        refit = refit,
        seed = if (!is.null(seed)) seed + i * 10000L else NULL,
        verbose = FALSE,
        minor_threshold = minor_threshold,
        ...
      )
    }, error = function(e) {
      warning("Imputation ", i, " failed: ", e$message)
      return(NULL)
    })
  }

  successful_fits <- Filter(Negate(is.null), fits)
  n_successful <- length(successful_fits)

  if (n_successful == 0) {
    stop("All imputations failed")
  }

  if (n_successful < m) {
    warning(m - n_successful, " imputations failed. Pooling ", n_successful, " successful fits.")
  }

  if (verbose) {
    message("Pooling results from ", n_successful, " imputations...")
  }

  pooled <- pool_sglwqs_results(successful_fits, verbose = verbose)

  result <- list(
    fits = fits,
    pooled = pooled,
    m = m,
    n_successful = n_successful,
    exposure_vars = exposure_vars,
    outcome_var = outcome_var,
    covariate_vars = covariate_vars,
    groups = groups,
    family = family,
    bootstrap = bootstrap,
    validation = validation,
    minor_threshold = minor_threshold,
    penalize_covariates = penalize_covariates,
    group_by_compound = group_by_compound,
    keep_boot_matrices = keep_boot_matrices,
    call = match.call()
  )

  class(result) <- "sglwqs_mids"
  result$diagnostics <- compute_diagnostics(result)
  .emit_diagnostics(result$diagnostics, verbose = verbose)

  return(result)
}

#' @keywords internal
.resolve_mids_var_spec <- function(spec, data_names, arg_name) {
  if (is.null(spec)) {
    return(NULL)
  }

  if (is.character(spec)) {
    vars <- spec
  } else if (is.numeric(spec)) {
    if (any(!is.finite(spec)) || any(spec < 1) || any(spec > length(data_names)) ||
        any(spec != as.integer(spec))) {
      stop(arg_name, " must contain valid column indices for `mids_obj`.")
    }
    vars <- data_names[as.integer(spec)]
  } else if (is.logical(spec)) {
    if (length(spec) != length(data_names)) {
      stop(arg_name, " must have length ", length(data_names), " when supplied as a logical selector.")
    }
    vars <- data_names[spec]
  } else if (is.data.frame(spec) || is.matrix(spec)) {
    if (is.null(colnames(spec))) {
      stop(arg_name, " must have column names when supplied as a matrix/data.frame.")
    }
    vars <- colnames(spec)
  } else {
    stop(arg_name, " must be specified using variable names, numeric indices, a logical selector, or a named matrix/data.frame.")
  }

  vars <- unique(as.character(vars))
  missing_vars <- setdiff(vars, data_names)
  if (length(missing_vars) > 0) {
    stop(arg_name, " contains variables not found in `mids_obj`: ", paste(missing_vars, collapse = ", "))
  }

  vars
}

#' @keywords internal
.resolve_single_mids_var <- function(spec, data_names, arg_name) {
  vars <- .resolve_mids_var_spec(spec, data_names, arg_name)
  if (length(vars) != 1) {
    stop(arg_name, " must identify exactly one variable.")
  }
  vars
}


#' Parse sglwqs formula
#'
#' @param formula A formula: y ~ x1 + x2 | cov1 + cov2
#' @return A list with outcome, exposures, and covariates
#' @keywords internal
parse_sglwqs_formula <- function(formula) {
  # Parse y ~ x1 + x2 + x3 | cov1 + cov2 format
  formula_str <- paste(deparse(formula, width.cutoff = 500L), collapse = "")
  
  # Split by |
  if (grepl("\\|", formula_str)) {
    parts <- strsplit(formula_str, "\\|")[[1]]
    main_part <- trimws(parts[1])
    cov_part <- trimws(parts[2])
    
    # Extract covariates
    covariates <- trimws(strsplit(cov_part, "\\+")[[1]])
  } else {
    main_part <- formula_str
    covariates <- NULL
  }
  
  # Extract from y ~ x1 + x2
  main_formula <- as.formula(main_part)
  outcome <- as.character(main_formula[[2]])
  
  # Extract variables from right-hand side
  rhs <- as.character(main_formula[[3]])
  if (length(rhs) == 1) {
    exposures <- rhs
  } else {
    exposures <- all.vars(main_formula[[3]])
  }
  
  return(list(
    outcome = outcome,
    exposures = exposures,
    covariates = covariates
  ))
}


#' Normalize group metadata for equality comparisons
#'
#' Sorts group names alphabetically and sorts each group's member variables
#' so that equivalent group specs compare identically even when names or
#' members arrive in different order.
#'
#' @param groups A named list of character vectors, or \code{NULL}.
#' @return A normalized named list, or \code{NULL}.
#' @keywords internal
.normalize_group_metadata <- function(groups) {
  if (is.null(groups)) {
    return(NULL)
  }

  group_names <- sort(names(groups))
  out <- groups[group_names]
  out <- lapply(out, function(vars) sort(unique(as.character(vars))))
  out
}

.row_var_or_zero <- function(x) {
  x <- as.matrix(x)
  if (ncol(x) <= 1L) {
    return(rep(0, nrow(x)))
  }
  apply(x, 1, var)
}

#' Pool sglwqs Results Using Rubin's Rules
#'
#' @param fits List of sglwqs objects
#' @param verbose Logical
#' @return A list with pooled results
#' @keywords internal
pool_sglwqs_results <- function(fits, verbose = TRUE) {

  m <- length(fits)
  first_fit <- fits[[1]]
  var_names <- first_fit$var_names
  p <- first_fit$p

  # ----- 0. Validate metadata consistency across fits -----
  ref_groups <- .normalize_group_metadata(first_fit$groups)
  for (i in seq_along(fits)[-1]) {
    if (!identical(sort(fits[[i]]$var_names), sort(var_names))) {
      stop("Variable names differ between imputation 1 and imputation ", i, ".",
           call. = FALSE)
    }
    fi_groups <- .normalize_group_metadata(fits[[i]]$groups)
    if (!identical(ref_groups, fi_groups)) {
      stop("Group metadata differs between imputation 1 (",
           paste(names(ref_groups), collapse = ", "), ") and imputation ", i, " (",
           paste(names(fi_groups), collapse = ", "), ").",
           call. = FALSE)
    }
  }

  # ----- 1. Pool weights (simple mean) -----
  # Align by name to be robust against ordering differences
  pos_weights_mat <- sapply(fits, function(f) f$pos_weights[var_names])
  neg_weights_mat <- sapply(fits, function(f) f$neg_weights[var_names])
  
  pooled_pos_weights <- rowMeans(pos_weights_mat)
  pooled_neg_weights <- rowMeans(neg_weights_mat)
  names(pooled_pos_weights) <- var_names
  names(pooled_neg_weights) <- var_names
  
  # Weight variation (between imputations)
  pos_weights_var <- .row_var_or_zero(pos_weights_mat)
  neg_weights_var <- .row_var_or_zero(neg_weights_mat)
  
  # ----- 2. Pool coefficients -----
  pos_coef_mat <- sapply(fits, function(f) f$pos_coef[var_names])
  neg_coef_mat <- sapply(fits, function(f) f$neg_coef[var_names])
  
  pooled_pos_coef <- rowMeans(pos_coef_mat)
  pooled_neg_coef <- rowMeans(neg_coef_mat)
  names(pooled_pos_coef) <- var_names
  names(pooled_neg_coef) <- var_names
  
  # ----- 3. Pool index sums -----
  pos_index_sums <- sapply(fits, function(f) f$pos_index_sum)
  neg_index_sums <- sapply(fits, function(f) f$neg_index_sum)
  
  pooled_pos_index_sum <- mean(pos_index_sums)
  pooled_neg_index_sum <- mean(neg_index_sums)
  
  # ----- 4. Selection frequency (proportion non-zero across imputations) -----
  mi_selection_freq_pos <- rowMeans(pos_coef_mat > 0)
  mi_selection_freq_neg <- rowMeans(neg_coef_mat > 0)
  names(mi_selection_freq_pos) <- var_names
  names(mi_selection_freq_neg) <- var_names
  
  # ----- 5. Pool bootstrap info (if available) -----
  if (first_fit$bootstrap) {
    # Average selection frequency across imputations (align by name)
    boot_sel_freq_pos <- sapply(fits, function(f) f$boot_info$selection_freq_pos[var_names])
    boot_sel_freq_neg <- sapply(fits, function(f) f$boot_info$selection_freq_neg[var_names])
    
    # Weighted mean by n_successful bootstrap draws per imputation
    boot_n_success <- vapply(fits, function(f) {
      ns <- f$boot_info$n_successful
      if (is.null(ns) || !is.finite(ns) || ns < 0) 0 else ns
    }, numeric(1))
    boot_total <- sum(boot_n_success)
    if (boot_total > 0) {
      pooled_boot_sel_freq_pos <- as.numeric((boot_sel_freq_pos %*% boot_n_success) / boot_total)
      pooled_boot_sel_freq_neg <- as.numeric((boot_sel_freq_neg %*% boot_n_success) / boot_total)
    } else {
      pooled_boot_sel_freq_pos <- rowMeans(boot_sel_freq_pos)
      pooled_boot_sel_freq_neg <- rowMeans(boot_sel_freq_neg)
    }
    names(pooled_boot_sel_freq_pos) <- var_names
    names(pooled_boot_sel_freq_neg) <- var_names
    
    # Pool SE info (Rubin's rules: T = W_bar + (1 + 1/m) * B)
    boot_se_pos <- sapply(fits, function(f) f$boot_info$se_pos_coef[var_names])
    boot_se_neg <- sapply(fits, function(f) f$boot_info$se_neg_coef[var_names])

    # Within-imputation variance: mean of SE^2
    W_bar_pos <- rowMeans(boot_se_pos^2)
    W_bar_neg <- rowMeans(boot_se_neg^2)

    # Between-imputation variance: variance of point estimates
    boot_mean_pos <- sapply(fits, function(f) f$boot_info$mean_pos_coef[var_names])
    boot_mean_neg <- sapply(fits, function(f) f$boot_info$mean_neg_coef[var_names])
    B_pos <- .row_var_or_zero(boot_mean_pos)
    B_neg <- .row_var_or_zero(boot_mean_neg)

    # Total variance (Rubin's rules)
    pooled_boot_se_pos <- sqrt(W_bar_pos + (1 + 1 / m) * B_pos)
    pooled_boot_se_neg <- sqrt(W_bar_neg + (1 + 1 / m) * B_neg)
    
    pooled_bootstrap <- list(
      selection_freq_pos = pooled_boot_sel_freq_pos,
      selection_freq_neg = pooled_boot_sel_freq_neg,
      se_pos_coef = pooled_boot_se_pos,
      se_neg_coef = pooled_boot_se_neg
    )
  } else {
    pooled_bootstrap <- NULL
  }
  
  # ----- 6. Pool downstream inference (Rubin's rules) -----
  has_refit <- !is.null(first_fit$refit_info)
  has_validation <- isTRUE(first_fit$validation)
  if (has_refit || has_validation) {
    pooled_inference <- pool_refit_rubin(fits, source = "auto", verbose = verbose)
    pooled_covariates <- pool_covariate_rubin(fits, source = "auto", verbose = verbose)
  } else {
    pooled_inference <- NULL
    pooled_covariates <- NULL
  }

  # ----- 6b. Pool bootstrap-based second-stage approximations -----
  if (first_fit$bootstrap) {
    pooled_bootstrap_inference <- pool_bootstrap_inference(fits, verbose = verbose)
  } else {
    pooled_bootstrap_inference <- NULL
  }
  
  # ----- 7. Pool group info -----
  if (!is.null(first_fit$groups)) {
    # Pool index sums by group
    group_names <- names(first_fit$groups)
    
    pooled_pos_index_by_group <- list()
    pooled_neg_index_by_group <- list()
    
    for (grp in group_names) {
      pos_grp_sums <- sapply(fits, function(f) f$pos_index_sum_by_group[[grp]])
      neg_grp_sums <- sapply(fits, function(f) f$neg_index_sum_by_group[[grp]])
      
      pooled_pos_index_by_group[[grp]] <- mean(pos_grp_sums)
      pooled_neg_index_by_group[[grp]] <- mean(neg_grp_sums)
    }
  } else {
    pooled_pos_index_by_group <- NULL
    pooled_neg_index_by_group <- NULL
  }
  
  # Compile results
  pooled <- list(
    # Weights
    pos_weights = pooled_pos_weights,
    neg_weights = pooled_neg_weights,
    pos_weights_var = pos_weights_var,
    neg_weights_var = neg_weights_var,
    
    # Coefficients
    pos_coef = pooled_pos_coef,
    neg_coef = pooled_neg_coef,
    
    # Index sums
    pos_index_sum = pooled_pos_index_sum,
    neg_index_sum = pooled_neg_index_sum,
    pos_index_sum_by_group = pooled_pos_index_by_group,
    neg_index_sum_by_group = pooled_neg_index_by_group,
    
    # Selection frequency (between imputations)
    mi_selection_freq_pos = mi_selection_freq_pos,
    mi_selection_freq_neg = mi_selection_freq_neg,
    
    # Bootstrap info
    bootstrap = pooled_bootstrap,
    bootstrap_inference = pooled_bootstrap_inference,
    
    # Downstream inference (Rubin's rules)
    inference = pooled_inference,
    covariates = pooled_covariates,
    validation = if (!is.null(pooled_inference) &&
                     identical(pooled_inference$source_used, "validation_info")) {
      pooled_inference
    } else {
      NULL
    },
    
    # Meta information
    m = m,
    var_names = var_names,
    groups = first_fit$groups
  )
  
  return(pooled)
}


#' Single source of truth for groups in a sglwqs_mids object
#'
#' Priority: top-level groups > pooled groups > first successful fit.
#' All MICE-side S3 methods should use this instead of accessing
#' \code{object$groups} or \code{object$pooled$groups} directly.
#'
#' @param object A \code{sglwqs_mids} object.
#' @return A named list of groups, or \code{NULL}.
#' @keywords internal
.get_effective_groups <- function(object) {
  if (isTRUE(object$group_by_compound == FALSE)) return(NULL)
  if (!is.null(object$groups)) return(object$groups)

  if (!is.null(object$pooled$groups)) return(object$pooled$groups)

  fits <- Filter(function(f) inherits(f, "sglwqs"), object$fits)
  if (length(fits) > 0) return(fits[[1]]$groups)

  NULL
}

#' @keywords internal
.resolve_mi_inference_source <- function(fit,
                                         source = c("auto", "refit_info", "validation_info")) {
  source <- match.arg(source)

  if (identical(source, "refit_info")) {
    return(if (!is.null(fit$refit_info)) "refit_info" else NULL)
  }
  if (identical(source, "validation_info")) {
    return(if (!is.null(fit$validation_info)) "validation_info" else NULL)
  }

  if (isTRUE(fit$validation) && !is.null(fit$validation_info)) {
    return("validation_info")
  }
  if (!is.null(fit$refit_info)) {
    return("refit_info")
  }
  if (!is.null(fit$validation_info)) {
    return("validation_info")
  }
  NULL
}

#' @keywords internal
.resolve_mi_inference_info <- function(fit,
                                       source = c("auto", "refit_info", "validation_info")) {
  source_used <- .resolve_mi_inference_source(fit, source = source)
  if (is.null(source_used)) {
    return(NULL)
  }
  fit[[source_used]]
}

#' @keywords internal
.get_mi_pooled_inference <- function(object) {
  object$pooled$inference %||% object$pooled$validation
}


#' Pool Refit/Validation Results Using Rubin's Rules
#'
#' @param fits List of \code{sglwqs} objects with downstream GLM results.
#' @param source Character. One of \code{"auto"}, \code{"refit_info"}, or
#'   \code{"validation_info"}. \code{"auto"} prefers \code{refit_info} over
#'   \code{validation_info} when available.
#' @param verbose Logical. Whether to print pooling progress messages.
#' @return A list with pooled downstream-inference statistics.
#' @keywords internal
pool_refit_rubin <- function(fits,
                             source = c("auto", "refit_info", "validation_info"),
                             verbose = TRUE) {
  source <- match.arg(source)
  m <- length(fits)
  first_info <- .resolve_mi_inference_info(fits[[1]], source = source)

  if (is.null(first_info)) {
    if (verbose) {
      message("  No refit_info or validation_info available; skipping pooled inference.")
    }
    return(NULL)
  }

  source_used <- .resolve_mi_inference_source(fits[[1]], source = source)
  has_group_results <- !is.null(first_info$group_results) &&
    length(first_info$group_results) > 0

  if (verbose) {
    message(
      "  Pooling inference from `", source_used, "`: ",
      ifelse(has_group_results, "group-level inference", "overall inference")
    )
  }

  if (has_group_results) {
    group_names <- unique(unlist(lapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      if (!is.null(f$groups)) {
        names(f$groups)
      } else if (!is.null(info$group_results)) {
        names(info$group_results)
      } else {
        NULL
      }
    })))
    group_results <- list()

    for (grp in group_names) {
      pos_estimates <- sapply(fits, function(f) {
        info <- .resolve_mi_inference_info(f, source = source_used)
        res <- info$group_results[[grp]]
        if (is.null(res) || is.na(res$pos_estimate)) NA_real_ else res$pos_estimate
      })
      pos_ses <- sapply(fits, function(f) {
        info <- .resolve_mi_inference_info(f, source = source_used)
        res <- info$group_results[[grp]]
        if (is.null(res) || is.na(res$pos_se)) NA_real_ else res$pos_se
      })

      neg_estimates <- sapply(fits, function(f) {
        info <- .resolve_mi_inference_info(f, source = source_used)
        res <- info$group_results[[grp]]
        if (is.null(res) || is.na(res$neg_estimate)) NA_real_ else res$neg_estimate
      })
      neg_ses <- sapply(fits, function(f) {
        info <- .resolve_mi_inference_info(f, source = source_used)
        res <- info$group_results[[grp]]
        if (is.null(res) || is.na(res$neg_se)) NA_real_ else res$neg_se
      })

      pos_valid <- !is.na(pos_estimates) & !is.na(pos_ses) & pos_ses > 0
      neg_valid <- !is.na(neg_estimates) & !is.na(neg_ses) & neg_ses > 0

      if (verbose) {
        message("    ", grp, ": pos valid=", sum(pos_valid), "/", m,
                ", neg valid=", sum(neg_valid), "/", m)
      }

      pos_pooled <- if (sum(pos_valid) >= 1) {
        rubin_pool(
          pos_estimates[pos_valid],
          pos_ses[pos_valid]^2,
          sum(pos_valid),
          df_complete = .mi_complete_df(fits, pos_valid)
        )
      } else {
        list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
             p_value = NA, fmi = NA, df = NA, t_stat = NA)
      }

      neg_pooled <- if (sum(neg_valid) >= 1) {
        rubin_pool(
          neg_estimates[neg_valid],
          neg_ses[neg_valid]^2,
          sum(neg_valid),
          df_complete = .mi_complete_df(fits, neg_valid)
        )
      } else {
        list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
             p_value = NA, fmi = NA, df = NA, t_stat = NA)
      }

      group_results[[grp]] <- list(
        positive = pos_pooled,
        negative = neg_pooled
      )
    }

    result <- list(
      group_results = group_results,
      wqs_pos = list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
                     p_value = NA, fmi = NA, df = NA, t_stat = NA),
      wqs_neg = list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
                     p_value = NA, fmi = NA, df = NA, t_stat = NA),
      has_group_results = TRUE
    )
  } else {
    pos_estimates <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$wqs_pos_estimate
    })
    pos_ses <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$wqs_pos_se
    })
    pos_valid <- !is.na(pos_estimates) & !is.na(pos_ses) & pos_ses > 0

    if (verbose) {
      message("    Overall: pos valid=", sum(pos_valid), "/", m)
    }

    pos_pooled <- if (sum(pos_valid) >= 1) {
      rubin_pool(
        pos_estimates[pos_valid],
        pos_ses[pos_valid]^2,
        sum(pos_valid),
        df_complete = .mi_complete_df(fits, pos_valid)
      )
    } else {
      list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
           p_value = NA, fmi = NA, df = NA, t_stat = NA)
    }

    neg_estimates <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$wqs_neg_estimate
    })
    neg_ses <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$wqs_neg_se
    })
    neg_valid <- !is.na(neg_estimates) & !is.na(neg_ses) & neg_ses > 0

    if (verbose) {
      message("    Overall: neg valid=", sum(neg_valid), "/", m)
    }

    neg_pooled <- if (sum(neg_valid) >= 1) {
      rubin_pool(
        neg_estimates[neg_valid],
        neg_ses[neg_valid]^2,
        sum(neg_valid),
        df_complete = .mi_complete_df(fits, neg_valid)
      )
    } else {
      list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
           p_value = NA, fmi = NA, df = NA, t_stat = NA)
    }

    result <- list(
      wqs_pos = pos_pooled,
      wqs_neg = neg_pooled,
      group_results = NULL,
      has_group_results = FALSE
    )
  }

  if (identical(source_used, "validation_info")) {
    result$n_train <- mean(sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$n_train
    }))
    result$n_val <- mean(sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$n_val
    }))
    result$n_obs <- NA_real_
  } else {
    result$n_train <- NA_real_
    result$n_val <- NA_real_
    result$n_obs <- mean(sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      info$n_obs
    }))
  }

  result$source_used <- source_used
  result
}

#' Pool Validation Results Using Rubin's Rules
#' @keywords internal
pool_validation_rubin <- function(fits, verbose = TRUE) {
  pool_refit_rubin(fits, source = "validation_info", verbose = verbose)
}

#' Pool covariate coefficients from downstream GLM coef_table via Rubin's rules
#' @keywords internal
pool_covariate_rubin <- function(fits,
                                 source = c("auto", "refit_info", "validation_info"),
                                 verbose = TRUE) {
  source <- match.arg(source)
  source_used <- .resolve_mi_inference_source(fits[[1]], source = source)
  if (is.null(source_used)) {
    return(NULL)
  }

  cov_names <- fits[[1]]$cov_names
  if (is.null(cov_names) || length(cov_names) == 0) {
    return(NULL)
  }

  results <- list()
  for (cov_nm in cov_names) {
    estimates <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      if (is.null(info) || is.null(info$coef_table)) return(NA_real_)
      if (cov_nm %in% rownames(info$coef_table)) info$coef_table[cov_nm, "Estimate"] else NA_real_
    })
    ses <- sapply(fits, function(f) {
      info <- .resolve_mi_inference_info(f, source = source_used)
      if (is.null(info) || is.null(info$coef_table)) return(NA_real_)
      if (cov_nm %in% rownames(info$coef_table)) info$coef_table[cov_nm, "Std. Error"] else NA_real_
    })

    valid <- !is.na(estimates) & !is.na(ses) & ses > 0
    if (sum(valid) >= 1) {
      results[[cov_nm]] <- rubin_pool(
        estimates[valid],
        ses[valid]^2,
        sum(valid),
        df_complete = .mi_complete_df(fits, valid)
      )
    } else {
      results[[cov_nm]] <- list(
        estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
        p_value = NA, fmi = NA, df = NA, t_stat = NA
      )
    }
  }

  if (verbose && length(results) > 0) {
    message("  Pooled ", length(results), " covariate coefficient(s) from `", source_used, "`")
  }

  results
}


#' Pool bootstrap-derived index-sum and covariate summaries via Rubin's rules
#' @keywords internal
pool_bootstrap_inference <- function(fits, verbose = TRUE) {
  if (length(fits) == 0 || !any(vapply(fits, function(f) !is.null(f$boot_info), logical(1)))) {
    return(NULL)
  }

  first_fit <- fits[[1]]
  has_groups <- !is.null(first_fit$groups)
  result <- list(
    source_used = "boot_info",
    has_group_results = has_groups
  )

  if (has_groups) {
    group_results <- list()
    for (grp in names(first_fit$groups)) {
      pos_est <- vapply(fits, function(f) {
        f$boot_info$mean_index_sum_by_group_pos[[grp]] %||% NA_real_
      }, numeric(1))
      pos_se <- vapply(fits, function(f) {
        f$boot_info$se_index_sum_by_group_pos[[grp]] %||% NA_real_
      }, numeric(1))
      neg_est <- vapply(fits, function(f) {
        f$boot_info$mean_index_sum_by_group_neg[[grp]] %||% NA_real_
      }, numeric(1))
      neg_se <- vapply(fits, function(f) {
        f$boot_info$se_index_sum_by_group_neg[[grp]] %||% NA_real_
      }, numeric(1))

      pos_valid <- is.finite(pos_est) & is.finite(pos_se) & pos_se >= 0
      neg_valid <- is.finite(neg_est) & is.finite(neg_se) & neg_se >= 0

      group_results[[grp]] <- list(
        positive = if (sum(pos_valid) >= 1) {
          rubin_pool(
            pos_est[pos_valid],
            pos_se[pos_valid]^2,
            sum(pos_valid),
            df_complete = .mi_complete_df(fits, pos_valid)
          )
        } else {
          list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
               p_value = NA, fmi = NA, df = NA, t_stat = NA)
        },
        negative = if (sum(neg_valid) >= 1) {
          rubin_pool(
            neg_est[neg_valid],
            neg_se[neg_valid]^2,
            sum(neg_valid),
            df_complete = .mi_complete_df(fits, neg_valid)
          )
        } else {
          list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
               p_value = NA, fmi = NA, df = NA, t_stat = NA)
        }
      )
    }
    result$group_results <- group_results
    result$wqs_pos <- NULL
    result$wqs_neg <- NULL
  } else {
    pos_est <- vapply(fits, function(f) f$boot_info$mean_index_sum_pos %||% NA_real_, numeric(1))
    pos_se <- vapply(fits, function(f) f$boot_info$se_index_sum_pos %||% NA_real_, numeric(1))
    neg_est <- vapply(fits, function(f) f$boot_info$mean_index_sum_neg %||% NA_real_, numeric(1))
    neg_se <- vapply(fits, function(f) f$boot_info$se_index_sum_neg %||% NA_real_, numeric(1))

    pos_valid <- is.finite(pos_est) & is.finite(pos_se) & pos_se >= 0
    neg_valid <- is.finite(neg_est) & is.finite(neg_se) & neg_se >= 0

    result$wqs_pos <- if (sum(pos_valid) >= 1) {
      rubin_pool(pos_est[pos_valid], pos_se[pos_valid]^2, sum(pos_valid),
                 df_complete = .mi_complete_df(fits, pos_valid))
    } else {
      list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
           p_value = NA, fmi = NA, df = NA, t_stat = NA)
    }
    result$wqs_neg <- if (sum(neg_valid) >= 1) {
      rubin_pool(neg_est[neg_valid], neg_se[neg_valid]^2, sum(neg_valid),
                 df_complete = .mi_complete_df(fits, neg_valid))
    } else {
      list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
           p_value = NA, fmi = NA, df = NA, t_stat = NA)
    }
    result$group_results <- NULL
  }

  cov_names <- first_fit$cov_names %||% character(0)
  if (length(cov_names) > 0) {
    result$covariates <- lapply(cov_names, function(cov_nm) {
      est <- vapply(fits, function(f) f$boot_info$mean_cov_coef[[cov_nm]] %||% NA_real_, numeric(1))
      se <- vapply(fits, function(f) f$boot_info$se_cov_coef[[cov_nm]] %||% NA_real_, numeric(1))
      valid <- is.finite(est) & is.finite(se) & se >= 0
      if (sum(valid) >= 1) {
        rubin_pool(est[valid], se[valid]^2, sum(valid),
                   df_complete = .mi_complete_df(fits, valid))
      } else {
        list(estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
             p_value = NA, fmi = NA, df = NA, t_stat = NA)
      }
    })
    names(result$covariates) <- cov_names
  } else {
    result$covariates <- NULL
  }

  if (verbose) {
    message("  Pooled bootstrap-derived summaries across ", length(fits), " imputations")
  }

  result
}


#' Rubin's Rules for Pooling
#'
#' @param estimates Vector of estimates from each imputation
#' @param variances Vector of variances from each imputation
#' @param m Number of imputations
#' @param conf_level Confidence level for intervals (default: 0.95)
#' @param df_complete Optional complete-data residual degrees of freedom used
#'   for Barnard-Rubin adjustment.
#' @return A list with pooled estimate, SE, CI, and p-value
#' @keywords internal
rubin_pool <- function(estimates, variances, m, conf_level = 0.95,
                       df_complete = NA_real_) {
  
  # If m=1, return simply
  if (m == 1) {
    se <- sqrt(variances[1])
    Q_bar <- estimates[1]
    alpha <- 1 - conf_level
    z_crit <- qnorm(1 - alpha / 2)
    
    return(list(
      estimate = Q_bar,
      se = se,
      ci_lower = Q_bar - z_crit * se,
      ci_upper = Q_bar + z_crit * se,
      conf_level = conf_level,
      df = Inf,  # Normal distribution approximation
      t_stat = Q_bar / se,
      p_value = 2 * pnorm(-abs(Q_bar / se)),
      fmi = NA,  # FMI cannot be calculated for m=1
      m = 1
    ))
  }
  
  # Pooled estimate
  Q_bar <- mean(estimates)
  
  # Within-imputation variance (average)
  U_bar <- mean(variances)
  
  # Between-imputation variance
  B <- var(estimates)
  
  # Total variance
  T <- U_bar + (1 + 1/m) * B
  
  # Standard errors
  se <- sqrt(T)
  
  # Degrees of freedom (Barnard-Rubin)
  r <- (1 + 1/m) * B / U_bar
  
  if (is.finite(r) && r > 0 && U_bar > 0) {
    # Adjusted degrees of freedom
    df_old <- (m - 1) * (1 + 1/r)^2
    
    if (is.finite(df_complete) && df_complete > 0 && is.finite(T) && T > 0) {
      lambda <- (1 + 1/m) * B / T
      df_obs <- (df_complete + 1) / (df_complete + 3) * df_complete * (1 - lambda)
      df <- (df_old * df_obs) / (df_old + df_obs)
    } else {
      df <- df_old
    }
    df <- max(df, 1)  # Minimum 1
  } else {
    df <- m - 1
    df <- max(df, 1)
  }
  
  # Confidence interval (t distribution)
  alpha <- 1 - conf_level
  t_crit <- qt(1 - alpha / 2, df)
  ci_lower <- Q_bar - t_crit * se
  ci_upper <- Q_bar + t_crit * se
  
  # p-value (two-sided test)
  t_stat <- Q_bar / se
  p_value <- 2 * pt(-abs(t_stat), df)
  
  # Fraction of Missing Information (FMI)
  if (is.finite(r) && r >= 0) {
    fmi <- (r + 2 / (df + 3)) / (r + 1)
  } else {
    fmi <- NA
  }
  
  return(list(
    estimate = Q_bar,
    se = se,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    conf_level = conf_level,
    df = df,
    t_stat = t_stat,
    p_value = p_value,
    fmi = fmi,
    m = m
  ))
}


#' @keywords internal
.complete_data_df_from_sglwqs_fit <- function(fit) {
  info <- fit$validation_info
  if (is.null(info)) {
    info <- fit$refit_info
  }
  if (is.null(info)) {
    return(NA_real_)
  }
  
  fit_obj <- if (!is.null(info$refit_fit)) info$refit_fit else info$glm_fit
  if (is.null(fit_obj)) {
    return(NA_real_)
  }
  
  df <- tryCatch(stats::df.residual(fit_obj), error = function(e) NA_real_)
  if (!is.finite(df) || df <= 0) {
    return(NA_real_)
  }
  
  as.numeric(df)
}


#' @keywords internal
.mi_complete_df <- function(fits, valid = NULL) {
  fit_subset <- if (is.null(valid)) fits else fits[valid]
  if (length(fit_subset) == 0) {
    return(NA_real_)
  }
  
  dfs <- vapply(fit_subset, .complete_data_df_from_sglwqs_fit, numeric(1))
  dfs <- dfs[is.finite(dfs) & dfs > 0]
  
  if (length(dfs) == 0) {
    return(NA_real_)
  }
  
  min(dfs)
}


#' Model-Aware CI Multiplier for sglwqs_mids
#'
#' Extracts a refit/validation GLM from the first successful fit in the
#' sglwqs_mids object and delegates to \code{.prediction_interval_multiplier()}.
#' Falls back to \code{qnorm} when no model object is available.
#'
#' @param object A \code{sglwqs_mids} object.
#' @param conf_level Confidence level in (0, 1).
#' @return A single numeric critical value.
#' @keywords internal
.mi_ci_multiplier <- function(object, conf_level) {
  # Try to find a refit/validation GLM from any successful fit
  fits <- Filter(Negate(is.null), object$fits)
  for (f in fits) {
    refit_fit <- NULL
    if (!is.null(f$refit_info$refit_fit)) {
      refit_fit <- f$refit_info$refit_fit
    } else if (!is.null(f$validation_info$refit_fit)) {
      refit_fit <- f$validation_info$refit_fit
    } else if (!is.null(f$validation_info$glm_fit)) {
      refit_fit <- f$validation_info$glm_fit
    }
    if (!is.null(refit_fit)) {
      return(.prediction_interval_multiplier(refit_fit, conf_level))
    }
  }
  # Fallback: no model object available
  stats::qnorm(1 - (1 - conf_level) / 2)
}


#' Print sglwqs_mids Object
#'
#' @param x An object of class \code{"sglwqs_mids"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.sglwqs_mids <- function(x, ...) {
  cat("SGL-WQS with Multiple Imputation\n")
  cat(strrep("=", 50), "\n\n")
  
  cat("Imputations: ", x$n_successful, "/", x$m, " successful\n", sep = "")
  cat("Exposures:", length(x$exposure_vars), "\n")
  cat("Family:", x$family, "\n")
  cat("Bootstrap:", ifelse(x$bootstrap, "Yes", "No"), "\n")
  cat("Validation:", ifelse(x$validation, "Yes", "No"), "\n")
  
  eff_groups <- .get_effective_groups(x)
  if (!is.null(eff_groups)) {
    cat("Groups:", length(eff_groups), "\n")
  }

  cat("\n--- Pooled Index Sums ---\n")
  cat("Positive:", round(x$pooled$pos_index_sum, 4), "\n")
  cat("Negative:", round(x$pooled$neg_index_sum, 4), "\n")
  
  # Number of selected variables
  n_pos_selected <- sum(x$pooled$pos_weights > 0)
  n_neg_selected <- sum(x$pooled$neg_weights > 0)
  cat("\n--- Selected Variables (pooled) ---\n")
  cat("Positive direction:", n_pos_selected, "/", length(x$pooled$pos_weights), "\n")
  cat("Negative direction:", n_neg_selected, "/", length(x$pooled$neg_weights), "\n")
  
  inf <- .get_mi_pooled_inference(x)
  if (!is.null(inf)) {
    header <- if (identical(inf$source_used, "refit_info")) {
      "--- Pooled Refit Inference (Rubin's Rules) ---"
    } else {
      "--- Pooled Validation Inference (Rubin's Rules) ---"
    }
    cat("\n", header, "\n", sep = "")
    val <- inf
    
    if (isTRUE(val$has_group_results) && !is.null(val$group_results)) {
      # Group-specific results
      for (grp in names(val$group_results)) {
        cat("\n", grp, ":\n", sep = "")
        pos <- val$group_results[[grp]]$positive
        neg <- val$group_results[[grp]]$negative
        
        if (!is.na(pos$estimate)) {
          cat(sprintf("  Positive: estimate=%.4f, SE=%.4f, p=%s\n",
                      pos$estimate, pos$se, format.pval(pos$p_value, digits = 3)))
        } else {
          cat("  Positive: not estimated (coefficients = 0)\n")
        }
        
        if (!is.na(neg$estimate)) {
          cat(sprintf("  Negative: estimate=%.4f, SE=%.4f, p=%s\n",
                      neg$estimate, neg$se, format.pval(neg$p_value, digits = 3)))
        } else {
          cat("  Negative: not estimated (coefficients = 0)\n")
        }
      }
    } else {
      # Overall results
      pos <- val$wqs_pos
      neg <- val$wqs_neg
      
      if (!is.na(pos$estimate)) {
        cat(sprintf("WQS Positive: estimate=%.4f, SE=%.4f, p=%s\n",
                    pos$estimate, pos$se, format.pval(pos$p_value, digits = 3)))
      } else {
        cat("WQS Positive: not estimated\n")
      }
      if (!is.na(neg$estimate)) {
        cat(sprintf("WQS Negative: estimate=%.4f, SE=%.4f, p=%s\n",
                    neg$estimate, neg$se, format.pval(neg$p_value, digits = 3)))
      } else {
        cat("WQS Negative: not estimated\n")
      }
    }
  } else if (!is.null(x$pooled$bootstrap_inference)) {
    cat("\n--- Pooled Bootstrap Inference (Rubin's Rules) ---\n")
    boot_inf <- x$pooled$bootstrap_inference
    if (isTRUE(boot_inf$has_group_results) && !is.null(boot_inf$group_results)) {
      for (grp in names(boot_inf$group_results)) {
        cat("\n", grp, ":\n", sep = "")
        pos <- boot_inf$group_results[[grp]]$positive
        neg <- boot_inf$group_results[[grp]]$negative
        cat(sprintf("  Positive: estimate=%.4f, SE=%.4f\n", pos$estimate, pos$se))
        cat(sprintf("  Negative: estimate=%.4f, SE=%.4f\n", neg$estimate, neg$se))
      }
    } else {
      if (!is.null(boot_inf$wqs_pos)) {
        cat(sprintf("WQS Positive: estimate=%.4f, SE=%.4f\n", boot_inf$wqs_pos$estimate, boot_inf$wqs_pos$se))
      }
      if (!is.null(boot_inf$wqs_neg)) {
        cat(sprintf("WQS Negative: estimate=%.4f, SE=%.4f\n", boot_inf$wqs_neg$estimate, boot_inf$wqs_neg$se))
      }
    }
  }

  if (!is.null(x$diagnostics)) {
    cat("\n")
    print(x$diagnostics)
    cat("See `compute_diagnostics(fit)` for access to these values programmatically.\n")
  }
  
  invisible(x)
}


#' Summary of sglwqs_mids Object
#'
#' @param object An object of class \code{"sglwqs_mids"}.
#' @param ... Additional arguments passed to \code{summary()}.
#' @export
summary.sglwqs_mids <- function(object, ...) {
  cat("SGL-WQS with Multiple Imputation - Summary\n")
  cat(strrep("=", 60), "\n\n")

  cat("Model Information:\n")
  cat("  Imputations:", object$n_successful, "/", object$m, "\n")
  cat("  Exposures:", length(object$exposure_vars), "\n")
  cat("  Outcome:", object$outcome_var, "\n")
  if (!is.null(object$covariate_vars)) {
    cat("  Covariates:", paste(object$covariate_vars, collapse = ", "), "\n")
  }
  cat("  Family:", object$family, "\n")
  cat("  Bootstrap:", ifelse(object$bootstrap, "Yes", "No"), "\n")
  cat("  Validation:", ifelse(object$validation, "Yes", "No"), "\n")
  if (!is.null(object$penalize_covariates)) {
    cat("  Penalize covariates:", ifelse(isTRUE(object$penalize_covariates), "Yes", "No"), "\n")
  }

  eff_groups <- .get_effective_groups(object)
  if (!is.null(eff_groups)) {
    cat("\nChemical Groups:\n")
    for (grp_name in names(eff_groups)) {
      cat(sprintf("  %s: %d variables\n", grp_name, length(eff_groups[[grp_name]])))
    }
  }

  if (!is.null(object$diagnostics)) {
    cat("\n")
    print(object$diagnostics)
    cat("See `compute_diagnostics(fit)` for access to these values programmatically.\n")
  }

  coef_df <- if (!is.null(.get_mi_pooled_inference(object))) {
    .build_mi_validation_coef_table(object, conf_level = 0.95)
  } else if (!is.null(object$pooled$bootstrap_inference)) {
    boot_df <- summary_inference(object, conf_level = 0.95)
    data.frame(
      Term = boot_df$term,
      Estimate = boot_df$estimate,
      `Std.Error` = boot_df$std_error,
      `t value` = ifelse(is.na(boot_df$std_error) | boot_df$std_error == 0,
                         NA_real_, boot_df$estimate / boot_df$std_error),
      DF = NA_real_,
      `Pr(>|t|)` = boot_df$p_value,
      FMI = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    .build_mi_training_coef_table(object)
  }

  if (nrow(coef_df) > 0) {
    coef_df$Signif <- ""
    coef_df$Signif[!is.na(coef_df$`Pr(>|t|)`) & coef_df$`Pr(>|t|)` < 0.001] <- "***"
    coef_df$Signif[!is.na(coef_df$`Pr(>|t|)`) & coef_df$`Pr(>|t|)` >= 0.001 & coef_df$`Pr(>|t|)` < 0.01] <- "**"
    coef_df$Signif[!is.na(coef_df$`Pr(>|t|)`) & coef_df$`Pr(>|t|)` >= 0.01 & coef_df$`Pr(>|t|)` < 0.05] <- "*"
    coef_df$Signif[!is.na(coef_df$`Pr(>|t|)`) & coef_df$`Pr(>|t|)` >= 0.05 & coef_df$`Pr(>|t|)` < 0.1] <- "."

    print_df <- coef_df
    print_df$Estimate <- ifelse(is.na(coef_df$Estimate), "---", sprintf("%.6f", coef_df$Estimate))
    print_df$`Std.Error` <- ifelse(is.na(coef_df$`Std.Error`), "---", sprintf("%.6f", coef_df$`Std.Error`))
    print_df$`t value` <- ifelse(is.na(coef_df$`t value`), "---", sprintf("%.3f", coef_df$`t value`))
    print_df$DF <- ifelse(is.na(coef_df$DF), "---",
                          ifelse(is.infinite(coef_df$DF), "Inf", sprintf("%.1f", coef_df$DF)))
    print_df$`Pr(>|t|)` <- ifelse(is.na(coef_df$`Pr(>|t|)`), "---", format.pval(coef_df$`Pr(>|t|)`, digits = 3))
    print_df$FMI <- ifelse(is.na(coef_df$FMI), "---", sprintf("%.1f%%", coef_df$FMI * 100))

    cat("\n", strrep("=", 40), "\n", sep = "")
    if (!is.null(.get_mi_pooled_inference(object))) {
      inf <- .get_mi_pooled_inference(object)
      label <- if (identical(inf$source_used, "refit_info")) {
        "Rubin-Pooled Refit Coefficients"
      } else {
        "Rubin-Pooled Validation Coefficients"
      }
      cat(label, "\n")
    } else if (!is.null(object$pooled$bootstrap_inference)) {
      cat("Rubin-Pooled Bootstrap Inference\n")
    } else {
      cat("Pooled Training Coefficients\n")
    }
    cat(strrep("=", 40), "\n\n")
    print(print_df, row.names = FALSE, right = FALSE)

    if (!is.null(.get_mi_pooled_inference(object))) {
      inf <- .get_mi_pooled_inference(object)
      cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
      cat(
        "Note: Rubin's rules are applied term-wise to the ",
        if (identical(inf$source_used, "refit_info")) "refit" else "validation",
        " GLM coefficients, including covariates.\n",
        sep = ""
      )
      if (any(is.na(coef_df$Estimate) | is.na(coef_df$`Std.Error`))) {
        cat("      '---' indicates terms not estimable across imputations (often regularized to zero).\n")
      }
    } else {
      cat("\nNote: Without validation, coefficients are pooled point summaries from the training fit.\n")
      cat("      Standard errors for training coefficients are not Rubin-pooled in this summary.\n")
    }
  }

  cat("\n", strrep("=", 40), "\n", sep = "")
  cat("Pooled Weights by Group\n")
  cat(strrep("=", 40), "\n")

  pooled <- object$pooled

  if (!is.null(eff_groups)) {
    for (grp_name in names(eff_groups)) {
      grp_vars <- eff_groups[[grp_name]]
      pos_index_val <- pooled$pos_index_sum_by_group[[grp_name]]
      neg_index_val <- pooled$neg_index_sum_by_group[[grp_name]]
      if (!is.numeric(pos_index_val) || length(pos_index_val) != 1 || !is.finite(pos_index_val)) {
        pos_index_val <- object$pooled$bootstrap_inference$group_results[[grp_name]]$positive$estimate %||% NA_real_
      }
      if (!is.numeric(neg_index_val) || length(neg_index_val) != 1 || !is.finite(neg_index_val)) {
        neg_index_val <- object$pooled$bootstrap_inference$group_results[[grp_name]]$negative$estimate %||% NA_real_
      }

      cat("\n--- ", grp_name, " ---\n", sep = "")
      cat("Positive Index Sum:", round(pos_index_val, 4), "\n")

      pos_w <- pooled$pos_weights[grp_vars]
      pos_nonzero <- pos_w[pos_w > 0]
      if (length(pos_nonzero) > 0) {
        pos_sorted <- sort(pos_nonzero, decreasing = TRUE)
        for (nm in names(pos_sorted)) {
          line <- sprintf("  %-18s: %.4f", nm, pos_sorted[nm])
          mi_freq <- pooled$mi_selection_freq_pos[nm]
          line <- paste0(line, sprintf(" (MI sel: %.0f%%)", mi_freq * 100))
          cat(line, "\n")
        }
      }

      cat("Negative Index Sum:", round(neg_index_val, 4), "\n")

      neg_w <- pooled$neg_weights[grp_vars]
      neg_nonzero <- neg_w[neg_w > 0]
      if (length(neg_nonzero) > 0) {
        neg_sorted <- sort(neg_nonzero, decreasing = TRUE)
        for (nm in names(neg_sorted)) {
          line <- sprintf("  %-18s: %.4f", nm, neg_sorted[nm])
          mi_freq <- pooled$mi_selection_freq_neg[nm]
          line <- paste0(line, sprintf(" (MI sel: %.0f%%)", mi_freq * 100))
          cat(line, "\n")
        }
      }
    }
  } else {
    cat("\nPositive Direction (Index Sum =", round(pooled$pos_index_sum, 4), "):\n")
    pos_nonzero <- pooled$pos_weights[pooled$pos_weights > 0]
    if (length(pos_nonzero) > 0) {
      pos_sorted <- sort(pos_nonzero, decreasing = TRUE)
      for (nm in names(pos_sorted)) {
        cat(sprintf("  %-18s: %.4f (MI sel: %.0f%%)\n",
                    nm, pos_sorted[nm], pooled$mi_selection_freq_pos[nm] * 100))
      }
    }

    cat("\nNegative Direction (Index Sum =", round(pooled$neg_index_sum, 4), "):\n")
    neg_nonzero <- pooled$neg_weights[pooled$neg_weights > 0]
    if (length(neg_nonzero) > 0) {
      neg_sorted <- sort(neg_nonzero, decreasing = TRUE)
      for (nm in names(neg_sorted)) {
        cat(sprintf("  %-18s: %.4f (MI sel: %.0f%%)\n",
                    nm, neg_sorted[nm], pooled$mi_selection_freq_neg[nm] * 100))
      }
    }
  }

  invisible(object)
}


#' Print Method for sglwqs_mids_inference
#'
#' @param x An object of class \code{"sglwqs_mids_inference"}.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#' @export
print.sglwqs_mids_inference <- function(x, ...) {
  cat("SGL-WQS Pooled Inference Table (Multiple Imputation)\n")
  cat(strrep("=", 60), "\n\n")
  
  m <- attr(x, "m")
  n_successful <- attr(x, "n_successful")
  cat("Imputations: ", n_successful, "/", m, "\n", sep = "")
  
  # Validation results
  val_type <- attr(x, "validation_type")
  if (!is.null(val_type)) {
    cat("\n--- Pooled Mixture Effect Estimates (Rubin's Rules) ---\n")
    
    if (val_type == "group") {
      group_results <- attr(x, "group_results")
      
      # Display in table format
      cat("\n")
      cat(sprintf("%-12s %-10s %10s %10s %12s %8s\n", 
                  "Group", "Direction", "Estimate", "SE", "P-value", "FMI"))
      cat(strrep("-", 70), "\n")
      
      for (grp in names(group_results)) {
        pos <- group_results[[grp]]$positive
        neg <- group_results[[grp]]$negative
        
        if (!is.na(pos$estimate)) {
          cat(sprintf("%-12s %-10s %10.4f %10.4f %12s %7.1f%%\n",
                      grp, "Positive", pos$estimate, pos$se, 
                      format.pval(pos$p_value, digits = 3),
                      ifelse(is.na(pos$fmi), NA, pos$fmi * 100)))
        } else {
          cat(sprintf("%-12s %-10s %10s %10s %12s %8s\n",
                      grp, "Positive", "---", "---", "---", "---"))
        }
        
        if (!is.na(neg$estimate)) {
          cat(sprintf("%-12s %-10s %10.4f %10.4f %12s %7.1f%%\n",
                      grp, "Negative", neg$estimate, neg$se,
                      format.pval(neg$p_value, digits = 3),
                      ifelse(is.na(neg$fmi), NA, neg$fmi * 100)))
        } else {
          cat(sprintf("%-12s %-10s %10s %10s %12s %8s\n",
                      grp, "Negative", "---", "---", "---", "---"))
        }
      }
      cat("\nNote: '---' indicates coefficients = 0 due to regularization\n")
      
    } else {
      pos <- attr(x, "wqs_pos")
      neg <- attr(x, "wqs_neg")
      
      cat("\n")
      cat(sprintf("%-12s %10s %10s %12s %8s\n", 
                  "Direction", "Estimate", "SE", "P-value", "FMI"))
      cat(strrep("-", 55), "\n")
      
      if (!is.na(pos$estimate)) {
        cat(sprintf("%-12s %10.4f %10.4f %12s %7.1f%%\n",
                    "Positive", pos$estimate, pos$se,
                    format.pval(pos$p_value, digits = 3),
                    ifelse(is.na(pos$fmi), NA, pos$fmi * 100)))
      }
      if (!is.na(neg$estimate)) {
        cat(sprintf("%-12s %10.4f %10.4f %12s %7.1f%%\n",
                    "Negative", neg$estimate, neg$se,
                    format.pval(neg$p_value, digits = 3),
                    ifelse(is.na(neg$fmi), NA, neg$fmi * 100)))
      }
    }
    
    cat("\nTraining n:", round(attr(x, "n_train")), 
        "| Validation n:", round(attr(x, "n_val")), "\n")
  }
  
  cat("\n--- Variable Weights ---\n\n")
  # Display dataframe part
  print.data.frame(x, row.names = FALSE, ...)
  
  invisible(x)
}



#' Compute Unified Variable Order for sglwqs_mids Plots
#' @keywords internal
.compute_mi_unified_var_order <- function(object,
                                          top_n = 10,
                                          select_by = c("combined", "weight", "mi_selection_freq"),
                                          order_by = select_by,
                                          freq_col = "mi_selection_freq",
                                          require_nonzero_weight = TRUE,
                                          respect_groups = TRUE) {

  metric_choices <- c("combined", "weight", "mi_selection_freq")
  select_by <- match.arg(select_by, choices = metric_choices)
  order_by <- match.arg(order_by, choices = metric_choices)

  empty_result <- list(
    var_order = character(0),
    df = data.frame(stringsAsFactors = FALSE)
  )

  weights_df <- extract_weights(object, direction = "both", by_group = TRUE)

  if (is.null(weights_df) || nrow(weights_df) == 0) {
    return(empty_result)
  }

  if (!freq_col %in% names(weights_df)) {
    fallback_col <- if ("mi_selection_freq" %in% names(weights_df)) "mi_selection_freq" else NULL
    if (is.null(fallback_col)) {
      weights_df[[freq_col]] <- 1
    } else {
      weights_df[[freq_col]] <- weights_df[[fallback_col]]
    }
  }

  if (require_nonzero_weight) {
    weights_df <- weights_df[weights_df$weight > 0, , drop = FALSE]
  }
  if (nrow(weights_df) == 0) {
    return(empty_result)
  }

  pos_df <- weights_df[
    weights_df$direction == "positive",
    c("variable", "weight", freq_col),
    drop = FALSE
  ]
  names(pos_df) <- c("variable", "pos_weight", "pos_freq")

  neg_df <- weights_df[
    weights_df$direction == "negative",
    c("variable", "weight", freq_col),
    drop = FALSE
  ]
  names(neg_df) <- c("variable", "neg_weight", "neg_freq")

  var_df <- merge(pos_df, neg_df, by = "variable", all = TRUE, sort = FALSE)

  for (nm in c("pos_weight", "neg_weight", "pos_freq", "neg_freq")) {
    if (!nm %in% names(var_df)) {
      var_df[[nm]] <- 0
    }
    var_df[[nm]][is.na(var_df[[nm]])] <- 0
  }

  groups_to_use <- NULL
  if (respect_groups) {
    groups_to_use <- .get_effective_groups(object)
  }

  if (!is.null(groups_to_use) && "group" %in% names(weights_df)) {
    group_map <- unique(weights_df[, c("variable", "group"), drop = FALSE])
    var_df <- merge(var_df, group_map, by = "variable", all.x = TRUE, sort = FALSE)
    var_df$group[is.na(var_df$group)] <- "Other"
    group_order <- .group_display_order(groups_to_use, unique(var_df$group))
  } else {
    var_df$group <- "All"
    group_order <- "All"
  }

  var_df$weight_max <- pmax(var_df$pos_weight, var_df$neg_weight, na.rm = TRUE)
  var_df$freq_max <- pmax(var_df$pos_freq, var_df$neg_freq, na.rm = TRUE)

  score_value <- function(type, df) {
    switch(type,
      combined = df$weight_max * df$freq_max,
      weight = df$weight_max,
      mi_selection_freq = df$freq_max
    )
  }

  var_df$select_val <- score_value(select_by, var_df)
  var_df$order_val <- score_value(order_by, var_df)

  select_and_order_group <- function(g) {
    if (nrow(g) == 0) {
      return(g)
    }

    g <- g[order(-g$select_val, -g$order_val, -g$weight_max, -g$freq_max, g$variable), , drop = FALSE]

    if (!is.null(top_n)) {
      g <- head(g, top_n)
    }

    g <- g[order(-g$order_val, -g$select_val, -g$weight_max, -g$freq_max, g$variable), , drop = FALSE]
    g
  }

  df_list <- lapply(group_order, function(grp) {
    select_and_order_group(var_df[var_df$group == grp, , drop = FALSE])
  })
  df_list <- Filter(function(x) nrow(x) > 0, df_list)

  if (length(df_list) == 0) {
    return(empty_result)
  }

  df_top <- do.call(rbind, df_list)
  rownames(df_top) <- NULL

  list(
    var_order = df_top$variable,
    df = df_top
  )
}


#' @rdname plot_selection_frequency
#' @param source For \code{sglwqs_mids} objects, selection-frequency source:
#'   \code{"mi"} (default) for across-imputation selection frequency or
#'   \code{"bootstrap"} for pooled bootstrap selection frequency.
#' @param pos_color Color for positive direction.
#' @param neg_color Color for negative direction.
#' @param base_size Base font size.
#' @param title Optional plot title. When \code{NULL}, a default title is used.
#' @export
plot_selection_frequency.sglwqs_mids <- function(object,
                                                 top_n = 15,
                                                 min_freq = 0.1,
                                                 direction = c("both", "positive", "negative"),
                                                 sort_by = c("frequency", "combined", "weight"),
                                                 facet_by_group = NULL,
                                                 show_threshold = TRUE,
                                                 source = c("mi", "bootstrap"),
                                                 pos_color = "#0072B2",
                                                 neg_color = "#E69F00",
                                                 base_size = 11,
                                                 title = NULL,
                                                 ...) {

  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }

  direction <- match.arg(direction)
  sort_by <- match.arg(sort_by)
  source <- match.arg(source)

  has_groups <- !is.null(.get_effective_groups(object))
  if (is.null(facet_by_group)) {
    facet_by_group <- has_groups
  }

  freq_col <- if (source == "bootstrap") "boot_selection_freq" else "mi_selection_freq"
  weights_df <- extract_weights(object, direction = "both", by_group = TRUE)

  if (is.null(weights_df) || nrow(weights_df) == 0) {
    message("No pooled selection-frequency data available to plot")
    return(invisible(NULL))
  }

  if (source == "bootstrap" && !freq_col %in% names(weights_df)) {
    stop("Pooled bootstrap selection frequencies are not available. Fit with bootstrap = TRUE.")
  }
  if (!freq_col %in% names(weights_df)) {
    weights_df[[freq_col]] <- 1
  }

  weights_df$selection_freq <- weights_df[[freq_col]]
  plot_data <- weights_df[weights_df$selection_freq >= min_freq, , drop = FALSE]

  if (direction == "positive") {
    plot_data <- plot_data[plot_data$direction == "positive", , drop = FALSE]
  } else if (direction == "negative") {
    plot_data <- plot_data[plot_data$direction == "negative", , drop = FALSE]
  }

  if (nrow(plot_data) == 0) {
    message("No variables with selection frequency >= ", min_freq * 100, "%")
    return(invisible(NULL))
  }

  if (!facet_by_group || !has_groups) {
    plot_data$group <- "All"
  }

  score_df <- do.call(rbind, lapply(split(plot_data, plot_data$variable), function(df) {
    data.frame(
      variable = df$variable[[1]],
      group = df$group[[1]],
      weight_max = max(df$weight, na.rm = TRUE),
      freq_max = max(df$selection_freq, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  score_df$sort_val <- switch(
    sort_by,
    frequency = score_df$freq_max,
    combined = score_df$freq_max * score_df$weight_max,
    weight = score_df$weight_max
  )

  if (facet_by_group) {
    group_order <- unique(score_df$group)
    ordered_list <- lapply(group_order, function(grp) {
      g <- score_df[score_df$group == grp, , drop = FALSE]
      g <- g[order(-g$sort_val, -g$weight_max, -g$freq_max, g$variable), , drop = FALSE]
      if (!is.null(top_n)) {
        g <- head(g, top_n)
      }
      g
    })
    ordered_list <- Filter(function(x) nrow(x) > 0, ordered_list)
    if (length(ordered_list) == 0) {
      message("No pooled selection-frequency data available to plot")
      return(invisible(NULL))
    }
    score_df <- do.call(rbind, ordered_list)
  } else {
    score_df <- score_df[order(-score_df$sort_val, -score_df$weight_max, -score_df$freq_max, score_df$variable), , drop = FALSE]
    if (!is.null(top_n)) {
      score_df <- head(score_df, top_n)
    }
  }

  var_info <- list(
    var_order = score_df$variable,
    df = score_df
  )

  if (is.null(title)) {
    title <- if (source == "bootstrap") {
      paste0("Pooled Bootstrap Selection Frequency (m=", object$n_successful, ")")
    } else {
      paste0("MI Selection Frequency (m=", object$n_successful, ")")
    }
  }

  subtitle <- if (source == "bootstrap") {
    "Average bootstrap selection frequency across successful imputations"
  } else {
    "Proportion of successful imputations with non-zero directional coefficient"
  }

  p <- .plot_mi_selection_butterfly(
    object = object,
    top_n = top_n,
    sort_by = "mi_selection_freq",
    var_info = var_info,
    pos_color = pos_color,
    neg_color = neg_color,
    base_size = base_size,
    freq_col = freq_col,
    y_label = if (source == "bootstrap") "Pooled Bootstrap Selection Frequency" else "MI Selection Frequency",
    show_threshold = show_threshold,
    show_legend = TRUE,
    title = title,
    subtitle = subtitle
  )

  if (direction != "both") {
    data_dir <- if (direction == "positive") "positive" else "negative"
    p$data <- p$data[p$data$direction == data_dir, , drop = FALSE]
  }

  return(p)
}


#' Summarize Pooled Bootstrap Results for sglwqs_mids
#' @keywords internal
.summarize_pooled_bootstrap <- function(object, conf_level = 0.95) {

  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }
  if (!isTRUE(object$bootstrap) || is.null(object$pooled$bootstrap)) {
    stop("Pooled bootstrap information is not available. Fit with bootstrap = TRUE.")
  }

  pooled <- object$pooled
  boot <- pooled$bootstrap

  var_names <- pooled$var_names
  if (is.null(var_names) || length(var_names) == 0) {
    var_names <- names(pooled$pos_weights)
  }
  if (is.null(var_names) || length(var_names) == 0) {
    var_names <- names(pooled$neg_weights)
  }
  if (is.null(var_names) || length(var_names) == 0) {
    stop("No pooled variable names found.")
  }

  n_vars <- length(var_names)
  safe_numeric <- function(x) {
    out <- as.numeric(x)
    if (is.null(out) || length(out) == 0) {
      stop("Missing pooled bootstrap vector in .summarize_pooled_bootstrap().")
    }
    if (length(out) == 1L) return(rep(out, n_vars))
    if (length(out) != n_vars) {
      stop("Length mismatch in .summarize_pooled_bootstrap(): expected ",
           n_vars, ", got ", length(out), ".")
    }
    out
  }

  # Model-aware CI multiplier: use refit/validation fit if available
  z_val <- .mi_ci_multiplier(object, conf_level)

  pos_df <- data.frame(
    variable = var_names,
    direction = "positive",
    weight = safe_numeric(pooled$pos_weights),
    mean_coef = safe_numeric(pooled$pos_coef),
    se_coef = safe_numeric(boot$se_pos_coef),
    selection_freq = if (!is.null(boot$selection_freq_pos)) safe_numeric(boot$selection_freq_pos) else safe_numeric(pooled$mi_selection_freq_pos),
    stringsAsFactors = FALSE
  )
  pos_df$ci_lower <- pos_df$mean_coef - z_val * pos_df$se_coef
  pos_df$ci_upper <- pos_df$mean_coef + z_val * pos_df$se_coef

  neg_df <- data.frame(
    variable = var_names,
    direction = "negative",
    weight = safe_numeric(pooled$neg_weights),
    mean_coef = safe_numeric(pooled$neg_coef),
    se_coef = safe_numeric(boot$se_neg_coef),
    selection_freq = if (!is.null(boot$selection_freq_neg)) safe_numeric(boot$selection_freq_neg) else safe_numeric(pooled$mi_selection_freq_neg),
    stringsAsFactors = FALSE
  )
  neg_df$ci_lower <- neg_df$mean_coef - z_val * neg_df$se_coef
  neg_df$ci_upper <- neg_df$mean_coef + z_val * neg_df$se_coef

  out <- rbind(pos_df, neg_df)

  groups_to_use <- .get_effective_groups(object)
  if (!is.null(groups_to_use)) {
    out$group <- .assign_group_labels(out$variable, groups_to_use)
  }

  attr(out, "n_imputations") <- object$n_successful
  attr(out, "conf_level") <- conf_level
  class(out) <- c("sglwqs_mids_bootstrap_summary", "data.frame")
  out
}


#' @rdname plot_bootstrap_ci
#' @export
plot_bootstrap_ci.sglwqs_mids <- function(object, top_n = 15, conf_level = 0.95,
                                          direction = c("both", "positive", "negative"),
                                          pos_color = "#2166AC",
                                          neg_color = "#B2182B",
                                          base_size = 11,
                                          title = NULL,
                                          ...) {

  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }
  if (!isTRUE(object$bootstrap) || is.null(object$pooled$bootstrap)) {
    stop("Pooled bootstrap information is not available. Fit with bootstrap = TRUE.")
  }

  direction <- match.arg(direction)
  plot_data <- .summarize_pooled_bootstrap(object, conf_level = conf_level)

  if (direction == "positive") {
    plot_data <- plot_data[plot_data$direction == "positive", , drop = FALSE]
  } else if (direction == "negative") {
    plot_data <- plot_data[plot_data$direction == "negative", , drop = FALSE]
  }

  plot_data <- plot_data[abs(plot_data$mean_coef) > 1e-6, , drop = FALSE]

  if (nrow(plot_data) == 0) {
    message("No non-zero coefficients to plot")
    return(invisible(NULL))
  }

  plot_data <- head(plot_data[order(-abs(plot_data$mean_coef)), , drop = FALSE], top_n)
  plot_data$variable <- factor(plot_data$variable,
                               levels = rev(unique(plot_data$variable[order(plot_data$mean_coef)])))

  if (is.null(title)) {
    title <- paste0("Pooled Bootstrap Coefficient Estimates (", conf_level * 100, "% CI)")
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = mean_coef)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lower, ymax = ci_upper, color = direction),
                             size = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Normal-approximation CI using imputation-averaged bootstrap SEs (m=", object$n_successful, ")"),
      caption = "Intervals are based on pooled bootstrap summaries rather than a Rubin-pooled bootstrap distribution.",
      x = NULL,
      y = "Coefficient"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::scale_color_manual(values = c("positive" = pos_color, "negative" = neg_color),
                                name = "Direction") +
    ggplot2::theme(
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(hjust = 0, size = base_size * 0.8,
                                           color = "gray40", face = "italic")
    )

  if (direction == "both") {
    p <- p + ggplot2::facet_wrap(~ direction, scales = "free")
  }

  return(p)
}


#' @keywords internal
.get_mi_inference_fits <- function(object,
                                   source = c("auto", "refit_info", "validation_info")) {
  source <- match.arg(source)
  Filter(function(f) {
    inherits(f, "sglwqs") && {
      info <- .resolve_mi_inference_info(f, source = source)
      !is.null(info) && !is.null(info$coef_table)
    }
  }, object$fits)
}

#' @keywords internal
.get_mi_inference_covariate_terms <- function(fits,
                                              source = c("auto", "refit_info", "validation_info")) {
  source <- match.arg(source)
  if (length(fits) == 0) {
    return(character())
  }

  cov_terms <- unique(unlist(lapply(fits, function(f) f$cov_names)))
  cov_terms <- cov_terms[!is.na(cov_terms) & nzchar(cov_terms)]

  if (length(cov_terms) == 0) {
    return(character())
  }

  ordered_terms <- unique(unlist(lapply(fits, function(f) {
    info <- .resolve_mi_inference_info(f, source = source)
    if (is.null(info) || is.null(info$coef_table)) {
      return(character())
    }
    rn <- rownames(info$coef_table)
    rn[rn %in% cov_terms]
  })))

  ordered_terms[!is.na(ordered_terms) & nzchar(ordered_terms)]
}

#' @keywords internal
.pool_mi_inference_term <- function(fits, term,
                                    source = c("auto", "refit_info", "validation_info"),
                                    conf_level = 0.95) {
  source <- match.arg(source)
  est <- sapply(fits, function(f) {
    info <- .resolve_mi_inference_info(f, source = source)
    if (is.null(info) || is.null(info$coef_table)) return(NA_real_)
    ct <- info$coef_table
    if (term %in% rownames(ct)) ct[term, "Estimate"] else NA_real_
  })
  se <- sapply(fits, function(f) {
    info <- .resolve_mi_inference_info(f, source = source)
    if (is.null(info) || is.null(info$coef_table)) return(NA_real_)
    ct <- info$coef_table
    if (term %in% rownames(ct)) ct[term, "Std. Error"] else NA_real_
  })

  valid <- is.finite(est) & is.finite(se) & se > 0
  m_used <- sum(valid)

  if (m_used == 0) {
    return(list(
      term = term,
      estimate = NA_real_,
      se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      conf_level = conf_level,
      df = NA_real_,
      t_stat = NA_real_,
      p_value = NA_real_,
      fmi = NA_real_,
      m = 0,
      m_used = 0
    ))
  }

  pooled <- rubin_pool(
    est[valid],
    se[valid]^2,
    m_used,
    conf_level = conf_level,
    df_complete = .mi_complete_df(fits, valid)
  )
  pooled$term <- term
  pooled$m_used <- m_used
  pooled
}

#' @keywords internal
.pool_mi_inference_terms_df <- function(fits, terms,
                                        source = c("auto", "refit_info", "validation_info"),
                                        conf_level = 0.95) {
  source <- match.arg(source)
  terms <- unique(as.character(terms))
  terms <- terms[!is.na(terms) & nzchar(terms)]

  if (length(fits) == 0 || length(terms) == 0) {
    return(data.frame())
  }

  pooled_rows <- lapply(terms, function(term) {
    pooled <- .pool_mi_inference_term(fits, term, source = source, conf_level = conf_level)
    if (pooled$m_used == 0) {
      return(NULL)
    }

    data.frame(
      term = term,
      estimate = pooled$estimate,
      se = pooled$se,
      t_stat = pooled$t_stat,
      df = pooled$df,
      p_value = pooled$p_value,
      fmi = pooled$fmi,
      ci_lower = pooled$ci_lower,
      ci_upper = pooled$ci_upper,
      m_used = pooled$m_used,
      stringsAsFactors = FALSE
    )
  })

  pooled_rows <- Filter(Negate(is.null), pooled_rows)
  if (length(pooled_rows) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, pooled_rows)
  rownames(out) <- NULL
  out
}

#' @keywords internal
.get_mi_validation_fits <- function(object) {
  .get_mi_inference_fits(object, source = "validation_info")
}

#' @keywords internal
.get_mi_validation_covariate_terms <- function(fits) {
  .get_mi_inference_covariate_terms(fits, source = "validation_info")
}

#' @keywords internal
.pool_mi_validation_term <- function(fits, term, conf_level = 0.95) {
  .pool_mi_inference_term(fits, term, source = "validation_info", conf_level = conf_level)
}

#' @keywords internal
.pool_mi_validation_terms_df <- function(fits, terms, conf_level = 0.95) {
  .pool_mi_inference_terms_df(fits, terms, source = "validation_info", conf_level = conf_level)
}

#' @keywords internal
.pooled_validation_field <- function(x, name, default = NA_real_) {
  value <- x[[name]]
  if (is.null(value) || length(value) == 0) default else value
}

#' @keywords internal
.build_mi_training_coef_table <- function(object) {
  fits <- Filter(function(f) inherits(f, "sglwqs"), object$fits)
  if (length(fits) == 0) {
    return(data.frame())
  }

  rows <- list()
  intercepts <- sapply(fits, function(f) as.numeric(f$intercept)[1])
  intercept_est <- if (all(is.na(intercepts))) NA_real_ else mean(intercepts, na.rm = TRUE)

  rows[[length(rows) + 1]] <- data.frame(
    Term = "(Intercept)",
    Estimate = intercept_est,
    `Std.Error` = NA_real_,
    `t value` = NA_real_,
    DF = NA_real_,
    `Pr(>|t|)` = NA_real_,
    FMI = NA_real_,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  eff_groups <- .get_effective_groups(object)
  if (!is.null(eff_groups)) {
    for (grp in names(eff_groups)) {
      rows[[length(rows) + 1]] <- data.frame(
        Term = paste0(grp, " (positive)"),
        Estimate = object$pooled$pos_index_sum_by_group[[grp]],
        `Std.Error` = NA_real_,
        `t value` = NA_real_,
        DF = NA_real_,
        `Pr(>|t|)` = NA_real_,
        FMI = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        Term = paste0(grp, " (negative)"),
        Estimate = object$pooled$neg_index_sum_by_group[[grp]],
        `Std.Error` = NA_real_,
        `t value` = NA_real_,
        DF = NA_real_,
        `Pr(>|t|)` = NA_real_,
        FMI = NA_real_,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  } else {
    rows[[length(rows) + 1]] <- data.frame(
      Term = "WQS Positive (index sum)",
      Estimate = object$pooled$pos_index_sum,
      `Std.Error` = NA_real_,
      `t value` = NA_real_,
      DF = NA_real_,
      `Pr(>|t|)` = NA_real_,
      FMI = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      Term = "WQS Negative (index sum)",
      Estimate = object$pooled$neg_index_sum,
      `Std.Error` = NA_real_,
      `t value` = NA_real_,
      DF = NA_real_,
      `Pr(>|t|)` = NA_real_,
      FMI = NA_real_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (!is.null(object$pooled$covariates) && length(object$pooled$covariates) > 0) {
    cov_rows <- lapply(names(object$pooled$covariates), function(term) {
      p <- object$pooled$covariates[[term]]
      if (is.null(p) || is.na(p$estimate)) {
        return(NULL)
      }
      data.frame(
        Term = term,
        Estimate = p$estimate,
        `Std.Error` = p$se,
        `t value` = p$t_stat,
        DF = p$df,
        `Pr(>|t|)` = p$p_value,
        FMI = p$fmi,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    cov_rows <- Filter(Negate(is.null), cov_rows)
    if (length(cov_rows) > 0) {
      rows[[length(rows) + 1]] <- do.call(rbind, cov_rows)
    }
  } else {
    get_cov_for_pool <- function(f) {
      if (!is.null(f$boot_info$mean_cov_coef)) {
        return(f$boot_info$mean_cov_coef)
      }
      f$cov_coef
    }

    cov_terms <- unique(unlist(lapply(fits, function(f) names(get_cov_for_pool(f)))))
    cov_terms <- cov_terms[!is.na(cov_terms) & nzchar(cov_terms)]
    if (length(cov_terms) > 0) {
      cov_rows <- lapply(cov_terms, function(term) {
        est <- sapply(fits, function(f) {
          cv <- get_cov_for_pool(f)
          if (!is.null(cv) && term %in% names(cv)) cv[[term]] else NA_real_
        })
        if (all(is.na(est))) {
          return(NULL)
        }
        data.frame(
          Term = term,
          Estimate = mean(est, na.rm = TRUE),
          `Std.Error` = NA_real_,
          `t value` = NA_real_,
          DF = NA_real_,
          `Pr(>|t|)` = NA_real_,
          FMI = NA_real_,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      })
      cov_rows <- Filter(Negate(is.null), cov_rows)
      if (length(cov_rows) > 0) {
        rows[[length(rows) + 1]] <- do.call(rbind, cov_rows)
      }
    }
  }

  coef_df <- do.call(rbind, rows)
  rownames(coef_df) <- NULL
  coef_df
}

#' @keywords internal
.build_mi_validation_coef_table <- function(object, conf_level = 0.95) {
  inf <- .get_mi_pooled_inference(object)
  if (!inherits(object, "sglwqs_mids") || is.null(inf)) {
    return(data.frame())
  }

  source_used <- inf$source_used %||% "validation_info"
  fits <- .get_mi_inference_fits(object, source = source_used)
  if (length(fits) == 0) {
    return(data.frame())
  }

  rows <- list()

  if (isTRUE(inf$has_group_results) && !is.null(inf$group_results)) {
    for (grp in names(inf$group_results)) {
      pos <- inf$group_results[[grp]]$positive
      neg <- inf$group_results[[grp]]$negative

      rows[[length(rows) + 1]] <- data.frame(
        Term = paste0(grp, " (positive)"),
        Estimate = .pooled_validation_field(pos, "estimate"),
        `Std.Error` = .pooled_validation_field(pos, "se"),
        `t value` = .pooled_validation_field(pos, "t_stat"),
        DF = .pooled_validation_field(pos, "df"),
        `Pr(>|t|)` = .pooled_validation_field(pos, "p_value"),
        FMI = .pooled_validation_field(pos, "fmi"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      rows[[length(rows) + 1]] <- data.frame(
        Term = paste0(grp, " (negative)"),
        Estimate = .pooled_validation_field(neg, "estimate"),
        `Std.Error` = .pooled_validation_field(neg, "se"),
        `t value` = .pooled_validation_field(neg, "t_stat"),
        DF = .pooled_validation_field(neg, "df"),
        `Pr(>|t|)` = .pooled_validation_field(neg, "p_value"),
        FMI = .pooled_validation_field(neg, "fmi"),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  } else {
    rows[[length(rows) + 1]] <- data.frame(
      Term = "WQS Positive",
      Estimate = .pooled_validation_field(inf$wqs_pos, "estimate"),
      `Std.Error` = .pooled_validation_field(inf$wqs_pos, "se"),
      `t value` = .pooled_validation_field(inf$wqs_pos, "t_stat"),
      DF = .pooled_validation_field(inf$wqs_pos, "df"),
      `Pr(>|t|)` = .pooled_validation_field(inf$wqs_pos, "p_value"),
      FMI = .pooled_validation_field(inf$wqs_pos, "fmi"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      Term = "WQS Negative",
      Estimate = .pooled_validation_field(inf$wqs_neg, "estimate"),
      `Std.Error` = .pooled_validation_field(inf$wqs_neg, "se"),
      `t value` = .pooled_validation_field(inf$wqs_neg, "t_stat"),
      DF = .pooled_validation_field(inf$wqs_neg, "df"),
      `Pr(>|t|)` = .pooled_validation_field(inf$wqs_neg, "p_value"),
      FMI = .pooled_validation_field(inf$wqs_neg, "fmi"),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (!is.null(object$pooled$covariates) && length(object$pooled$covariates) > 0) {
    cov_rows <- lapply(names(object$pooled$covariates), function(term) {
      p <- object$pooled$covariates[[term]]
      if (is.null(p) || is.na(p$estimate)) {
        return(NULL)
      }
      data.frame(
        Term = term,
        Estimate = p$estimate,
        `Std.Error` = p$se,
        `t value` = p$t_stat,
        DF = p$df,
        `Pr(>|t|)` = p$p_value,
        FMI = p$fmi,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    cov_rows <- Filter(Negate(is.null), cov_rows)
    if (length(cov_rows) > 0) {
      rows[[length(rows) + 1]] <- do.call(rbind, cov_rows)
    }
  } else {
    cov_terms <- .get_mi_inference_covariate_terms(fits, source = source_used)
    cov_df <- .pool_mi_inference_terms_df(
      fits, cov_terms, source = source_used, conf_level = conf_level
    )
    if (nrow(cov_df) > 0) {
      cov_rows <- data.frame(
        Term = cov_df$term,
        Estimate = cov_df$estimate,
        `Std.Error` = cov_df$se,
        `t value` = cov_df$t_stat,
        DF = cov_df$df,
        `Pr(>|t|)` = cov_df$p_value,
        FMI = cov_df$fmi,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      rows[[length(rows) + 1]] <- cov_rows
    }
  }

  if (length(rows) == 0) {
    return(data.frame())
  }

  coef_df <- do.call(rbind, rows)
  coef_df <- coef_df[coef_df$Term != "(Intercept)", , drop = FALSE]
  rownames(coef_df) <- NULL
  coef_df
}

#' @keywords internal
.build_mi_validation_summary_table <- function(object, conf_level = 0.95) {
  coef_df <- .build_mi_validation_coef_table(object, conf_level = conf_level)
  if (nrow(coef_df) == 0) {
    stop("No pooled inference results available.")
  }

  group <- ifelse(grepl(" \\(positive\\)$", coef_df$Term), sub(" \\(positive\\)$", "", coef_df$Term),
                  ifelse(grepl(" \\(negative\\)$", coef_df$Term), sub(" \\(negative\\)$", "", coef_df$Term),
                         ifelse(grepl("^WQS ", coef_df$Term), "Overall",
                                ifelse(coef_df$Term %in% (object$covariate_vars %||% character(0)), "Covariate",
                                       coef_df$Term))))
  direction <- ifelse(grepl(" \\(positive\\)$", coef_df$Term) | coef_df$Term == "WQS Positive",
                      "positive",
                      ifelse(grepl(" \\(negative\\)$", coef_df$Term) | coef_df$Term == "WQS Negative",
                             "negative", NA_character_))
  type <- ifelse(is.na(direction), "Covariate", "WQS Index")

  out <- data.frame(
    term = coef_df$Term,
    group = group,
    direction = direction,
    estimate = coef_df$Estimate,
    std_error = coef_df$`Std.Error`,
    ci_lower = ifelse(is.na(coef_df$Estimate) | is.na(coef_df$`Std.Error`), NA_real_,
                      coef_df$Estimate - qt(1 - (1 - conf_level) / 2, pmax(coef_df$DF, 1), lower.tail = TRUE) * coef_df$`Std.Error`),
    ci_upper = ifelse(is.na(coef_df$Estimate) | is.na(coef_df$`Std.Error`), NA_real_,
                      coef_df$Estimate + qt(1 - (1 - conf_level) / 2, pmax(coef_df$DF, 1), lower.tail = TRUE) * coef_df$`Std.Error`),
    p_value = coef_df$`Pr(>|t|)`,
    type = type,
    stringsAsFactors = FALSE
  )
  out$signif <- ""
  out$signif[!is.na(out$p_value) & out$p_value < 0.1] <- "."
  out$signif[!is.na(out$p_value) & out$p_value < 0.05] <- "*"
  out$signif[!is.na(out$p_value) & out$p_value < 0.01] <- "**"
  out$signif[!is.na(out$p_value) & out$p_value < 0.001] <- "***"

  inf <- .get_mi_pooled_inference(object)
  attr(out, "n_train") <- inf$n_train %||% NA_real_
  attr(out, "n_val") <- inf$n_val %||% NA_real_
  attr(out, "conf_level") <- conf_level
  attr(out, "family") <- object$family
  attr(out, "group_inference") <- isTRUE(inf$has_group_results)
  attr(out, "formula") <- NULL
  attr(out, "source") <- inf$source_used %||% "validation_info"
  class(out) <- c("sglwqs_validation_summary", "data.frame")
  out
}

#' @keywords internal
.pool_mi_validation_covariates <- function(object, conf_level = 0.95) {
  inf <- .get_mi_pooled_inference(object)
  if (is.null(inf)) {
    return(data.frame())
  }

  if (!is.null(object$pooled$covariates) && length(object$pooled$covariates) > 0) {
    rows <- lapply(names(object$pooled$covariates), function(term) {
      p <- object$pooled$covariates[[term]]
      if (is.null(p) || is.na(p$estimate)) {
        return(NULL)
      }
      data.frame(
        term = term,
        group = "Covariate",
        direction = NA_character_,
        estimate = p$estimate,
        se = p$se,
        t_stat = p$t_stat,
        p_value = p$p_value,
        fmi = p$fmi,
        df = p$df,
        ci_lower = p$ci_lower,
        ci_upper = p$ci_upper,
        m_used = NA_real_,
        type = "Covariate",
        stringsAsFactors = FALSE
      )
    })
  } else {
    source_used <- inf$source_used %||% "validation_info"
    fits <- .get_mi_inference_fits(object, source = source_used)
    cov_terms <- .get_mi_inference_covariate_terms(fits, source = source_used)
    pooled_df <- .pool_mi_inference_terms_df(
      fits, cov_terms, source = source_used, conf_level = conf_level
    )
    if (nrow(pooled_df) == 0) {
      return(data.frame())
    }
    rows <- lapply(seq_len(nrow(pooled_df)), function(i) {
      data.frame(
        term = pooled_df$term[i],
        group = "Covariate",
        direction = NA_character_,
        estimate = pooled_df$estimate[i],
        se = pooled_df$se[i],
        t_stat = pooled_df$t_stat[i],
        p_value = pooled_df$p_value[i],
        fmi = pooled_df$fmi[i],
        df = pooled_df$df[i],
        ci_lower = pooled_df$ci_lower[i],
        ci_upper = pooled_df$ci_upper[i],
        m_used = pooled_df$m_used[i],
        type = "Covariate",
        stringsAsFactors = FALSE
      )
    })
  }
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) {
    return(data.frame())
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname plot_validation_results
#' @param pos_color Color for positive direction (default: "#0072B2").
#' @param neg_color Color for negative direction (default: "#E69F00").
#' @param show_significance For \code{sglwqs_mids} objects, significance annotation
#'   mode: \code{"none"}, \code{"asterisk"}, or \code{"pvalue"}.
#' @param base_size Base font size (default: 11).
#' @param title Optional plot title for the pooled validation forest plot.
#'   Use \code{NULL} to suppress the title.
#' @export
plot_validation_results.sglwqs_mids <- function(object,
                                                conf_level = 0.95,
                                                include_covariates = FALSE,
                                                pos_color = "#0072B2",
                                                neg_color = "#E69F00",
                                                show_significance = c("none", "asterisk", "pvalue"),
                                                base_size = 11,
                                                title = "Mixture Effect Estimates (Rubin's Rules)",
                                                ...) {

  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }

  if (is.null(.get_mi_pooled_inference(object))) {
    stop("Pooled validation information not available. Fit with validation = TRUE or refit = \"full\".")
  }

  show_significance <- match.arg(show_significance)

  p <- .plot_pooled_validation(
    object,
    pos_color = pos_color,
    neg_color = neg_color,
    base_size = base_size,
    conf_level = conf_level,
    include_covariates = include_covariates,
    show_significance = show_significance
  )

  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}


#' Combined Results Plot for sglwqs_mids (Publication-Ready)
#'
#' Creates a combined visualization consistent with plot_combined_results() style.
#' The MI selection-frequency panel and pooled validation forest plot can also be
#' generated separately with \code{plot_selection_frequency()} and
#' \code{plot_validation_results()}, respectively.
#'
#' @param object A sglwqs_mids object.
#' @param top_n Number of top variables to show per group (default: 10).
#' @param sort_by Character. Variable-selection criterion: "combined", "weight",
#'   or "mi_selection_freq". In the combined figure, panel order is aligned
#'   strictly to pooled weights so the MI frequency and weight panels share the
#'   same variable order.
#' @param pos_color Color for positive direction (default: "#0072B2").
#' @param neg_color Color for negative direction (default: "#E69F00").
#' @param layout Layout option for \code{plot_combined_results.sglwqs_mids}:
#'   \code{"auto"}, \code{"vertical"}, or \code{"horizontal"}.
#' @param show_significance Significance annotation mode for the pooled
#'   validation panel: \code{"none"}, \code{"asterisk"}, or
#'   \code{"pvalue"}.
#' @param conf_level Confidence level for the pooled validation panel.
#' @param include_covariates Logical. Include pooled covariate coefficients in
#'   the pooled validation panel when available.
#' @param base_size Base font size (default: 11).
#' @param ... Additional arguments passed to methods.
#'
#' @return A combined ggplot object (if patchwork is available) or list of plots.
#'
#' @export
plot_combined_results.sglwqs_mids <- function(object,
                                               top_n = 10,
                                               sort_by = c("combined", "weight", "mi_selection_freq"),
                                               pos_color = "#0072B2",
                                               neg_color = "#E69F00",
                                               layout = c("auto", "vertical", "horizontal"),
                                               show_significance = c("none", "asterisk", "pvalue"),
                                               conf_level = 0.95,
                                               include_covariates = FALSE,
                                               base_size = 11,
                                               ...) {
  
  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }
  
  sort_by <- match.arg(sort_by)
  layout <- match.arg(layout)
  show_significance <- match.arg(show_significance)
  
  # Select variables by the requested criterion, but keep the display order
  # strictly anchored to pooled weights so the frequency and weight panels align.
  var_info <- .compute_mi_unified_var_order(
    object,
    top_n = top_n,
    select_by = sort_by,
    order_by = "weight"
  )
  
  if (length(var_info$var_order) == 0) {
    message("No non-zero pooled weights to plot")
    return(invisible(NULL))
  }
  
  plots <- list()
  
  # A) MI selection-frequency plot
  plots$mi_selection <- .plot_mi_selection_butterfly(
    object,
    top_n = top_n,
    sort_by = sort_by,
    var_info = var_info,
    pos_color = pos_color,
    neg_color = neg_color,
    base_size = base_size
  ) + ggplot2::ggtitle("A) MI Selection Frequency")
  
  # B) Weight plot
  plots$weights <- plot_weights(
    object,
    top_n = top_n,
    sort_by = sort_by,
    pos_color = pos_color,
    neg_color = neg_color,
    show_freq = FALSE,
    base_size = base_size,
    var_info = var_info
  ) + ggplot2::ggtitle("B) Pooled Variable Weights")
  
  # C) Inference results (when available)
  if (!is.null(.get_mi_pooled_inference(object)) || !is.null(object$pooled$bootstrap_inference)) {
    plots$validation <- plot_inference_results(
      object,
      conf_level = conf_level,
      include_covariates = include_covariates,
      pos_color = pos_color,
      neg_color = neg_color,
      show_significance = show_significance,
      base_size = base_size,
      title = NULL
    ) + ggplot2::ggtitle("C) Mixture Effect Estimates (Rubin's Rules)")
  }
  
  if (layout == "auto") {
    if (length(plots) == 3) {
      layout <- "vertical"
    } else {
      layout <- "horizontal"
    }
  }
  
  # Combine with patchwork when available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    if (length(plots) >= 2) {
      if (layout == "vertical" && length(plots) == 3) {
        top_row <- patchwork::wrap_plots(plots$mi_selection, plots$weights, ncol = 2)
        combined <- top_row / plots$validation +
          patchwork::plot_layout(heights = c(2, 1))
      } else {
        combined <- patchwork::wrap_plots(plots, ncol = 2)
      }
      return(combined)
    }
  }
  
  message("Install 'patchwork' package for combined plot layout.")
  return(plots)
}

#' Plot MI Selection Frequency (Butterfly Chart)
#' @keywords internal
.plot_mi_selection_butterfly <- function(object, top_n = 10, sort_by = "combined",
                                         var_info = NULL,
                                         pos_color = "#0072B2", neg_color = "#E69F00",
                                         base_size = 11,
                                         freq_col = "mi_selection_freq",
                                         y_label = "MI Selection Frequency",
                                         show_threshold = FALSE,
                                         show_legend = FALSE,
                                         title = NULL,
                                         subtitle = NULL) {

  if (is.null(var_info)) {
    var_info <- .compute_mi_unified_var_order(
      object,
      top_n = top_n,
      select_by = sort_by,
      order_by = sort_by,
      freq_col = freq_col,
      require_nonzero_weight = FALSE
    )
  }

  if (length(var_info$var_order) == 0) {
    message("No pooled selection-frequency data available to plot")
    return(invisible(NULL))
  }

  weights_df <- extract_weights(object, direction = "both", by_group = TRUE)

  if (!freq_col %in% names(weights_df)) {
    fallback_col <- if ("mi_selection_freq" %in% names(weights_df)) "mi_selection_freq" else NULL
    if (is.null(fallback_col)) {
      weights_df[[freq_col]] <- 1
    } else {
      weights_df[[freq_col]] <- weights_df[[fallback_col]]
    }
  }

  plot_data <- weights_df[weights_df$variable %in% var_info$var_order, , drop = FALSE]

  if (nrow(plot_data) == 0) {
    message("No pooled selection-frequency data available to plot")
    return(invisible(NULL))
  }

  plot_data$selection_freq <- plot_data[[freq_col]]
  plot_data$freq_butterfly <- ifelse(
    plot_data$direction == "negative",
    -plot_data$selection_freq,
    plot_data$selection_freq
  )

  has_groups <- "group" %in% names(plot_data) && length(unique(var_info$df$group)) > 1

  if (has_groups) {
    group_order <- unique(var_info$df$group)
    plot_data$group <- factor(plot_data$group, levels = group_order)
    df_list <- lapply(split(plot_data, plot_data$group, drop = TRUE), function(g) {
      g_order <- var_info$var_order[var_info$var_order %in% g$variable]
      g$variable <- factor(g$variable, levels = rev(g_order))
      g
    })
    plot_data <- do.call(rbind, df_list)
  } else {
    plot_data$variable <- factor(plot_data$variable, levels = rev(var_info$var_order))
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$freq_butterfly))

  if (has_groups) {
    p <- p + 
      ggplot2::geom_bar(ggplot2::aes(fill = .data$direction), stat = "identity", width = 0.7) +
      ggplot2::facet_wrap(~ group, scales = "free_y", ncol = 1, strip.position = "right")
  } else {
    p <- p + 
      ggplot2::geom_bar(ggplot2::aes(fill = .data$direction), stat = "identity", width = 0.7)
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c("positive" = pos_color, "negative" = neg_color),
      labels = c("positive" = "Positive", "negative" = "Negative"),
      name = "Direction"
    ) +
    ggplot2::scale_y_continuous(labels = function(x) paste0(abs(x) * 100, "%")) +
    ggplot2::labs(
      x = NULL,
      y = y_label,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = if (show_legend) "bottom" else "none",
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )

  if (show_threshold) {
    p <- p + ggplot2::geom_hline(
      yintercept = c(-0.5, 0.5),
      linetype = "dashed",
      color = "gray60",
      linewidth = 0.3
    )
  }

  return(p)
}


#' Plot Pooled Validation Results (Forest Plot)
#' @keywords internal
.plot_pooled_validation <- function(object, pos_color = "#0072B2", neg_color = "#E69F00",
                                     base_size = 11, conf_level = 0.95,
                                     include_covariates = FALSE,
                                     show_significance = c("none", "asterisk", "pvalue")) {
  
  val <- .get_mi_pooled_inference(object)
  show_significance <- match.arg(show_significance)
  
  if (isTRUE(val$has_group_results) && !is.null(val$group_results)) {
    plot_data <- do.call(rbind, lapply(names(val$group_results), function(grp) {
      r <- val$group_results[[grp]]
      data.frame(
        term = paste0(grp, c(" (positive)", " (negative)")),
        group = grp,
        direction = c("Positive", "Negative"),
        estimate = c(r$positive$estimate, r$negative$estimate),
        se = c(r$positive$se, r$negative$se),
        p_value = c(r$positive$p_value, r$negative$p_value),
        fmi = c(r$positive$fmi, r$negative$fmi),
        df = c(r$positive$df, r$negative$df),
        type = "WQS Index",
        stringsAsFactors = FALSE
      )
    }))
  } else {
    plot_data <- data.frame(
      term = c("WQS Positive", "WQS Negative"),
      group = "Overall",
      direction = c("Positive", "Negative"),
      estimate = c(val$wqs_pos$estimate, val$wqs_neg$estimate),
      se = c(val$wqs_pos$se, val$wqs_neg$se),
      p_value = c(val$wqs_pos$p_value, val$wqs_neg$p_value),
      fmi = c(val$wqs_pos$fmi, val$wqs_neg$fmi),
      df = c(val$wqs_pos$df, val$wqs_neg$df),
      type = "WQS Index",
      stringsAsFactors = FALSE
    )
  }
  
  if (isTRUE(include_covariates)) {
    cov_df <- .pool_mi_validation_covariates(object, conf_level = conf_level)
    if (is.data.frame(cov_df) && nrow(cov_df) > 0) {
      shared_cols <- intersect(names(plot_data), names(cov_df))
      plot_data <- rbind(
        plot_data[, shared_cols, drop = FALSE],
        cov_df[, shared_cols, drop = FALSE]
      )
    }
  }
  
  na_effects <- plot_data[is.na(plot_data$estimate), , drop = FALSE]
  caption_text <- NULL
  
  if (nrow(na_effects) > 0) {
    pos_na_groups <- na_effects$group[na_effects$direction == "Positive"]
    neg_na_groups <- na_effects$group[na_effects$direction == "Negative"]
    
    caption_parts <- c()
    if (length(pos_na_groups) > 0) {
      caption_parts <- c(caption_parts,
        paste0("Positive effects not shown in ", paste(pos_na_groups, collapse = " and ")))
    }
    if (length(neg_na_groups) > 0) {
      caption_parts <- c(caption_parts,
        paste0("Negative effects not shown in ", paste(neg_na_groups, collapse = " and ")))
    }
    if (length(caption_parts) > 0) {
      caption_text <- paste0(paste(caption_parts, collapse = "; "),
                             "\n(coefficients = 0 due to regularization)")
    }
  }
  
  alpha <- 1 - conf_level
  plot_data$crit_value <- vapply(plot_data$df, function(df_val) {
    if (!is.na(df_val) && is.finite(df_val) && df_val > 0) {
      stats::qt(1 - alpha / 2, df = df_val)
    } else {
      stats::qnorm(1 - alpha / 2)
    }
  }, numeric(1))
  
  if (!"ci_lower" %in% names(plot_data)) plot_data$ci_lower <- NA_real_
  if (!"ci_upper" %in% names(plot_data)) plot_data$ci_upper <- NA_real_
  need_ci <- !is.finite(plot_data$ci_lower) | !is.finite(plot_data$ci_upper)
  plot_data$ci_lower[need_ci] <- plot_data$estimate[need_ci] - plot_data$crit_value[need_ci] * plot_data$se[need_ci]
  plot_data$ci_upper[need_ci] <- plot_data$estimate[need_ci] + plot_data$crit_value[need_ci] * plot_data$se[need_ci]
  
  if (show_significance == "asterisk") {
    plot_data$sig_label <- vapply(plot_data$p_value, function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      ""
    }, character(1))
  } else if (show_significance == "pvalue") {
    plot_data$sig_label <- vapply(plot_data$p_value, function(p) {
      if (is.na(p)) return("")
      format.pval(p, digits = 2)
    }, character(1))
  } else {
    plot_data$sig_label <- ""
  }
  
  plot_data$label <- ifelse(
    plot_data$type == "Covariate",
    plot_data$term,
    paste(plot_data$group, plot_data$direction, sep = "\n")
  )
  plot_data$label <- factor(plot_data$label, levels = rev(unique(plot_data$label)))
  plot_data$color <- ifelse(plot_data$type == "Covariate", "gray35",
                            ifelse(plot_data$direction == "Positive", pos_color, neg_color))
  
  plot_data_valid <- plot_data[!is.na(plot_data$estimate), , drop = FALSE]
  
  if (nrow(plot_data_valid) == 0) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "No mixture effects detected\n(all coefficients = 0)",
                        size = 4, color = "gray50") +
      ggplot2::theme_void()
    return(p)
  }
  
  p <- ggplot2::ggplot(plot_data_valid, ggplot2::aes(x = .data$estimate, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      orientation = "y",
      width = 0.2, linewidth = 0.8, color = "gray30"
    ) +
    ggplot2::geom_point(size = 3, color = plot_data_valid$color)
  
  if (any(plot_data_valid$sig_label != "", na.rm = TRUE)) {
    x_range <- max(plot_data_valid$ci_upper, na.rm = TRUE) - min(plot_data_valid$ci_lower, na.rm = TRUE)
    sig_offset <- if (is.finite(x_range) && x_range > 0) x_range * 0.08 else 0.1
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = .data$ci_upper + sig_offset, label = .data$sig_label),
      hjust = 0, vjust = 0.5, size = 4, fontface = "bold"
    )
  }
  
  p <- p +
    ggplot2::labs(
      x = paste0("Coefficient (", format(100 * conf_level, trim = TRUE, scientific = FALSE), "% CI)"),
      y = NULL,
      caption = caption_text
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(hjust = 0, size = base_size * 0.8,
                                            color = "gray40", face = "italic")
    )
  
  return(p)
}
