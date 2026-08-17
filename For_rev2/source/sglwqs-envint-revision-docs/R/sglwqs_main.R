#' Fit WQS Model with Sparse Group Lasso
#'
#' Fits a Weighted Quantile Sum (WQS) regression model using Sparse Group Lasso
#' for simultaneous variable selection and weight estimation. Supports bootstrap
#' aggregation for stable weights and train/validation split for exploratory
#' held-out conditional association summaries.
#'
#' @param X A data frame or matrix of exposure variables, or (when \code{data} is supplied) a selector such as variable names, numeric indices, or a logical column mask.
#' @param y A numeric outcome vector, or (when \code{data} is supplied) a single outcome variable name/index.
#' @param covariates Optional data frame or matrix of adjustment covariates.
#'   Factor and character variables are automatically dummy-coded. When \code{data}
#'   is supplied, this can also be specified via variable names, numeric indices,
#'   or a logical column mask.
#' @param groups Optional named list specifying variable groups for Group Lasso penalty.
#'   Example: list(PCBs = c("var1", "var2"), Phthalates = c("var3", "var4"))
#' @param n_quantiles Integer. Number of quantiles for transformation (default: 4).
#' @param family Character. Family for GLM ("gaussian" or "binomial", default: "gaussian").
#' @param lambda Character or numeric. Lambda selection method ("lambda.min" or "lambda.1se") 
#'   or a specific lambda value (default: "lambda.min").
#' @param lambda_path Optional numeric vector passed to the backend
#'   \code{cv.sparsegl(lambda = ...)}. This is an experimental escape hatch for
#'   survey-weighted selection diagnostics; \code{lambda} still controls the
#'   coefficient extraction point. When supplied, the backend uses this explicit
#'   sequence instead of generating an automatic path from \code{nlambda}.
#' @param nfolds Integer. Number of folds for cross-validation (default: 10).
#' @param penalize_covariates Logical. Whether to apply sparse penalty to covariates (default: FALSE).
#' @param group_by_compound Logical. Whether to group variables by chemical groups for
#'   Group Lasso penalty (default: TRUE when groups is specified).
#' @param group_structure Character. Fixed to "direction" - groups positive and negative 
#'   effects within each chemical group separately.
#' @param bootstrap Logical. Whether to use bootstrap aggregation for stable weights (default: FALSE).
#' @param n_boot Integer. Number of bootstrap iterations (default: 100).
#' @param boot_method Character. Bootstrap weighting method. \code{"auto"}
#'   uses ordinary row bootstrap without \code{survey_design} and survey
#'   replicate weights with \code{survey_design}; \code{"naive"} uses ordinary
#'   row bootstrap, and \code{"svrep"} uses replicate weights from
#'   \code{survey_design}.
#' @param svrep_type Character. Replicate design type passed to
#'   \code{survey::as.svrepdesign()} when \code{boot_method = "svrep"} and
#'   \code{survey_design} is not already a replicate design. \code{"auto"}
#'   uses \code{"bootstrap"} replicate weights; specify \code{"JK1"} or
#'   \code{"JKn"} explicitly when jackknife replicate weights are intended.
#' @param svrep_args Optional list of additional arguments passed to
#'   \code{survey::as.svrepdesign()}.
#' @param parallel Logical. Whether to use parallel processing for bootstrap (default: FALSE).
#'   Requires the 'future.apply' package. Set up parallel backend with future::plan() before calling.
#' @param keep_boot_matrices Logical. Whether to keep full bootstrap coefficient matrices
#'   in the output (default: FALSE). Setting to TRUE increases memory usage but allows
#'   custom analysis of bootstrap distributions.
#' @param checkpoint_dir Character. Directory to save intermediate bootstrap results 
#'   (default: NULL, no checkpointing). When specified, results are saved periodically,
#'   allowing recovery from crashes and reducing memory pressure in large-scale analyses.
#' @param checkpoint_interval Integer. Number of bootstrap iterations between checkpoints 
#'   (default: 50). Smaller values provide more frequent saves but may slow computation.
#' @param cleanup_checkpoint Logical. Whether to delete checkpoint files after successful
#'   completion (default: TRUE).
#' @param future_globals_max_size Numeric. Maximum size (in bytes) of global variables 
#'   allowed to be exported to parallel workers (default: NULL, auto-detect).
#'   When NULL and data size > 100MB, automatically calculates an appropriate size.
#'   For manual control, specify a value like \code{1000 * 1024^2} for 1GB.
#'   If the user has already set a larger value via \code{options()}, that value is respected.
#'   This option is only used when \code{parallel = TRUE}.
#' @param stratified_bootstrap Logical. Whether to use stratified bootstrap for binomial family 
#'   (default: TRUE). When TRUE and family = "binomial", bootstrap samples maintain the original 
#'   case/control ratio, which improves convergence stability for imbalanced binary outcomes.
#' @param validation Logical. Whether to split data for exploratory
#'   train/validation summaries (default: FALSE).
#' @param train_prop Numeric. Proportion of data for training (default: 0.6).
#' @param seed Integer. Random seed for reproducibility (default: NULL).
#' @param verbose Logical. Whether to show progress messages (default: TRUE).
#' @param obs_weights Optional numeric vector of observation weights used in
#'   the sparse-group selection loss.
#' @param quantile_weights Optional numeric vector of weights used to compute
#'   quantile cutpoints. When omitted, defaults to \code{survey_design}
#'   sampling weights when a survey design is supplied, then \code{obs_weights},
#'   otherwise unweighted.
#' @param refit Character. One of \code{"none"}, \code{"full"}, or
#'   \code{"validation"}.
#' @param refit_engine Character. One of \code{"glm"} or \code{"svyglm"}.
#' @param survey_design Optional pre-constructed \code{svydesign} or
#'   \code{svrepdesign} object used when \code{refit_engine = "svyglm"}.
#'   Downstream survey refit is supported; survey-weighted sparse-group
#'   selection uses the current \pkg{sparsegl} weighted backend and records
#'   \code{selection_diagnostics}, including all-zero exposure and
#'   lambda-path-boundary flags. These flags are backend diagnostics rather
#'   than fit failures.
#' @param analysis_id Optional analysis-row identifier used to verify alignment
#'   between the modeling data and \code{survey_design}. If \code{data} is
#'   supplied, a single column name may be used.
#' @param data Optional data frame or matrix used to resolve \code{X}, \code{y},
#'   \code{covariates}, \code{exposure_vars}, \code{outcome_var}, and
#'   \code{covariate_vars} by column name/index.
#' @param formula Optional formula of the form \code{y ~ x1 + x2 | z1 + z2}.
#'   Requires \code{data}; the left-hand side defines the outcome, the right-hand
#'   side before \code{|} defines exposures, and the optional part after \code{|}
#'   defines covariates.
#' @param exposure_vars Optional alias for exposure variable selection when
#'   \code{data} is supplied.
#' @param outcome_var Optional alias for outcome variable selection when
#'   \code{data} is supplied.
#' @param covariate_vars Optional alias for covariate selection when
#'   \code{data} is supplied.
#' @param minor_threshold Numeric (0-1). When groups are specified and validation/refit
#'   is performed, group-direction indices whose coefficient mass ratio
#'   (minor direction / major direction) falls below this threshold are excluded
#'   from the validation GLM to prevent suppressor effects. Set to 0 to disable.
#'   Default: 0.10.
#' @param ... Additional arguments passed to cv.sparsegl (e.g., \code{nlambda},
#'   \code{maxit}, \code{trace_it}).
#'
#' @return An object of class \code{"sglwqs"} with point estimates plus any
#'   requested inference layers. Common components include:
#'   \itemize{
#'   \item \code{pos_weights}, \code{neg_weights}: WQS weights.
#'   \item \code{cov_coef}: Penalized covariate coefficients from the sparse-group
#'     selection stage.
#'   \item \code{boot_info}: Present when \code{bootstrap = TRUE}; contains
#'     bootstrap summaries for exposure coefficients, covariates, and WQS
#'     index-sum summaries.
#'   \item \code{validation_info}: Present when \code{validation = TRUE}; contains
#'     downstream held-out GLM inference.
#'   \item \code{refit_info}: Present when \code{refit = "full"}; contains
#'     downstream full-data GLM inference.
#'   \item \code{diagnostics}: Fact-only diagnostic summaries returned by
#'     \code{compute_diagnostics()}.
#'   }
#'
#' @details
#' \strong{Bootstrap Aggregation (bootstrap = TRUE):}
#' Performs B bootstrap iterations to stabilize weight estimates. This is especially
#' recommended when the sample size is small relative to the number of exposures.
#' Returns selection frequencies and standard errors for coefficients.
#'
#' \strong{Parallel Processing (parallel = TRUE):}
#' When using bootstrap, parallel processing can significantly speed up computation.
#' Requires the 'future.apply' package. Set up parallel backend before calling:
#' \preformatted{
#'   library(future)
#'   plan(multisession, workers = 4)  # Use 4 cores
#'   fit <- sglwqs(..., bootstrap = TRUE, parallel = TRUE)
#'   plan(sequential)  # Reset to sequential
#' }
#'
#' \strong{Train/Validation Split (validation = TRUE):}
#' Splits data into training and validation sets. Weights are estimated on the training
#' set, then fixed weights are used to construct WQS indices for GLM on the validation set.
#' This provides exploratory validation-stage summaries conditional on
#' training-estimated indices. The validation GLM treats the constructed
#' indices as fixed and does not propagate uncertainty from regularized index
#' estimation, lambda selection, bootstrap aggregation, or minor-direction
#' filtering; p-values are therefore conditional summaries, not formal
#' post-selection inference.
#'
#' \strong{Two-Stage Inference Paths:}
#' \itemize{
#'   \item \code{bootstrap = FALSE, refit = "none"} returns point estimates only.
#'   \item \code{bootstrap = TRUE, refit = "none"} returns bootstrap summaries for
#'     weight stability, covariate coefficients, and WQS index-sum summaries.
#'   \item \code{refit = "full"} runs an in-sample downstream GLM after weights are
#'     fixed.
#'   \item \code{validation = TRUE} (equivalently \code{refit = "validation"})
#'     runs held-out downstream GLM inference.
#' }
#' Use \code{summary_inference()} or \code{plot_inference_results()} to access
#' whichever downstream or bootstrap-only inference source is active.
#'
#' \strong{Flexible Input Specification:}
#' In addition to the traditional \code{X}/\code{y}/\code{covariates} interface,
#' \code{sglwqs()} can also resolve variable names from a supplied \code{data}
#' object via \code{X}, \code{y}, \code{covariates},
#' \code{exposure_vars}/\code{outcome_var}/\code{covariate_vars}, or
#' \code{formula}. This mirrors the \code{sglwqs_mice()} interface for improved
#' consistency across complete-data and multiply-imputed workflows.
#'
#' \strong{Group Structure:}
#' \itemize{
#'   \item "direction": Positive and negative directions are separate groups.
#' }
#'
#' @examples
#' # Quick runnable example with built-in data
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:200, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:200],
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' summary(fit)
#'
#' \dontrun{
#' library(gWQS)
#' data("wqs_data")
#' 
#' PCBs <- names(wqs_data)[grepl("^LBX", names(wqs_data))]
#' phthalates <- names(wqs_data)[grepl("^URX", names(wqs_data))]
#' 
#' # Basic fit
#' fit1 <- sglwqs(
#'   X = wqs_data[, c(PCBs, phthalates)],
#'   y = wqs_data$y,
#'   covariates = wqs_data["sex"],
#'   groups = list(PCBs = PCBs, Phthalates = phthalates)
#' )
#' 
#' # With bootstrap for stable weights (parallel)
#' library(future)
#' plan(multisession, workers = 4)
#' fit2 <- sglwqs(
#'   X = wqs_data[, c(PCBs, phthalates)],
#'   y = wqs_data$y,
#'   covariates = wqs_data["sex"],
#'   groups = list(PCBs = PCBs, Phthalates = phthalates),
#'   bootstrap = TRUE,
#'   n_boot = 100,
#'   parallel = TRUE,
#'   seed = 123
#' )
#' plan(sequential)
#' 
#' # With validation for exploratory conditional summaries
#' fit3 <- sglwqs(
#'   X = wqs_data[, c(PCBs, phthalates)],
#'   y = wqs_data$y,
#'   covariates = wqs_data["sex"],
#'   groups = list(PCBs = PCBs, Phthalates = phthalates),
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   seed = 123
#' )
#' }
#'
#' @export
sglwqs <- function(X = NULL, y = NULL, covariates = NULL, groups = NULL, n_quantiles = 4, 
                   family = c("gaussian", "binomial"),
                   lambda = "lambda.min", nfolds = 10,
                   lambda_path = NULL,
                   penalize_covariates = FALSE, 
                   group_by_compound = NULL,
                   group_structure = "direction",
                   bootstrap = FALSE,
                   n_boot = 100,
                   boot_method = c("auto", "naive", "svrep"),
                   svrep_type = c("auto", "JK1", "JKn", "BRR", "Fay",
                                  "bootstrap", "subbootstrap", "mrbbootstrap"),
                   svrep_args = list(),
                   parallel = FALSE,
                   keep_boot_matrices = FALSE,
                   checkpoint_dir = NULL,
                   checkpoint_interval = 50,
                   cleanup_checkpoint = TRUE,
                   future_globals_max_size = NULL,
                   stratified_bootstrap = TRUE,
                   validation = FALSE,
                   train_prop = 0.6,
                   seed = NULL,
                   verbose = TRUE,
                   obs_weights = NULL,
                   quantile_weights = NULL,
                   refit = c("none", "full", "validation"),
                   refit_engine = c("glm", "svyglm"),
                   survey_design = NULL,
                   analysis_id = NULL,
                   minor_threshold = 0.10,
                   data = NULL,
                   formula = NULL,
                   exposure_vars = NULL,
                   outcome_var = NULL,
                   covariate_vars = NULL,
                   ...) {
  
  refit_missing <- missing(refit)
  refit_engine_missing <- missing(refit_engine)
  obs_weights_supplied <- !missing(obs_weights) && !is.null(obs_weights)
  quantile_weights_supplied <- !missing(quantile_weights) && !is.null(quantile_weights)
  family <- match.arg(family)
  lambda_path <- .validate_lambda_path(lambda_path)
  refit <- if (refit_missing) "none" else match.arg(refit)
  refit_requested <- refit
  if (validation && !refit_missing && !identical(refit_requested, "validation")) {
    stop(
      "`validation = TRUE` is only compatible with `refit = \"validation\"`.",
      call. = FALSE
    )
  }
  if (validation) {
    refit <- "validation"
  } else if (identical(refit, "validation")) {
    validation <- TRUE
  }
  boot_method <- match.arg(boot_method)
  svrep_type <- match.arg(svrep_type)
  refit_engine <- if (refit_engine_missing) {
    if (!is.null(survey_design) && !identical(refit, "none")) "svyglm" else "glm"
  } else {
    match.arg(refit_engine)
  }
  survey_mode <- !is.null(survey_design)
  survey_refit_mode <- identical(refit_engine, "svyglm")
  
  if (!is.numeric(minor_threshold) || length(minor_threshold) != 1 ||
      !is.finite(minor_threshold) || minor_threshold < 0 || minor_threshold > 1) {
    stop("`minor_threshold` must be a single numeric value between 0 and 1.", call. = FALSE)
  }

  if (survey_mode && !identical(refit, "none") && !survey_refit_mode) {
    stop(
      "`survey_design` with downstream refit requires `refit_engine = \"svyglm\"`.",
      call. = FALSE
    )
  }
  if (survey_refit_mode && is.null(survey_design)) {
    stop("`survey_design` must be provided when `refit_engine = \"svyglm\"`.", call. = FALSE)
  }
  if (survey_refit_mode && !identical(refit, "full") && !identical(refit, "validation")) {
    stop(
      "Survey refit is supported only for `refit = \"full\"` or `validation = TRUE`.",
      call. = FALSE
    )
  }
  if (survey_mode || survey_refit_mode) {
    if (!family %in% c("gaussian", "binomial")) {
      stop(
        "Survey refit is currently supported only for `family = \"gaussian\"` ",
        "or `family = \"binomial\"`.",
        call. = FALSE
      )
    }
  }
  if (bootstrap && identical(boot_method, "auto")) {
    boot_method <- if (!is.null(survey_design)) "svrep" else "naive"
  }
  if (bootstrap && identical(boot_method, "svrep") && is.null(survey_design)) {
    stop("`boot_method = \"svrep\"` requires `survey_design`.", call. = FALSE)
  }
  if (bootstrap && identical(boot_method, "naive") && !is.null(survey_design)) {
    warning(
      "`boot_method = \"naive\"` with `survey_design` ignores survey ",
      "cluster/strata structure; only design weights are reflected. Use ",
      "`boot_method = \"svrep\"` for replicate-weight survey bootstrap.",
      call. = FALSE
    )
  }
  if (bootstrap && identical(boot_method, "svrep") && identical(family, "binomial") &&
      isTRUE(stratified_bootstrap)) {
    message(
      "Outcome stratification (`stratified_bootstrap = TRUE`) is ignored when ",
      "`boot_method = \"svrep\"`; survey replicate weights define the resampling."
    )
    stratified_bootstrap <- FALSE
  }

  .advise_on_settings(
    bootstrap = bootstrap,
    validation = validation,
    refit = refit,
    is_mids = FALSE,
    verbose = verbose
  )
  
  # group_structure is fixed to "direction"
  if (!identical(group_structure, "direction")) {
    warning("group_structure = '", group_structure, "' is no longer supported; using 'direction'.",
            call. = FALSE)
    group_structure <- "direction"
  }

  data_frame <- NULL
  data_names <- NULL
  if (!is.null(data)) {
    if (!is.data.frame(data) && !is.matrix(data)) {
      stop("`data` must be a data.frame or matrix.")
    }
    data_frame <- as.data.frame(data)
    data_names <- names(data_frame)
  }

  if (!is.null(formula)) {
    if (is.null(data_frame)) {
      stop("`formula` requires `data`.")
    }
    formula_parts <- parse_sglwqs_formula(formula)
    if (is.null(exposure_vars) && is.null(X)) exposure_vars <- formula_parts$exposures
    if (is.null(outcome_var) && is.null(y)) outcome_var <- formula_parts$outcome
    if (is.null(covariate_vars) && is.null(covariates)) covariate_vars <- formula_parts$covariates
  }

  if (!is.null(exposure_vars)) {
    if (is.null(data_frame)) {
      stop("`exposure_vars` requires `data`.")
    }
    exposure_vars <- .resolve_data_var_spec(exposure_vars, data_names, "exposure_vars")
  }

  if (!is.null(outcome_var)) {
    if (is.null(data_frame)) {
      stop("`outcome_var` requires `data`.")
    }
    outcome_var <- .resolve_single_data_var(outcome_var, data_names, "outcome_var")
  }

  if (!is.null(covariate_vars)) {
    if (is.null(data_frame)) {
      stop("`covariate_vars` requires `data`.")
    }
    covariate_vars <- .resolve_data_var_spec(covariate_vars, data_names, "covariate_vars")
  }

  if (!is.null(analysis_id) && !is.null(data_frame) &&
      .is_single_data_var_selector(analysis_id, length(data_names))) {
    analysis_id_var <- .resolve_single_data_var(analysis_id, data_names, "analysis_id")
    analysis_id <- data_frame[[analysis_id_var]]
  }

  if (!is.null(X) && !is.null(data_frame) && .is_data_var_selector(X, length(data_names))) {
    X_alias <- .resolve_data_var_spec(X, data_names, "X")
    if (is.null(exposure_vars)) {
      exposure_vars <- X_alias
    } else if (!identical(exposure_vars, X_alias)) {
      stop("`X` and `exposure_vars` refer to different variables.")
    }
    X <- NULL
  }

  if (!is.null(y) && !is.null(data_frame) && .is_single_data_var_selector(y, length(data_names))) {
    y_alias <- .resolve_single_data_var(y, data_names, "y")
    if (is.null(outcome_var)) {
      outcome_var <- y_alias
    } else if (!identical(outcome_var, y_alias)) {
      stop("`y` and `outcome_var` refer to different variables.")
    }
    y <- NULL
  }

  if (!is.null(covariates) && !is.null(data_frame) && .is_data_var_selector(covariates, length(data_names))) {
    cov_alias <- .resolve_data_var_spec(covariates, data_names, "covariates")
    if (is.null(covariate_vars)) {
      covariate_vars <- cov_alias
    } else if (!identical(covariate_vars, cov_alias)) {
      stop("`covariates` and `covariate_vars` refer to different variables.")
    }
    covariates <- NULL
  }

  if (is.null(X) && !is.null(exposure_vars)) {
    X <- data_frame[, exposure_vars, drop = FALSE]
  }

  if (is.null(y) && !is.null(outcome_var)) {
    y <- data_frame[[outcome_var]]
  }

  if (is.null(covariates) && !is.null(covariate_vars)) {
    covariates <- data_frame[, covariate_vars, drop = FALSE]
  }

  if (is.null(X) || is.null(y)) {
    stop("Specify exposures/outcome via `X` + `y`, `exposure_vars` + `outcome_var`, or `formula` (with `data`).")
  }

  # Input validation
  if (is.data.frame(X)) X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix or data.frame.", call. = FALSE)
  }
  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector.", call. = FALSE)
  }
  if (length(y) != nrow(X)) {
    stop("`y` (length ", length(y), ") must have the same length as `nrow(X)` (", nrow(X), ").",
         call. = FALSE)
  }
  if (any(is.infinite(y))) {
    stop("`y` contains Inf or -Inf values.", call. = FALSE)
  }
  if (all(is.na(y))) {
    stop("`y` is entirely NA.", call. = FALSE)
  }
  if (family == "binomial") {
    y_unique <- unique(y[!is.na(y)])
    if (!all(y_unique %in% c(0, 1))) {
      stop("`y` must be binary (0/1) when `family = \"binomial\"`.", call. = FALSE)
    }
  }
  if (!is.numeric(train_prop) || length(train_prop) != 1L ||
      !is.finite(train_prop) || train_prop <= 0 || train_prop >= 1) {
    stop("`train_prop` must be a single numeric value strictly between 0 and 1.",
         call. = FALSE)
  }

  # Set future.globals.maxSize (only for parallel processing)
  if (parallel) {
    original_max_size <- getOption("future.globals.maxSize")
    
    # Calculate recommended size based on data size
    data_size_mb <- object.size(X) / 1024^2
    if (!is.null(covariates)) {
      data_size_mb <- data_size_mb + object.size(covariates) / 1024^2
    }
    # Safety margin (2x) + base 500MB
    auto_recommended_size <- max(500, ceiling(data_size_mb * 2 + 500)) * 1024^2
    
    # Logic for determining future_globals_max_size
    if (is.null(future_globals_max_size)) {
      # If NULL: auto-adjust for large data
      if (data_size_mb > 100) {  # Data exceeding 100MB
        effective_max_size <- auto_recommended_size
        if (verbose) {
          message("Large dataset detected (", round(data_size_mb), " MB). ",
                  "Auto-setting future.globals.maxSize to ", 
                  round(effective_max_size / 1024^2), " MB")
        }
        options(future.globals.maxSize = effective_max_size)
      }
      # Use future's default (500MB) for data <= 100MB
    } else {
      # If explicitly specified
      effective_max_size <- future_globals_max_size
      
      # Respect user's existing larger setting
      if (!is.null(original_max_size) && original_max_size > future_globals_max_size) {
        effective_max_size <- original_max_size
        if (verbose) {
          message("Using existing future.globals.maxSize: ", 
                  round(original_max_size / 1024^2), " MB (larger than requested)")
        }
      } else {
        options(future.globals.maxSize = future_globals_max_size)
        if (verbose) {
          message("Setting future.globals.maxSize to ", 
                  round(future_globals_max_size / 1024^2), " MB")
        }
      }
    }
    
    # Restore original setting on function exit
    on.exit({
      if (!is.null(original_max_size)) {
        options(future.globals.maxSize = original_max_size)
      } else {
        options(future.globals.maxSize = NULL)
      }
    }, add = TRUE)
  }
  
  # Set seed
  if (!is.null(seed)) set.seed(seed)

  # Default setting for group_by_compound
  if (is.null(group_by_compound)) {
    group_by_compound <- !is.null(groups)
  }

  # Get variable names
  if (is.null(colnames(X))) {
    var_names <- paste0("X", seq_len(ncol(X)))
  } else {
    var_names <- colnames(X)
  }

  # Validate variable and group names (warn about special characters)
  if (verbose) {
    validate_names(var_names, groups)
  }

  # Validate/complete groups and ensure group_by_compound consistency

  if (!isTRUE(group_by_compound) && !is.null(groups)) {
    warning("`groups` is ignored when `group_by_compound = FALSE`.", call. = FALSE)
    groups_active <- NULL
  } else {
    groups_active <- .validate_and_complete_groups(var_names, groups)
  }
  # Use groups_active from here on
  groups <- groups_active
  
  # Process covariates (automatic dummy coding)
  cov_processed <- process_covariates(covariates)
  cov_matrix <- cov_processed$matrix
  cov_names <- cov_processed$names
  
  n <- nrow(X)
  p <- ncol(X)
  q <- if (is.null(cov_matrix)) 0 else ncol(cov_matrix)
  analysis_id <- .validate_analysis_id(analysis_id, n)
  obs_weights <- .validate_obs_weights(obs_weights, n)
  analysis_weights <- NULL
  quantile_weights_source <- "unweighted"
  survey_info <- NULL

  if (survey_mode) {
    .validate_survey_design(survey_design)
    alignment_df <- data.frame(.row = seq_len(n))
    rownames(alignment_df) <- rownames(X)
    if (is.null(analysis_id) && !.has_automatic_rownames(alignment_df)) {
      analysis_id <- rownames(alignment_df)
    }
    .check_survey_alignment(
      survey_design = survey_design,
      glm_data = alignment_df,
      analysis_id = analysis_id,
      context = "analysis"
    )
    design_weights <- .extract_analysis_weights(survey_design)
    if (is.null(design_weights)) {
      stop(
        "Could not extract sampling weights from `survey_design`. ",
        "Build the design with explicit sampling weights before using survey mode.",
        call. = FALSE
      )
    }
    design_weights <- .validate_obs_weights(design_weights, n)
    analysis_weights <- design_weights
    if (obs_weights_supplied && !isTRUE(all.equal(obs_weights, design_weights))) {
      warning(
        "`survey_design` sampling weights are used for survey-mode selection; ",
        "explicit `obs_weights` are ignored.",
        call. = FALSE
      )
    }
    obs_weights <- design_weights / mean(design_weights, na.rm = TRUE)
    survey_info <- list(
      refit_engine = refit_engine,
      degf = .survey_design_degf(survey_design),
      weights_normalized = TRUE
    )
  }

  if (quantile_weights_supplied) {
    quantile_weights <- .validate_quantile_weights(quantile_weights, n)
    quantile_weights_source <- "explicit"
  } else if (survey_mode) {
    quantile_weights <- analysis_weights
    quantile_weights_source <- "survey_design"
  } else if (!is.null(obs_weights)) {
    quantile_weights <- obs_weights
    quantile_weights_source <- "obs_weights"
  }
  
  # ----- Train/Validation Split -----
  if (validation) {
    if (verbose) message("Splitting data: ", round(train_prop * 100), "% training, ", 
                          round((1 - train_prop) * 100), "% validation")
    
    # Stratified split: maintain class ratio for binomial
    if (identical(family, "binomial")) {
      y_non_na <- y[!is.na(y)]
      y_levels <- sort(unique(y_non_na))
      if (length(y_levels) == 2) {
        train_idx <- integer(0)
        for (lvl in y_levels) {
          idx_lvl <- which(y == lvl)
          n_lvl <- length(idx_lvl)
          n_train_lvl <- floor(n_lvl * train_prop)
          if (n_lvl >= 2) {
            n_train_lvl <- min(max(n_train_lvl, 1L), n_lvl - 1L)
          }
          train_idx <- c(train_idx, sample(idx_lvl, n_train_lvl))
        }
        train_idx <- sort(train_idx)
      } else {
        train_idx <- sample(n, floor(n * train_prop))
      }
    } else {
      train_idx <- sample(n, floor(n * train_prop))
    }
    if (survey_mode) {
      train_idx <- sort(train_idx)
    }
    val_idx <- setdiff(seq_len(n), train_idx)
    
    # Quantile transform: learn breaks on train data, apply to val data (prevent leakage)
    X_train_raw <- X[train_idx, , drop = FALSE]
    X_val_raw <- X[val_idx, , drop = FALSE]
    quantile_weights_train <- if (!is.null(quantile_weights)) quantile_weights[train_idx] else NULL
    obs_weights_train <- if (!is.null(obs_weights)) obs_weights[train_idx] else NULL
    
    train_qt_result <- quantile_transform(
      X_train_raw,
      n_quantiles = n_quantiles,
      var_names = var_names,
      weights = quantile_weights_train
    )
    X_train <- train_qt_result$q
    if (!is.null(rownames(X_train_raw))) {
      rownames(X_train) <- rownames(X_train_raw)
    }
    train_breaks <- train_qt_result$breaks
    
    # Validation: apply breaks learned on train data
    val_qt_result <- quantile_transform(X_val_raw, n_quantiles = n_quantiles, 
                                         var_names = var_names, breaks_list = train_breaks)
    X_val <- val_qt_result$q
    if (!is.null(rownames(X_val_raw))) {
      rownames(X_val) <- rownames(X_val_raw)
    }
    
    y_train <- y[train_idx]
    y_val <- y[val_idx]
    cov_train <- if (!is.null(cov_matrix)) cov_matrix[train_idx, , drop = FALSE] else NULL
    cov_val <- if (!is.null(cov_matrix)) cov_matrix[val_idx, , drop = FALSE] else NULL
    obs_weights_val <- if (!is.null(obs_weights)) obs_weights[val_idx] else NULL
    train_survey_design <- if (survey_mode) {
      .subset_survey_design(survey_design, train_idx, context = "training")
    } else {
      NULL
    }
    
    # Quantile transform all data (for output) - apply train breaks to full data
    full_qt_result <- quantile_transform(X, n_quantiles = n_quantiles, 
                                          var_names = var_names, breaks_list = train_breaks)
    X_quantile <- full_qt_result$q
    if (!is.null(rownames(X))) {
      rownames(X_quantile) <- rownames(X)
    }
    quantile_breaks <- train_breaks
    
  } else {
    # Quantile transform (computed on full data)
    qt_result <- quantile_transform(
      X,
      n_quantiles = n_quantiles,
      var_names = var_names,
      weights = quantile_weights
    )
    X_quantile <- qt_result$q
    if (!is.null(rownames(X))) {
      rownames(X_quantile) <- rownames(X)
    }
    quantile_breaks <- qt_result$breaks
    
    X_train <- X_quantile
    y_train <- y
    cov_train <- cov_matrix
    obs_weights_train <- obs_weights
    train_survey_design <- survey_design
    
    X_val <- NULL
    y_val <- NULL
    cov_val <- NULL
  }
  
  # ----- Weight Estimation -----
  if (bootstrap) {
    # Bootstrap aggregation
    boot_result <- bootstrap_sgl(
      X_quantile = X_train,
      y = y_train,
      cov_matrix = cov_train,
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
      n_boot = n_boot,
      seed = NULL,  # Seed already set
      verbose = verbose,
      parallel = parallel,
      checkpoint_dir = checkpoint_dir,
      checkpoint_interval = checkpoint_interval,
      cleanup_checkpoint = cleanup_checkpoint,
      stratified = stratified_bootstrap,
      obs_weights = obs_weights_train,
      boot_method = boot_method,
      survey_design = train_survey_design,
      svrep_type = svrep_type,
      svrep_args = svrep_args,
      ...
    )
    
    pos_coef <- boot_result$mean_pos_coef
    neg_coef <- boot_result$mean_neg_coef
    
    boot_info <- list(
      mean_pos_coef = boot_result$mean_pos_coef,
      mean_neg_coef = boot_result$mean_neg_coef,
      se_pos_coef = boot_result$se_pos_coef,
      se_neg_coef = boot_result$se_neg_coef,
      selection_freq_pos = boot_result$selection_freq_pos,
      selection_freq_neg = boot_result$selection_freq_neg,
      mean_index_sum_pos = boot_result$mean_index_sum_pos,
      mean_index_sum_neg = boot_result$mean_index_sum_neg,
      se_index_sum_pos = boot_result$se_index_sum_pos,
      se_index_sum_neg = boot_result$se_index_sum_neg,
      mean_index_sum_by_group_pos = boot_result$mean_index_sum_by_group_pos,
      mean_index_sum_by_group_neg = boot_result$mean_index_sum_by_group_neg,
      se_index_sum_by_group_pos = boot_result$se_index_sum_by_group_pos,
      se_index_sum_by_group_neg = boot_result$se_index_sum_by_group_neg,
      mean_cov_coef = boot_result$mean_cov_coef,
      se_cov_coef = boot_result$se_cov_coef,
      ci_lower_cov = boot_result$ci_lower_cov,
      ci_upper_cov = boot_result$ci_upper_cov,
      boot_success = boot_result$boot_success,
      boot_error_msg = boot_result$boot_error_msg,
      boot_error_class = boot_result$boot_error_class,
      boot_error_counts = boot_result$boot_error_counts,
      parallel_batch_errors = boot_result$parallel_batch_errors,
      n_parallel_batch_failures = boot_result$n_parallel_batch_failures,
      n_successful = boot_result$n_successful,
      n_failed = boot_result$n_failed,
      method = boot_result$method,
      svrep_type_used = boot_result$svrep_type_used,
      svrep_args = boot_result$svrep_args,
      bootstrap_design = boot_result$bootstrap_design,
      survey_design_ignored = boot_result$survey_design_ignored,
      variance_scale = boot_result$variance_scale,
      variance_rscales = boot_result$variance_rscales,
      variance_df = boot_result$variance_df,
      svrep_center_pos_coef = boot_result$svrep_center_pos_coef,
      svrep_center_neg_coef = boot_result$svrep_center_neg_coef,
      svrep_center_cov_coef = boot_result$svrep_center_cov_coef,
      svrep_center_index_sum_pos = boot_result$svrep_center_index_sum_pos,
      svrep_center_index_sum_neg = boot_result$svrep_center_index_sum_neg,
      svrep_center_index_sum_by_group_pos = boot_result$svrep_center_index_sum_by_group_pos,
      svrep_center_index_sum_by_group_neg = boot_result$svrep_center_index_sum_by_group_neg,
      weights_normalized = boot_result$weights_normalized,
      n_boot_actual = boot_result$n_boot_actual,
      n_boot_requested = boot_result$n_boot_requested
    )

    # Memory-saving option: do not store matrices when keep_boot_matrices=FALSE
    if (keep_boot_matrices) {
      boot_info$boot_pos_coef <- boot_result$boot_pos_coef
      boot_info$boot_neg_coef <- boot_result$boot_neg_coef
      boot_info$boot_cov_coef <- boot_result$boot_cov_coef
      boot_info$boot_pos_index_sum <- boot_result$boot_pos_index_sum
      boot_info$boot_neg_index_sum <- boot_result$boot_neg_index_sum
      boot_info$boot_pos_index_sum_by_group <- boot_result$boot_pos_index_sum_by_group
      boot_info$boot_neg_index_sum_by_group <- boot_result$boot_neg_index_sum_by_group
    }
    
    # Final fit on full data (for cv.sparsegl object)
    fit_result <- fit_sgl_core(
      X_quantile = X_train,
      y = y_train,
      cov_matrix = cov_train,
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
      obs_weights = obs_weights_train,
      ...
    )
    
    cov_coef <- fit_result$cov_coef
    intercept <- fit_result$intercept
    fit <- fit_result$fit
    
  } else {
    # Single fit
    fit_result <- fit_sgl_core(
      X_quantile = X_train,
      y = y_train,
      cov_matrix = cov_train,
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
      obs_weights = obs_weights_train,
      ...
    )
    
    pos_coef <- fit_result$pos_coef
    neg_coef <- fit_result$neg_coef
    cov_coef <- fit_result$cov_coef
    intercept <- fit_result$intercept
    fit <- fit_result$fit
    boot_info <- NULL
  }
  
  # ----- Calculate Weights -----
  weight_result <- calculate_weights(
    pos_coef = pos_coef,
    neg_coef = neg_coef,
    var_names = var_names,
    groups = groups,
    handle_collinearity = "net"
  )
  selection_diagnostics <- fit_result$selection_diagnostics
  if (!is.null(selection_diagnostics)) {
    selection_diagnostics$positive_weight_sum <- sum(weight_result$pos_weights, na.rm = TRUE)
    selection_diagnostics$negative_weight_sum <- sum(weight_result$neg_weights, na.rm = TRUE)
    selection_diagnostics$weighted_selection <- !is.null(obs_weights_train)
    selection_diagnostics$survey_mode <- isTRUE(survey_mode)
    if ((isTRUE(survey_mode) || !is.null(obs_weights_train)) &&
        isTRUE(selection_diagnostics$all_zero_exposure)) {
      warning(
        "Weighted SGL selection returned all-zero exposure coefficients. ",
        "Downstream survey refit may therefore have no retained WQS exposure index. ",
        "Inspect `fit$selection_diagnostics` for the lambda path and backend diagnostics.",
        call. = FALSE
      )
    } else if ((isTRUE(survey_mode) || !is.null(obs_weights_train)) &&
               isTRUE(selection_diagnostics$selected_lambda_at_path_boundary)) {
      warning(
        "Weighted SGL selected lambda at the path boundary. ",
        "The lambda path may be too short or overly sparse for the weighted selection problem. ",
        "Inspect `fit$selection_diagnostics` and consider an explicit `lambda_path` sensitivity check.",
        call. = FALSE
      )
    }
  }

  # In bootstrap mode, overwrite with "mean of per-draw weights" (manuscript estimator)

  if (bootstrap && !identical(boot_result$method, "svrep") &&
      !is.null(boot_result$boot_pos_coef)) {
    agg_weights <- .aggregate_bootstrap_weights(
      boot_pos_coef = boot_result$boot_pos_coef,
      boot_neg_coef = boot_result$boot_neg_coef,
      boot_success = boot_result$boot_success,
      var_names = var_names,
      groups = groups,
      handle_collinearity = "net"
    )
    weight_result$pos_weights <- agg_weights$pos_weights
    weight_result$neg_weights <- agg_weights$neg_weights
    weight_result$pos_index_sum <- agg_weights$pos_index_sum
    weight_result$neg_index_sum <- agg_weights$neg_index_sum
    weight_result$pos_index_sum_by_group <- agg_weights$pos_index_sum_by_group
    weight_result$neg_index_sum_by_group <- agg_weights$neg_index_sum_by_group
  }

  # ----- Minor Direction Exclusion -----
  excluded_directions <- list()

  if (!is.null(groups) && minor_threshold > 0 && !identical(refit, "none")) {
    for (grp_name in names(groups)) {
      pos_mass <- weight_result$pos_index_sum_by_group[[grp_name]]
      neg_mass <- weight_result$neg_index_sum_by_group[[grp_name]]

      if (is.null(pos_mass)) pos_mass <- 0
      if (is.null(neg_mass)) neg_mass <- 0

      major <- max(pos_mass, neg_mass)
      minor <- min(pos_mass, neg_mass)

      if (major > 0 && (minor / major) < minor_threshold) {
        if (pos_mass < neg_mass) {
          excluded_directions[[grp_name]] <- "positive"
        } else {
          excluded_directions[[grp_name]] <- "negative"
        }

        if (verbose) {
          cat(sprintf("  Group '%s': excluding %s direction (mass ratio = %.1f%%)\n",
                      grp_name, excluded_directions[[grp_name]],
                      minor / major * 100))
        }
      }
    }
  }

  # ----- Refit / Inference -----
  refit_info <- NULL
  validation_info <- NULL
  
  if (identical(refit, "validation")) {
    validation_survey_design <- if (survey_mode) {
      .subset_survey_design(survey_design, val_idx, context = "validation")
    } else {
      NULL
    }
    val_result <- validation_glm(
      X_quantile_val = X_val,
      y_val = y_val,
      cov_matrix_val = cov_val,
      pos_weights = weight_result$pos_weights,
      neg_weights = weight_result$neg_weights,
      family = family,
      groups = groups,
      group_inference = !is.null(groups),
      obs_weights = obs_weights_val,
      excluded_directions = excluded_directions,
      engine = refit_engine,
      survey_design = validation_survey_design,
      analysis_id = if (!is.null(analysis_id)) analysis_id[val_idx] else NULL
    )
    
    refit_info <- list(
      refit_fit = val_result$refit_fit,
      refit_summary = val_result$refit_summary,
      coef_table = val_result$coef_table,
      p_col = val_result$p_col,
      wqs_pos_estimate = val_result$wqs_pos_estimate,
      wqs_pos_se = val_result$wqs_pos_se,
      wqs_pos_pvalue = val_result$wqs_pos_pvalue,
      wqs_neg_estimate = val_result$wqs_neg_estimate,
      wqs_neg_se = val_result$wqs_neg_se,
      wqs_neg_pvalue = val_result$wqs_neg_pvalue,
      group_results = val_result$group_results,
      group_inference = val_result$group_inference,
      wqs_indices = val_result$wqs_indices,
      formula = val_result$formula,
      engine = val_result$engine,
      curve = val_result$curve,
      n_train = length(train_idx),
      n_val = length(val_idx),
      n_obs = length(val_idx),
      glm_fit = val_result$glm_fit,
      glm_summary = val_result$glm_summary
    )
    validation_info <- refit_info
  } else if (identical(refit, "full")) {
    refit_result <- refit_model(
      X_quantile = X_quantile,
      y = y,
      cov_matrix = cov_matrix,
      pos_weights = weight_result$pos_weights,
      neg_weights = weight_result$neg_weights,
      family = family,
      groups = groups,
      group_inference = !is.null(groups),
      engine = refit_engine,
      survey_design = survey_design,
      analysis_id = analysis_id,
      obs_weights = obs_weights,
      excluded_directions = excluded_directions
    )
    
    refit_info <- c(
      refit_result,
      list(
        n_train = NA_integer_,
        n_val = NA_integer_
      )
    )
  }
  
  # ----- Build Result Object -----
  result <- list(
    # Model fit
    fit = fit,
    coefficients = fit_result$coefficients,
    
    # Weights and coefficients
    pos_weights = weight_result$pos_weights,
    neg_weights = weight_result$neg_weights,
    pos_coef = weight_result$pos_coef,
    neg_coef = weight_result$neg_coef,
    
    # Index sums
    pos_index_sum = weight_result$pos_index_sum,
    neg_index_sum = weight_result$neg_index_sum,
    pos_index_sum_by_group = weight_result$pos_index_sum_by_group,
    neg_index_sum_by_group = weight_result$neg_index_sum_by_group,
    
    # Covariate coefficients
    cov_coef = cov_coef,
    intercept = intercept,
    
    # Variable info
    var_names = var_names,
    cov_names = cov_names,
    exposure_vars = var_names,
    outcome_var = outcome_var,
    covariate_vars = cov_processed$original_names,
    
    # Group info
    groups = groups,
    group_by_compound = group_by_compound,
    group_structure = if (group_by_compound) group_structure else NA_character_,
    
    # Model settings
    lambda = lambda,
    lambda_path = lambda_path,
    n_quantiles = n_quantiles,
    family = family,
    refit = refit,
    refit_engine = refit_engine,
    survey_mode = survey_mode,
    survey_info = survey_info,
    boot_method = boot_method,
    svrep_type = svrep_type,
    svrep_args = svrep_args,
    
    # Data info
    X_quantile = X_quantile,
    cov_matrix = cov_matrix,
    y = y,
    quantile_breaks = quantile_breaks,
    n = n,
    p = p,
    
    # Bootstrap info
    bootstrap = bootstrap,
    boot_info = boot_info,
    
    # Validation info
    validation = validation,
    validation_info = validation_info,
    refit_info = refit_info,
    obs_weights = obs_weights,
    analysis_id = analysis_id,
    analysis_weights = analysis_weights,
    quantile_weights = quantile_weights,
    quantile_weights_source = quantile_weights_source,
    selection_diagnostics = selection_diagnostics,
    
    # Minor direction exclusion
    minor_threshold = minor_threshold,
    excluded_directions = excluded_directions,

    # Warnings
    collinearity_warning = weight_result$collinearity_warning
  )

  class(result) <- "sglwqs"
  result$diagnostics <- compute_diagnostics(result)
  .emit_diagnostics(result$diagnostics, verbose = verbose)
  return(result)
}


#' @keywords internal
.is_data_var_selector <- function(spec, n_data_cols) {
  if (is.character(spec)) {
    return(TRUE)
  }

  if (is.numeric(spec)) {
    return(length(spec) >= 1 && all(is.finite(spec)) &&
             all(spec == as.integer(spec)) && all(spec >= 1) && all(spec <= n_data_cols))
  }

  if (is.logical(spec)) {
    return(length(spec) == n_data_cols)
  }

  FALSE
}

#' @keywords internal
.is_single_data_var_selector <- function(spec, n_data_cols) {
  if (is.character(spec)) {
    return(length(spec) == 1)
  }

  if (is.numeric(spec)) {
    return(length(spec) == 1 && is.finite(spec) && spec == as.integer(spec) &&
             spec >= 1 && spec <= n_data_cols)
  }

  if (is.logical(spec)) {
    return(length(spec) == n_data_cols)
  }

  FALSE
}

#' @keywords internal
.resolve_data_var_spec <- function(spec, data_names, arg_name) {
  if (is.null(spec)) {
    return(NULL)
  }

  if (is.character(spec)) {
    vars <- spec
  } else if (is.numeric(spec)) {
    if (any(!is.finite(spec)) || any(spec < 1) || any(spec > length(data_names)) ||
        any(spec != as.integer(spec))) {
      stop(arg_name, " must contain valid column indices for `data`.")
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
    stop(arg_name, " contains variables not found in `data`: ", paste(missing_vars, collapse = ", "))
  }

  vars
}

#' @keywords internal
.resolve_single_data_var <- function(spec, data_names, arg_name) {
  vars <- .resolve_data_var_spec(spec, data_names, arg_name)
  if (length(vars) != 1) {
    stop(arg_name, " must identify exactly one variable.")
  }
  vars
}


#' @export
print.sglwqs <- function(x, ...) {
  refit_info <- .get_active_refit_info(x)
  
  cat("SGL-WQS Model (Sparse Group Lasso Weighted Quantile Sum)\n")
  cat("========================================================\n\n")
  cat("Number of observations:", x$n, "\n")
  cat("Number of exposures:", x$p, "\n")
  cat("Number of quantiles:", x$n_quantiles, "\n")
  cat("Family:", x$family, "\n")
  cat("Lambda:", x$lambda, "\n")
  
  if (!is.null(x$groups)) {
    cat("Chemical groups:", length(x$groups), "\n")
    cat("Group-based penalty:", ifelse(x$group_by_compound, "Yes", "No"), "\n")
    if (x$group_by_compound && !is.na(x$group_structure)) {
      cat("Group structure:", x$group_structure, "\n")
    }
  }
  
  # Bootstrap/Validation information
  if (x$bootstrap) {
    cat("\nBootstrap: Yes (", x$boot_info$n_successful, " iterations)\n", sep = "")
  }
  if (!is.null(refit_info) && identical(attr(refit_info, "source"), "validation_info")) {
    cat("Validation: Yes (train:", refit_info$n_train,
        ", val:", refit_info$n_val, ")\n", sep = "")
  } else if (!is.null(refit_info)) {
    cat("Refit: ", x$refit, " (engine: ", refit_info$engine, ")\n", sep = "")
  }
  
  # Number of selected variables (sparsity indicator)
  n_pos_selected <- sum(x$pos_weights > 0)
  n_neg_selected <- sum(x$neg_weights > 0)
  cat("\nSelected variables:\n")
  cat("  Positive direction:", n_pos_selected, "/", x$p, 
      "(", round(n_pos_selected / x$p * 100, 1), "%)\n")
  cat("  Negative direction:", n_neg_selected, "/", x$p, 
      "(", round(n_neg_selected / x$p * 100, 1), "%)\n")
  
  cat("\nPositive Index Sum:", round(x$pos_index_sum, 4), "\n")
  cat("Negative Index Sum:", round(x$neg_index_sum, 4), "\n")
  
  if (!is.null(refit_info) && is.null(x$groups)) {
    header <- if (identical(attr(refit_info, "source"), "validation_info")) {
      "Validation Inference:"
    } else {
      "Refit Inference:"
    }
    cat("\n", header, "\n", sep = "")
    cat("  WQS Positive: estimate =", round(refit_info$wqs_pos_estimate, 4),
        ", SE =", round(refit_info$wqs_pos_se, 4),
        ", p =", format.pval(refit_info$wqs_pos_pvalue, digits = 3), "\n")
    cat("  WQS Negative: estimate =", round(refit_info$wqs_neg_estimate, 4),
        ", SE =", round(refit_info$wqs_neg_se, 4),
        ", p =", format.pval(refit_info$wqs_neg_pvalue, digits = 3), "\n")
  }
  
  # Excluded directions
  if (length(x$excluded_directions) > 0) {
    cat("\nExcluded directions (minor/major < ", x$minor_threshold * 100, "%):\n", sep = "")
    for (grp in names(x$excluded_directions)) {
      cat("  ", grp, ": ", x$excluded_directions[[grp]], " direction\n", sep = "")
    }
  }

  # Collinearity warning
  if (x$collinearity_warning) {
    cat("\nNote: Collinearity was detected and handled (see warnings).\n")
  }

  if (!is.null(x$diagnostics)) {
    cat("\n")
    print(x$diagnostics)
    cat("See `compute_diagnostics(fit)` for access to these values programmatically.\n")
  }

  invisible(x)
}


#' @export
summary.sglwqs <- function(object, ...) {
  refit_info <- .get_active_refit_info(object)
  refit_coef_table <- if (!is.null(refit_info)) refit_info$coef_table else NULL
  refit_p_col <- if (!is.null(refit_info)) refit_info$p_col else NULL
  
  cat("\nCall:\n")
  cat("sglwqs(n =", object$n, ", p =", object$p, ", family =", object$family, ")\n\n")
  
  # ----- Coefficient Table (GLM Style) -----
  cat("Coefficients:\n")
  
  # Build coefficient table
  coef_rows <- list()
  
  if (!is.null(refit_coef_table) && "(Intercept)" %in% rownames(refit_coef_table)) {
    intercept_est <- refit_coef_table["(Intercept)", "Estimate"]
    intercept_se <- refit_coef_table["(Intercept)", "Std. Error"]
    intercept_p <- refit_coef_table["(Intercept)", refit_p_col]
  } else {
    intercept_est <- object$intercept
    intercept_se <- NA_real_
    intercept_p <- NA_real_
  }
  
  coef_rows[[1]] <- data.frame(
    Term = "(Intercept)",
    Estimate = intercept_est,
    `Std.Error` = intercept_se,
    `z/t value` = intercept_est / intercept_se,
    `Pr(>|z|)` = intercept_p,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # If refit results exist
  if (!is.null(refit_info)) {
    val_info <- refit_info
    
    # If group-specific results exist
    if (!is.null(val_info$group_results) && length(val_info$group_results) > 0) {
      for (grp_name in names(val_info$group_results)) {
        grp_res <- val_info$group_results[[grp_name]]

        # Positive
        pos_label <- paste0(grp_name, " (positive)")
        if (isTRUE(grp_res$pos_excluded)) pos_label <- paste0(pos_label, " [excluded]")
        coef_rows[[length(coef_rows) + 1]] <- data.frame(
          Term = pos_label,
          Estimate = grp_res$pos_estimate,
          `Std.Error` = grp_res$pos_se,
          `z/t value` = grp_res$pos_estimate / grp_res$pos_se,
          `Pr(>|z|)` = grp_res$pos_pvalue,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        # Negative
        neg_label <- paste0(grp_name, " (negative)")
        if (isTRUE(grp_res$neg_excluded)) neg_label <- paste0(neg_label, " [excluded]")
        coef_rows[[length(coef_rows) + 1]] <- data.frame(
          Term = neg_label,
          Estimate = grp_res$neg_estimate,
          `Std.Error` = grp_res$neg_se,
          `z/t value` = grp_res$neg_estimate / grp_res$neg_se,
          `Pr(>|z|)` = grp_res$neg_pvalue,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    } else {
      # Overall WQS results
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        Term = "WQS Positive",
        Estimate = val_info$wqs_pos_estimate,
        `Std.Error` = val_info$wqs_pos_se,
        `z/t value` = val_info$wqs_pos_estimate / val_info$wqs_pos_se,
        `Pr(>|z|)` = val_info$wqs_pos_pvalue,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        Term = "WQS Negative",
        Estimate = val_info$wqs_neg_estimate,
        `Std.Error` = val_info$wqs_neg_se,
        `z/t value` = val_info$wqs_neg_estimate / val_info$wqs_neg_se,
        `Pr(>|z|)` = val_info$wqs_neg_pvalue,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  } else {
    # No validation: display Index Sum
    if (!is.null(object$groups)) {
      for (grp_name in names(object$groups)) {
        coef_rows[[length(coef_rows) + 1]] <- data.frame(
          Term = paste0(grp_name, " (positive)"),
          Estimate = object$pos_index_sum_by_group[[grp_name]],
          `Std.Error` = NA,
          `z/t value` = NA,
          `Pr(>|z|)` = NA,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        coef_rows[[length(coef_rows) + 1]] <- data.frame(
          Term = paste0(grp_name, " (negative)"),
          Estimate = object$neg_index_sum_by_group[[grp_name]],
          `Std.Error` = NA,
          `z/t value` = NA,
          `Pr(>|z|)` = NA,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    } else {
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        Term = "WQS Positive (index sum)",
        Estimate = object$pos_index_sum,
        `Std.Error` = NA,
        `z/t value` = NA,
        `Pr(>|z|)` = NA,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        Term = "WQS Negative (index sum)",
        Estimate = object$neg_index_sum,
        `Std.Error` = NA,
        `z/t value` = NA,
        `Pr(>|z|)` = NA,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  
  # Covariates
  if (!is.null(object$cov_coef) && !is.null(refit_coef_table)) {
    for (cov_nm in names(object$cov_coef)) {
      if (cov_nm %in% rownames(refit_coef_table)) {
        coef_rows[[length(coef_rows) + 1]] <- data.frame(
          Term = cov_nm,
          Estimate = refit_coef_table[cov_nm, "Estimate"],
          `Std.Error` = refit_coef_table[cov_nm, "Std. Error"],
          `z/t value` = refit_coef_table[cov_nm, "Estimate"] / refit_coef_table[cov_nm, "Std. Error"],
          `Pr(>|z|)` = refit_coef_table[cov_nm, refit_p_col],
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
      }
    }
  } else if (!is.null(object$cov_coef)) {
    for (cov_nm in names(object$cov_coef)) {
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        Term = cov_nm,
        Estimate = object$cov_coef[cov_nm],
        `Std.Error` = NA,
        `z/t value` = NA,
        `Pr(>|z|)` = NA,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  
  # Combine tables
  coef_df <- do.call(rbind, coef_rows)
  rownames(coef_df) <- NULL
  
  # Add significance markers
  coef_df$Signif <- ""
  coef_df$Signif[!is.na(coef_df$`Pr(>|z|)`) & coef_df$`Pr(>|z|)` < 0.001] <- "***"
  coef_df$Signif[!is.na(coef_df$`Pr(>|z|)`) & coef_df$`Pr(>|z|)` >= 0.001 & coef_df$`Pr(>|z|)` < 0.01] <- "**"
  coef_df$Signif[!is.na(coef_df$`Pr(>|z|)`) & coef_df$`Pr(>|z|)` >= 0.01 & coef_df$`Pr(>|z|)` < 0.05] <- "*"
  coef_df$Signif[!is.na(coef_df$`Pr(>|z|)`) & coef_df$`Pr(>|z|)` >= 0.05 & coef_df$`Pr(>|z|)` < 0.1] <- "."
  
  # Format and output
  print_df <- coef_df
  print_df$Estimate <- sprintf("%.6f", print_df$Estimate)
  print_df$`Std.Error` <- ifelse(is.na(coef_df$`Std.Error`), "---", sprintf("%.6f", coef_df$`Std.Error`))
  print_df$`z/t value` <- ifelse(is.na(coef_df$`z/t value`), "---", sprintf("%.3f", coef_df$`z/t value`))
  print_df$`Pr(>|z|)` <- ifelse(is.na(coef_df$`Pr(>|z|)`), "---", format.pval(coef_df$`Pr(>|z|)`, digits = 3))
  
  # Adjust column widths for output
  col_widths <- c(
    max(nchar(print_df$Term), 20),
    12, 12, 10, 12, 5
  )
  
  # Header
  header <- sprintf("%-*s %*s %*s %*s %*s %s",
                    col_widths[1], "",
                    col_widths[2], "Estimate",
                    col_widths[3], "Std.Error",
                    col_widths[4], "z/t value",
                    col_widths[5], "Pr(>|z|)", "")
  cat(header, "\n")
  
  # Each row
  for (i in seq_len(nrow(print_df))) {
    row_str <- sprintf("%-*s %*s %*s %*s %*s %s",
                       col_widths[1], print_df$Term[i],
                       col_widths[2], print_df$Estimate[i],
                       col_widths[3], print_df$`Std.Error`[i],
                       col_widths[4], print_df$`z/t value`[i],
                       col_widths[5], print_df$`Pr(>|z|)`[i],
                       print_df$Signif[i])
    cat(row_str, "\n")
  }
  
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  # ----- Additional information -----
  cat("\n---\n")
  if (!is.null(refit_info) && identical(attr(refit_info, "source"), "validation_info")) {
    cat("Validation: training n =", refit_info$n_train,
        ", validation n =", refit_info$n_val, "\n")
  } else if (!is.null(refit_info)) {
    cat("Refit:", object$refit, "(engine =", refit_info$engine, ")\n")
  }
  if (object$bootstrap) {
    cat("Bootstrap:", object$boot_info$n_successful, "iterations\n")
  }
  if (!is.null(object$groups)) {
    cat("Groups:", paste(names(object$groups), collapse = ", "), "\n")
  }
  
  invisible(object)
}


#' @keywords internal
.get_active_refit_info <- function(object) {
  if (!is.null(object$validation_info)) {
    info <- object$validation_info
    attr(info, "source") <- "validation_info"
    return(info)
  }
  info <- object$refit_info
  if (!is.null(info)) {
    attr(info, "source") <- "refit_info"
  }
  info
}

#' @keywords internal
.get_active_inference_source <- function(object) {
  info <- .get_active_refit_info(object)
  attr(info, "source")
}
