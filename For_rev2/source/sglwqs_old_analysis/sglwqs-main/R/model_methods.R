#' Predict Method for sglwqs Objects
#'
#' Generates WQS indices or model-based predictions for new data using
#' the quantile breaks learned during training.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param newdata Optional data frame or matrix containing exposure variables.
#'   If \code{NULL}, training quantile-transformed exposures are used.
#' @param covariates Optional covariate data frame or matrix for \code{newdata}.
#'   If omitted, covariates are searched in non-exposure columns of \code{newdata}.
#' @param type Prediction type: \code{"response"} (default), \code{"link"},
#'   \code{"index"}, or \code{"class"}.
#' @param threshold Classification threshold used when \code{type = "class"}
#'   for binomial models (default: 0.5).
#' @param se.fit Logical. If \code{TRUE}, return standard errors for predictions
#'   when uncertainty information is available (default: FALSE).
#'   Note: when the source is bootstrap matrices, the uncertainty reflects
#'   variation in exposure coefficients only (intercept and covariate
#'   coefficients are held fixed). For full joint uncertainty, use
#'   validation-based inference.
#' @param interval Character. \code{"none"} (default), \code{"confidence"},
#'   or \code{"prediction"}.
#' @param level Confidence level used for intervals (default: 0.95).
#' @param include_group_index Logical. If \code{TRUE} and groups are available,
#'   returns group-specific WQS indices when \code{type = "index"}.
#' @param ... Unused.
#'
#' @return
#' \itemize{
#'   \item \code{type = "index"}: data frame with \code{wqs_pos}, \code{wqs_neg}
#'   and optionally group-specific indices.
#'   \item \code{type = "link"}: numeric vector of linear predictors.
#'   \item \code{type = "response"}: numeric vector of predicted means/probabilities.
#'   \item \code{type = "class"}: integer vector (0/1) for binomial models.
#'   \item If \code{se.fit = TRUE}: list with \code{$fit}, \code{$se.fit},
#'   and \code{$source}.
#'   \item If \code{interval != "none"}: data frame with \code{fit}, \code{lwr},
#'   \code{upr}, and \code{source}.
#' }
#'
#' @examples
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#'
#' fit <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   covariates = sglwqs_example[c("age", "sex")],
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#'
#' # WQS indices for external data
#' idx <- predict(fit, newdata = sglwqs_example[1:5, ], type = "index")
#'
#' # Continuous predictions
#' pred <- predict(fit, newdata = sglwqs_example[1:5, ], type = "response")
#'
#' # Prediction uncertainty (requires bootstrap or validation)
#' fit_boot <- sglwqs(
#'   X = sglwqs_example[1:200, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:200],
#'   bootstrap = TRUE,
#'   n_boot = 5,
#'   keep_boot_matrices = TRUE,
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' pred_se <- predict(
#'   fit_boot,
#'   newdata = sglwqs_example[1:5, ],
#'   type = "response",
#'   se.fit = TRUE
#' )
#' head(idx)
#' pred
#'
#' @export
predict.sglwqs <- function(object, newdata = NULL, covariates = NULL,
                           type = c("response", "link", "index", "class"),
                           threshold = 0.5,
                           se.fit = FALSE,
                           interval = c("none", "confidence", "prediction"),
                           level = 0.95,
                           include_group_index = FALSE,
                           ...) {
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'.")
  }
  
  type <- match.arg(type)
  interval <- match.arg(interval)
  
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single numeric value in (0, 1).")
  }
  
  if (is.null(newdata)) {
    if (is.null(object$X_quantile)) {
      stop("Training quantile data are not stored in the fitted object.")
    }
    X_quantile <- object$X_quantile
    cov_matrix <- object$cov_matrix
  } else {
    X_raw <- .extract_exposure_matrix(newdata, object$var_names)
    qt_result <- quantile_transform(
      data = X_raw,
      n_quantiles = object$n_quantiles,
      var_names = object$var_names,
      breaks_list = object$quantile_breaks
    )
    X_quantile <- qt_result$q
    cov_matrix <- .extract_covariate_matrix(newdata, covariates, object)
  }
  
  index_df <- create_wqs_index(X_quantile, object$pos_weights, object$neg_weights)
  
  if (type == "index") {
    if (se.fit || interval != "none") {
      warning("Uncertainty is not defined for `type = \"index\"`. Returning indices only.")
    }
    
    if (isTRUE(include_group_index) && !is.null(object$groups)) {
      group_df <- .create_group_wqs_index(
        X_quantile = X_quantile,
        pos_weights = object$pos_weights,
        neg_weights = object$neg_weights,
        groups = object$groups
      )
      return(cbind(index_df, group_df))
    }
    return(index_df)
  }
  
  refit_info <- .get_active_refit_info(object)
  if (!is.null(refit_info) && !is.null(refit_info$refit_fit)) {
    pred_data <- .build_validation_prediction_data(
      object = object,
      X_quantile = X_quantile,
      cov_matrix = cov_matrix
    )
    
    if (type == "class") {
      response_pred <- as.numeric(stats::predict(
        refit_info$refit_fit,
        newdata = pred_data,
        type = "response"
      ))
      if (se.fit || interval != "none") {
        warning("Uncertainty output is not available for `type = \"class\"`. Returning classes only.")
      }
      return(as.integer(response_pred >= threshold))
    }
    
    base_pred <- as.numeric(stats::predict(
      refit_info$refit_fit,
      newdata = pred_data,
      type = type
    ))
    
    fit_link_default <- if (type == "link") {
      base_pred
    } else {
      as.numeric(stats::predict(refit_info$refit_fit, newdata = pred_data, type = "link"))
    }
    fit_response_default <- if (type == "response") {
      base_pred
    } else if (object$family == "binomial") {
      as.numeric(stats::predict(refit_info$refit_fit, newdata = pred_data, type = "response"))
    } else {
      fit_link_default
    }
    
    if (!se.fit && interval == "none") {
      return(base_pred)
    }
    
    unc <- .compute_prediction_uncertainty(
      object = object,
      X_quantile = X_quantile,
      cov_matrix = cov_matrix,
      eta = fit_link_default,
      response = fit_response_default,
      newdata = newdata,
      level = level
    )
    
    if (is.null(unc)) {
      stop("Prediction uncertainty is unavailable from the active refit model.")
    }
    
    fit_link <- if (!is.null(unc$fit_link)) unc$fit_link else fit_link_default
    fit_response <- if (!is.null(unc$fit_response)) unc$fit_response else fit_response_default
    fit_out <- if (type == "link") fit_link else fit_response
    se_out <- if (type == "link") unc$se_link else unc$se_response
    
    if (interval == "none") {
      return(list(
        fit = as.numeric(fit_out),
        se.fit = as.numeric(se_out),
        level = level,
        source = unc$source
      ))
    }
    
    interval_eff <- interval
    if (interval_eff == "prediction" && object$family != "gaussian") {
      warning("`interval = \"prediction\"` is only available for gaussian models. Using confidence interval.")
      interval_eff <- "confidence"
    }
    
    crit <- .prediction_interval_multiplier(refit_info$refit_fit, level)
    
    if (interval_eff == "confidence") {
      if (type == "link" && !is.null(unc$ci_link)) {
        lwr <- unc$ci_link[, 1]
        upr <- unc$ci_link[, 2]
      } else if (type != "link" && !is.null(unc$ci_response)) {
        lwr <- unc$ci_response[, 1]
        upr <- unc$ci_response[, 2]
      } else {
        lwr <- fit_out - crit * se_out
        upr <- fit_out + crit * se_out
      }
    } else {
      sigma <- unc$residual_sd
      if (!is.finite(sigma)) {
        warning("Residual scale is unavailable; using confidence interval instead.")
        lwr <- fit_out - crit * se_out
        upr <- fit_out + crit * se_out
        interval_eff <- "confidence"
      } else {
        total_se <- sqrt(se_out^2 + sigma^2)
        lwr <- fit_out - crit * total_se
        upr <- fit_out + crit * total_se
      }
    }
    
    return(data.frame(
      fit = as.numeric(fit_out),
      lwr = as.numeric(lwr),
      upr = as.numeric(upr),
      source = unc$source,
      interval = interval_eff,
      stringsAsFactors = FALSE
    ))
  }
  
  eta <- as.numeric(object$intercept)[1] +
    as.numeric(X_quantile %*% as.numeric(object$pos_coef)) -
    as.numeric(X_quantile %*% as.numeric(object$neg_coef))
  
  if (!is.null(object$cov_coef)) {
    if (is.null(cov_matrix)) {
      stop(
        "This model includes covariates. Supply covariates via `newdata` ",
        "or the `covariates` argument."
      )
    }
    eta <- eta + as.numeric(cov_matrix %*% as.numeric(object$cov_coef))
  }
  
  if (type == "link") {
    base_pred <- as.numeric(eta)
  } else if (object$family == "binomial") {
    prob <- stats::plogis(eta)
    if (type == "class") {
      if (se.fit || interval != "none") {
        warning("Uncertainty output is not available for `type = \"class\"`. Returning classes only.")
      }
      return(as.integer(prob >= threshold))
    }
    base_pred <- as.numeric(prob)
  } else {
    if (type == "class") {
      warning("type = 'class' is only available for binomial models. Returning response scale.")
    }
    base_pred <- as.numeric(eta)
  }
  
  if (!se.fit && interval == "none") {
    return(base_pred)
  }
  
  unc <- .compute_prediction_uncertainty(
    object = object,
    X_quantile = X_quantile,
    cov_matrix = cov_matrix,
    eta = as.numeric(eta),
    response = if (object$family == "binomial") as.numeric(stats::plogis(eta)) else as.numeric(eta),
    newdata = newdata,
    level = level
  )
  
  if (is.null(unc)) {
    stop(
      "Uncertainty information is unavailable. Fit with bootstrap = TRUE ",
      "(preferably keep_boot_matrices = TRUE) or validation = TRUE."
    )
  }
  
  fit_link <- if (!is.null(unc$fit_link)) unc$fit_link else as.numeric(eta)
  fit_response <- if (!is.null(unc$fit_response)) unc$fit_response else if (object$family == "binomial") {
    as.numeric(stats::plogis(fit_link))
  } else {
    fit_link
  }
  
  fit_out <- if (type == "link") fit_link else fit_response
  se_out <- if (type == "link") unc$se_link else unc$se_response
  
  if (interval == "none") {
    return(list(
      fit = as.numeric(fit_out),
      se.fit = as.numeric(se_out),
      level = level,
      source = unc$source
    ))
  }
  
  interval_eff <- interval
  if (interval_eff == "prediction" && object$family != "gaussian") {
    warning("`interval = \"prediction\"` is only available for gaussian models. Using confidence interval.")
    interval_eff <- "confidence"
  }
  
  alpha <- 1 - level
  # Use t distribution for survey refit (important for few-PSU designs)
  refit_info_ci <- .get_active_refit_info(object)
  glm_fit_ci <- if (!is.null(refit_info_ci)) {
    if (!is.null(refit_info_ci$refit_fit)) refit_info_ci$refit_fit else refit_info_ci$glm_fit
  } else {
    NULL
  }
  z <- .prediction_interval_multiplier(glm_fit_ci, level)

  if (interval_eff == "confidence") {
    if (type == "link" && !is.null(unc$ci_link)) {
      lwr <- unc$ci_link[, 1]
      upr <- unc$ci_link[, 2]
    } else if (type != "link" && !is.null(unc$ci_response)) {
      lwr <- unc$ci_response[, 1]
      upr <- unc$ci_response[, 2]
    } else {
      lwr <- fit_out - z * se_out
      upr <- fit_out + z * se_out
    }
  } else {
    sigma <- unc$residual_sd
    if (!is.finite(sigma)) {
      warning("Residual scale is unavailable; using confidence interval instead.")
      lwr <- fit_out - z * se_out
      upr <- fit_out + z * se_out
      interval_eff <- "confidence"
    } else {
      total_se <- sqrt(se_out^2 + sigma^2)
      lwr <- fit_out - z * total_se
      upr <- fit_out + z * total_se
    }
  }
  
  data.frame(
    fit = as.numeric(fit_out),
    lwr = as.numeric(lwr),
    upr = as.numeric(upr),
    source = unc$source,
    interval = interval_eff,
    stringsAsFactors = FALSE
  )
}


#' Prediction Intervals for sglwqs Objects
#'
#' Convenience wrapper for interval predictions. When the uncertainty source
#' is bootstrap matrices, intervals reflect exposure coefficient variation only
#' (approximate); validation-based intervals account for full model uncertainty.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param newdata Optional predictor data.
#' @param covariates Optional covariate data.
#' @param type \code{"response"} (default) or \code{"link"}.
#' @param level Confidence level (default: 0.95).
#' @param interval \code{"confidence"} (default) or \code{"prediction"}.
#' @param ... Additional arguments passed to \code{predict.sglwqs()}.
#'
#' @return A data frame with \code{fit}, \code{lwr}, \code{upr}, and metadata.
#'
#' @examples
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#'
#' fit <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   bootstrap = TRUE,
#'   n_boot = 20,
#'   keep_boot_matrices = TRUE,
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#'
#' predict_interval(fit, newdata = sglwqs_example[1:5, ], level = 0.95)
#'
#' @export
predict_interval <- function(object, ...) {
  UseMethod("predict_interval")
}


#' @rdname predict_interval
#' @export
predict_interval.sglwqs <- function(object, newdata = NULL, covariates = NULL,
                                    type = c("response", "link"),
                                    level = 0.95,
                                    interval = c("confidence", "prediction"),
                                    ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  
  predict.sglwqs(
    object = object,
    newdata = newdata,
    covariates = covariates,
    type = type,
    se.fit = FALSE,
    interval = interval,
    level = level,
    ...
  )
}


#' Coefficients for sglwqs Objects
#'
#' Extracts directional coefficients or weights as named vectors.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param what Character. One of \code{"coefficients"}, \code{"weights"},
#'   or \code{"all"}.
#' @param direction Character. \code{"both"} (default), \code{"positive"},
#'   or \code{"negative"}.
#' @param include_intercept Logical. Include intercept when \code{what}
#'   includes coefficients (default: TRUE).
#' @param include_covariates Logical. Include covariate coefficients when
#'   \code{what} includes coefficients (default: TRUE).
#' @param ... Unused.
#'
#' @return A named numeric vector, or a list with \code{$coefficients} and
#'   \code{$weights} when \code{what = "all"}.
#'
#' @examples
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#'
#' fit <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#'
#' coef(fit, what = "weights", direction = "positive")
#' coef(fit, what = "coefficients", direction = "both")
#'
#' @export
coef.sglwqs <- function(object,
                        what = c("coefficients", "weights", "all"),
                        direction = c("both", "positive", "negative"),
                        include_intercept = TRUE,
                        include_covariates = TRUE,
                        ...) {
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'.")
  }
  
  what <- match.arg(what)
  direction <- match.arg(direction)
  
  if (what == "all") {
    return(list(
      coefficients = coef.sglwqs(
        object,
        what = "coefficients",
        direction = direction,
        include_intercept = include_intercept,
        include_covariates = include_covariates
      ),
      weights = coef.sglwqs(
        object,
        what = "weights",
        direction = direction,
        include_intercept = FALSE,
        include_covariates = FALSE
      )
    ))
  }
  
  if (what == "weights") {
    return(.combine_directional_vectors(
      pos = object$pos_weights,
      neg = object$neg_weights,
      direction = direction,
      pos_prefix = "pos_w.",
      neg_prefix = "neg_w."
    ))
  }
  
  coef_vec <- .combine_directional_vectors(
    pos = object$pos_coef,
    neg = object$neg_coef,
    direction = direction,
    pos_prefix = "pos_b.",
    neg_prefix = "neg_b."
  )
  
  if (isTRUE(include_intercept)) {
    coef_vec <- c("(Intercept)" = as.numeric(object$intercept)[1], coef_vec)
  }
  
  if (isTRUE(include_covariates) && !is.null(object$cov_coef)) {
    coef_vec <- c(coef_vec, object$cov_coef)
  }
  
  coef_vec
}


#' Confidence Intervals for sglwqs Objects
#'
#' Returns confidence intervals for validation estimates or bootstrap coefficients.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param parm Optional term names (or numeric indices) to extract.
#' @param level Confidence level (default: 0.95).
#' @param type Character. \code{"validation"} or \code{"bootstrap"}.
#'   If \code{NULL}, the method prefers \code{"validation"} when available,
#'   otherwise \code{"bootstrap"}.
#' @param direction Character. \code{"both"} (default), \code{"positive"},
#'   or \code{"negative"}; used for \code{type = "bootstrap"}.
#' @param ... Unused.
#'
#' @return A two-column matrix with lower and upper confidence limits.
#'
#' @examples
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#'
#' fit_boot <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   bootstrap = TRUE,
#'   n_boot = 20,
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#'
#' confint(fit_boot, type = "bootstrap", direction = "positive")
#'
#' @export
confint.sglwqs <- function(object, parm = NULL, level = 0.95, type = NULL,
                           direction = c("both", "positive", "negative"), ...) {
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'.")
  }
  
  direction <- match.arg(direction)
  
  if (is.null(type)) {
    if (!is.null(.get_active_refit_info(object))) {
      type <- "validation"
    } else if (isTRUE(object$bootstrap)) {
      type <- "bootstrap"
    } else {
      stop("No confidence interval source found. Fit with validation = TRUE or bootstrap = TRUE.")
    }
  }
  
  type <- match.arg(type, c("validation", "bootstrap"))
  
  if (type == "validation") {
    if (is.null(.get_active_refit_info(object))) {
      stop("No inference results available. Fit with `validation = TRUE` or `refit = \"full\"` to enable inference.")
    }
    val <- summary_validation(object, conf_level = level)
    ci <- as.matrix(val[, c("ci_lower", "ci_upper"), drop = FALSE])
    rownames(ci) <- val$term
  } else {
    if (!isTRUE(object$bootstrap)) {
      stop("Bootstrap confidence intervals require `bootstrap = TRUE`.")
    }
    boot <- summary_bootstrap(object, conf_level = level, min_freq = 0)
    if (direction != "both") {
      boot <- boot[boot$direction == direction, , drop = FALSE]
    }
    if (nrow(boot) == 0) {
      stop("No rows available after filtering by direction.")
    }
    ci <- as.matrix(boot[, c("ci_lower", "ci_upper"), drop = FALSE])
    rownames(ci) <- paste0(boot$variable, " [", boot$direction, "]")
  }
  
  colnames(ci) <- .confint_colnames(level)
  .subset_confint_matrix(ci, parm)
}


#' Cross-validation Curve for Lambda Selection
#'
#' Visualizes the \code{cv.sparsegl} cross-validation curve stored in a
#' fitted \code{sglwqs} object.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param show_lambda Character. Which reference lines to draw:
#'   \code{"both"} (default), \code{"min"}, \code{"1se"}, or \code{"none"}.
#' @param show_se Logical. Show \eqn{\pm 1 SE} ribbon (default: TRUE).
#' @param point_alpha Alpha for CV points (default: 0.6).
#' @param base_size Base font size (default: 11).
#' @param title Plot title.
#' @param ... Unused.
#'
#' @return A ggplot object.
#'
#' @examples
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#'
#' fit <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   groups = list(
#'     Metals = paste0("metal", 1:4),
#'     Pesticides = paste0("pesticide", 1:4)
#'   ),
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#'
#' plot_cv(fit)
#'
#' @export
plot_cv <- function(object, ...) {
  UseMethod("plot_cv")
}


#' @rdname plot_cv
#' @export
plot_cv.sglwqs <- function(object,
                           show_lambda = c("both", "min", "1se", "none"),
                           show_se = TRUE,
                           point_alpha = 0.6,
                           base_size = 11,
                           title = "Cross-validation Curve",
                           ...) {
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'.")
  }
  
  show_lambda <- match.arg(show_lambda)
  
  cv_fit <- object$fit
  if (is.null(cv_fit$lambda) || is.null(cv_fit$cvm)) {
    stop("No cross-validation information found in `object$fit`.")
  }
  
  cv_df <- data.frame(
    lambda = cv_fit$lambda,
    log_lambda = log(cv_fit$lambda),
    cvm = cv_fit$cvm,
    cvsd = if (is.null(cv_fit$cvsd)) NA_real_ else cv_fit$cvsd
  )
  
  p <- ggplot2::ggplot(cv_df, ggplot2::aes(x = .data$log_lambda, y = .data$cvm)) +
    ggplot2::geom_line(color = "#2C7FB8", linewidth = 0.8) +
    ggplot2::geom_point(color = "#2C7FB8", size = 1.2, alpha = point_alpha)
  
  if (show_se && all(is.finite(cv_df$cvsd))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$cvm - .data$cvsd, ymax = .data$cvm + .data$cvsd),
      fill = "#A6CEE3",
      alpha = 0.35
    )
  }
  
  lambda_marks <- data.frame(
    label = c("lambda.min", "lambda.1se"),
    lambda = c(cv_fit$lambda.min, cv_fit$lambda.1se),
    stringsAsFactors = FALSE
  )
  lambda_marks <- lambda_marks[is.finite(lambda_marks$lambda), , drop = FALSE]
  
  if (show_lambda != "none" && nrow(lambda_marks) > 0) {
    if (show_lambda == "min") {
      lambda_marks <- lambda_marks[lambda_marks$label == "lambda.min", , drop = FALSE]
    } else if (show_lambda == "1se") {
      lambda_marks <- lambda_marks[lambda_marks$label == "lambda.1se", , drop = FALSE]
    }
    
    p <- p +
      ggplot2::geom_vline(
        data = lambda_marks,
        ggplot2::aes(xintercept = log(.data$lambda), color = .data$label),
        linetype = "dashed",
        linewidth = 0.7,
        show.legend = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("lambda.min" = "#D7301F", "lambda.1se" = "#238B45"),
        name = "Selected lambda"
      )
  }
  
  n_pos <- sum(object$pos_weights > 0)
  n_neg <- sum(object$neg_weights > 0)
  group_text <- if (!is.null(object$groups)) {
    paste0(" | groups: ", length(object$groups))
  } else {
    ""
  }
  
  subtitle <- paste0(
    "Selected exposures (+/-): ", n_pos, "/", object$p, " , ", n_neg, "/", object$p,
    group_text
  )
  
  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "log(lambda)",
      y = "Cross-validation error"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}


#' @rdname plot_cv
#' @param show_imputations For \code{sglwqs_mids} objects, draw thin
#'   imputation-specific CV curves in the background (default: TRUE).
#' @export
plot_cv.sglwqs_mids <- function(object,
                                show_lambda = c("both", "min", "1se", "none"),
                                show_se = TRUE,
                                point_alpha = 0.6,
                                base_size = 11,
                                title = "Pooled Cross-validation Curve",
                                show_imputations = TRUE,
                                ...) {
  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'.")
  }

  show_lambda <- match.arg(show_lambda)
  cv_info <- .pool_mi_cv_curve(object)

  p <- ggplot2::ggplot(cv_info$pooled_df, ggplot2::aes(x = .data$log_lambda, y = .data$cvm))

  if (show_imputations && !is.null(cv_info$individual_df) && nrow(cv_info$individual_df) > 0) {
    p <- p +
      ggplot2::geom_line(
        data = cv_info$individual_df,
        ggplot2::aes(x = .data$log_lambda, y = .data$cvm, group = .data$imputation),
        inherit.aes = FALSE,
        color = "gray80",
        alpha = 0.7,
        linewidth = 0.4
      )
  }

  if (show_se && all(is.finite(cv_info$pooled_df$cvsd))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$cvm - .data$cvsd, ymax = .data$cvm + .data$cvsd),
      fill = "#A6CEE3",
      alpha = 0.35
    )
  }

  p <- p +
    ggplot2::geom_line(color = "#2C7FB8", linewidth = 0.8) +
    ggplot2::geom_point(color = "#2C7FB8", size = 1.2, alpha = point_alpha)

  lambda_marks <- data.frame(
    label = c("lambda.min", "lambda.1se"),
    lambda = c(cv_info$lambda.min, cv_info$lambda.1se),
    stringsAsFactors = FALSE
  )
  lambda_marks <- lambda_marks[is.finite(lambda_marks$lambda), , drop = FALSE]

  if (show_lambda != "none" && nrow(lambda_marks) > 0) {
    if (show_lambda == "min") {
      lambda_marks <- lambda_marks[lambda_marks$label == "lambda.min", , drop = FALSE]
    } else if (show_lambda == "1se") {
      lambda_marks <- lambda_marks[lambda_marks$label == "lambda.1se", , drop = FALSE]
    }

    p <- p +
      ggplot2::geom_vline(
        data = lambda_marks,
        ggplot2::aes(xintercept = log(.data$lambda), color = .data$label),
        linetype = "dashed",
        linewidth = 0.7,
        show.legend = TRUE
      ) +
      ggplot2::scale_color_manual(
        values = c("lambda.min" = "#D7301F", "lambda.1se" = "#238B45"),
        name = "Selected lambda"
      )
  }

  n_pos <- if (!is.null(object$pooled$pos_weights)) sum(object$pooled$pos_weights > 0) else NA_integer_
  n_neg <- if (!is.null(object$pooled$neg_weights)) sum(object$pooled$neg_weights > 0) else NA_integer_
  p_total <- if (!is.null(object$pooled$var_names)) {
    length(object$pooled$var_names)
  } else {
    tmp_total <- suppressWarnings(max(c(n_pos, n_neg), na.rm = TRUE))
    if (is.finite(tmp_total)) tmp_total else NA_integer_
  }
  group_text <- if (!is.null(object$groups)) {
    paste0(" | groups: ", length(object$groups))
  } else {
    ""
  }

  subtitle <- paste0(
    "Selected exposures (+/-): ", n_pos, "/", p_total, " , ", n_neg, "/", p_total,
    group_text
  )

  caption_text <- "Pointwise mean CV error across imputations after interpolation to a common log(lambda) grid."
  if (show_imputations && !is.null(cv_info$individual_df) && nrow(cv_info$individual_df) > 0) {
    caption_text <- paste0(caption_text, " Gray lines denote imputation-specific curves.")
  }

  p +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption_text,
      x = "log(lambda)",
      y = "Cross-validation error"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(hjust = 0, size = base_size * 0.8,
                                           color = "gray40", face = "italic")
    )
}


#' @keywords internal
.pool_mi_cv_curve <- function(object) {
  fits <- Filter(function(f) {
    inherits(f, "sglwqs") &&
      !is.null(f$fit) &&
      !is.null(f$fit$lambda) &&
      !is.null(f$fit$cvm)
  }, object$fits)

  if (length(fits) == 0) {
    stop("No cross-validation information found in successful imputations.")
  }

  extract_curve <- function(fit, id) {
    cv_fit <- fit$fit
    df <- data.frame(
      imputation = id,
      lambda = as.numeric(cv_fit$lambda),
      log_lambda = log(as.numeric(cv_fit$lambda)),
      cvm = as.numeric(cv_fit$cvm),
      cvsd = if (is.null(cv_fit$cvsd)) NA_real_ else as.numeric(cv_fit$cvsd),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$lambda) & is.finite(df$log_lambda) & is.finite(df$cvm), , drop = FALSE]
    df <- df[order(df$log_lambda), , drop = FALSE]
    df
  }

  cv_list <- Filter(function(df) nrow(df) > 1, Map(extract_curve, fits, seq_along(fits)))
  if (length(cv_list) == 0) {
    stop("No valid cross-validation curves found in successful imputations.")
  }

  n_grid <- max(vapply(cv_list, nrow, integer(1)))
  mins <- vapply(cv_list, function(df) min(df$log_lambda), numeric(1))
  maxs <- vapply(cv_list, function(df) max(df$log_lambda), numeric(1))

  grid_low <- max(mins)
  grid_high <- min(maxs)
  if (!is.finite(grid_low) || !is.finite(grid_high) || grid_low >= grid_high) {
    grid_low <- min(mins)
    grid_high <- max(maxs)
  }

  common_grid <- seq(grid_low, grid_high, length.out = n_grid)

  interp_mat <- vapply(cv_list, function(df) {
    stats::approx(x = df$log_lambda, y = df$cvm, xout = common_grid, rule = 2, ties = mean)$y
  }, numeric(length(common_grid)))

  se_mat <- vapply(cv_list, function(df) {
    if (!all(is.finite(df$cvsd))) {
      rep(NA_real_, length(common_grid))
    } else {
      stats::approx(
        x = df$log_lambda,
        y = df$cvsd,
        xout = common_grid,
        rule = 2,
        ties = mean
      )$y
    }
  }, numeric(length(common_grid)))

  pooled_cvm <- rowMeans(interp_mat, na.rm = TRUE)
  within_var <- if (all(!is.finite(se_mat))) {
    rep(0, length(common_grid))
  } else {
    rowMeans(se_mat^2, na.rm = TRUE)
  }
  within_var[!is.finite(within_var)] <- 0
  between_var <- apply(interp_mat, 1, function(x) {
    if (sum(is.finite(x)) > 1) stats::var(x, na.rm = TRUE) else 0
  })
  pooled_cvsd <- sqrt(pmax(within_var + between_var, 0))

  idx_min <- which.min(pooled_cvm)
  threshold <- pooled_cvm[idx_min] + pooled_cvsd[idx_min]
  eligible <- which(pooled_cvm <= threshold)
  idx_1se <- if (length(eligible) > 0) max(eligible) else idx_min

  pooled_df <- data.frame(
    lambda = exp(common_grid),
    log_lambda = common_grid,
    cvm = pooled_cvm,
    cvsd = pooled_cvsd,
    stringsAsFactors = FALSE
  )

  individual_df <- do.call(rbind, cv_list)
  rownames(individual_df) <- NULL

  list(
    pooled_df = pooled_df,
    individual_df = individual_df,
    lambda.min = exp(common_grid[idx_min]),
    lambda.1se = exp(common_grid[idx_1se])
  )
}




# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @keywords internal
.extract_exposure_matrix <- function(newdata, var_names) {
  if (is.data.frame(newdata)) {
    if (all(var_names %in% names(newdata))) {
      X <- newdata[, var_names, drop = FALSE]
    } else if (ncol(newdata) == length(var_names)) {
      X <- newdata
      names(X) <- var_names
    } else {
      stop(
        "`newdata` must include all exposure columns used in training: ",
        paste(var_names, collapse = ", ")
      )
    }
    
    is_num <- vapply(X, is.numeric, logical(1))
    if (!all(is_num)) {
      bad <- names(X)[!is_num]
      stop("Exposure columns must be numeric. Non-numeric columns: ", paste(bad, collapse = ", "))
    }
    X <- as.matrix(X)
    colnames(X) <- var_names
    return(X)
  }
  
  if (is.matrix(newdata)) {
    if (!is.numeric(newdata)) {
      stop("`newdata` matrix must be numeric.")
    }
    
    if (!is.null(colnames(newdata)) && all(var_names %in% colnames(newdata))) {
      return(newdata[, var_names, drop = FALSE])
    }
    
    if (ncol(newdata) == length(var_names)) {
      X <- newdata
      colnames(X) <- var_names
      return(X)
    }
  }
  
  stop("`newdata` must be a data.frame or numeric matrix with exposure columns.")
}


#' @keywords internal
.extract_covariate_matrix <- function(newdata, covariates, object) {
  expected <- object$cov_names
  if (is.null(expected)) {
    return(NULL)
  }
  
  if (!is.null(covariates)) {
    cov_processed <- process_covariates(covariates)
    return(.align_covariate_matrix(cov_processed$matrix, expected))
  }
  
  if (is.data.frame(newdata)) {
    if (all(expected %in% names(newdata))) {
      cov_matrix <- as.matrix(newdata[, expected, drop = FALSE])
      if (!is.numeric(cov_matrix)) {
        stop("Covariate columns extracted from `newdata` must be numeric.")
      }
      return(cov_matrix)
    }
    
    cov_candidates <- setdiff(names(newdata), object$var_names)
    if (length(cov_candidates) == 0) {
      stop(
        "This model includes covariates. Provide covariate columns in `newdata` ",
        "or use the `covariates` argument."
      )
    }
    
    cov_processed <- process_covariates(newdata[, cov_candidates, drop = FALSE])
    return(.align_covariate_matrix(cov_processed$matrix, expected))
  }
  
  if (is.matrix(newdata) && !is.null(colnames(newdata)) && all(expected %in% colnames(newdata))) {
    cov_matrix <- newdata[, expected, drop = FALSE]
    if (!is.numeric(cov_matrix)) {
      stop("Covariate columns extracted from `newdata` must be numeric.")
    }
    return(cov_matrix)
  }
  
  stop(
    "This model includes covariates. Provide covariate columns in `newdata` ",
    "or use the `covariates` argument."
  )
}


#' @keywords internal
.align_covariate_matrix <- function(cov_matrix, expected_names) {
  if (is.null(cov_matrix)) {
    stop("Covariates are required, but none were provided.")
  }
  
  cov_matrix <- as.matrix(cov_matrix)
  if (is.null(colnames(cov_matrix))) {
    stop("Covariate matrix must have column names for alignment.")
  }
  
  if (!is.numeric(cov_matrix)) {
    stop("Covariate matrix must be numeric after preprocessing.")
  }
  
  aligned <- matrix(0, nrow = nrow(cov_matrix), ncol = length(expected_names))
  colnames(aligned) <- expected_names
  rownames(aligned) <- rownames(cov_matrix)
  
  common <- intersect(expected_names, colnames(cov_matrix))
  if (length(common) == 0) {
    stop(
      "No covariate columns in new data match expected columns: ",
      paste(expected_names, collapse = ", ")
    )
  }
  
  aligned[, common] <- cov_matrix[, common, drop = FALSE]
  
  missing <- setdiff(expected_names, colnames(cov_matrix))
  if (length(missing) > 0) {
    warning(
      "Missing covariate columns were filled with 0: ",
      paste(missing, collapse = ", ")
    )
  }
  
  aligned
}


#' @keywords internal
.create_group_wqs_index <- function(X_quantile, pos_weights, neg_weights, groups) {
  var_names <- colnames(X_quantile)
  out <- data.frame(matrix(nrow = nrow(X_quantile), ncol = 0))
  
  for (grp_name in names(groups)) {
    grp_idx <- which(var_names %in% groups[[grp_name]])
    if (length(grp_idx) == 0) next
    
    grp_pos <- as.numeric(pos_weights[grp_idx])
    grp_neg <- as.numeric(neg_weights[grp_idx])
    
    if (sum(grp_pos) > 0) grp_pos <- grp_pos / sum(grp_pos)
    if (sum(grp_neg) > 0) grp_neg <- grp_neg / sum(grp_neg)
    
    X_grp <- X_quantile[, grp_idx, drop = FALSE]
    out[[paste0(grp_name, "_pos")]] <- as.numeric(X_grp %*% grp_pos)
    out[[paste0(grp_name, "_neg")]] <- as.numeric(X_grp %*% grp_neg)
  }
  
  out
}


#' @keywords internal
.compute_prediction_uncertainty <- function(object, X_quantile, cov_matrix,
                                            eta, response, newdata = NULL, level = 0.95) {
  if (!is.null(.get_active_refit_info(object))) {
    unc <- .prediction_uncertainty_from_validation_glm(
      object = object,
      X_quantile = X_quantile,
      cov_matrix = cov_matrix,
      newdata = newdata,
      level = level
    )
    if (!is.null(unc)) return(unc)
  }
  
  unc <- .prediction_uncertainty_from_bootstrap_matrices(
    object = object,
    X_quantile = X_quantile,
    cov_matrix = cov_matrix,
    level = level
  )
  if (!is.null(unc)) return(unc)
  
  unc <- .prediction_uncertainty_from_bootstrap_se(
    object = object,
    X_quantile = X_quantile,
    eta = eta,
    response = response,
    level = level
  )
  if (!is.null(unc)) return(unc)
  
  .prediction_uncertainty_from_validation_glm(
    object = object,
    X_quantile = X_quantile,
    cov_matrix = cov_matrix,
    newdata = newdata,
    level = level
  )
}


#' Note: this is an approximate prediction uncertainty that varies exposure
#' coefficients across bootstrap draws but uses the single-fit intercept and
#' covariate coefficients. It does not capture full joint parameter uncertainty.
#' @keywords internal
.prediction_uncertainty_from_bootstrap_matrices <- function(object, X_quantile, cov_matrix, level) {
  if (!isTRUE(object$bootstrap) || is.null(object$boot_info)) {
    return(NULL)
  }
  
  boot_pos <- object$boot_info$boot_pos_coef
  boot_neg <- object$boot_info$boot_neg_coef
  if (is.null(boot_pos) || is.null(boot_neg)) {
    return(NULL)
  }
  
  boot_pos <- as.matrix(boot_pos)
  boot_neg <- as.matrix(boot_neg)
  
  # Prioritize boot_success if available (includes all-zero coefficients as valid results)
  if (!is.null(object$boot_info$boot_success)) {
    valid <- object$boot_info$boot_success
  } else {
    valid <- rowSums(abs(boot_pos) + abs(boot_neg)) > 0
  }
  if (!any(valid)) {
    return(NULL)
  }

  boot_pos <- boot_pos[valid, , drop = FALSE]
  boot_neg <- boot_neg[valid, , drop = FALSE]
  
  eta_draws <- X_quantile %*% t(boot_pos) - X_quantile %*% t(boot_neg)
  eta_draws <- sweep(eta_draws, 1, as.numeric(object$intercept)[1], FUN = "+")
  
  if (!is.null(object$cov_coef) && !is.null(cov_matrix)) {
    cov_part <- as.numeric(cov_matrix %*% as.numeric(object$cov_coef))
    eta_draws <- sweep(eta_draws, 1, cov_part, FUN = "+")
  }
  
  alpha <- 1 - level
  ci_eta <- t(apply(eta_draws, 1, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2)))
  se_eta <- apply(eta_draws, 1, stats::sd)
  fit_eta <- rowMeans(eta_draws)
  
  if (object$family == "binomial") {
    resp_draws <- stats::plogis(eta_draws)
    ci_resp <- t(apply(resp_draws, 1, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2)))
    se_resp <- apply(resp_draws, 1, stats::sd)
    fit_resp <- rowMeans(resp_draws)
  } else {
    ci_resp <- ci_eta
    se_resp <- se_eta
    fit_resp <- fit_eta
  }
  
  list(
    source = "bootstrap_matrices",
    fit_link = as.numeric(fit_eta),
    fit_response = as.numeric(fit_resp),
    se_link = as.numeric(se_eta),
    se_response = as.numeric(se_resp),
    ci_link = ci_eta,
    ci_response = ci_resp,
    residual_sd = NA_real_
  )
}


#' @keywords internal
.prediction_uncertainty_from_bootstrap_se <- function(object, X_quantile, eta, response, level) {
  if (!isTRUE(object$bootstrap) || is.null(object$boot_info)) {
    return(NULL)
  }

  se_pos <- object$boot_info$se_pos_coef
  se_neg <- object$boot_info$se_neg_coef
  if (is.null(se_pos) || is.null(se_neg)) {
    return(NULL)
  }

  se_pos <- as.numeric(se_pos)
  se_neg <- as.numeric(se_neg)

  if (length(se_pos) != ncol(X_quantile) || length(se_neg) != ncol(X_quantile)) {
    return(NULL)
  }

  # Note: This assumes independence between pos and neg coefficients.

  # Covariance is ignored because only marginal SEs are available in this path.
  # For accurate SE accounting for covariance, use keep_boot_matrices = TRUE
  # (which routes to .prediction_uncertainty_from_bootstrap_matrices instead).
  var_eta <- as.numeric((X_quantile^2) %*% (se_pos^2 + se_neg^2))
  se_eta <- sqrt(pmax(var_eta, 0))

  if (object$family == "binomial") {
    se_resp <- abs(response * (1 - response) * se_eta)
  } else {
    se_resp <- se_eta
  }


  # Use t distribution if survey refit exists
  refit_info_bs <- .get_active_refit_info(object)
  glm_fit_bs <- if (!is.null(refit_info_bs)) {
    if (!is.null(refit_info_bs$refit_fit)) refit_info_bs$refit_fit else refit_info_bs$glm_fit
  } else {
    NULL
  }
  z <- .prediction_interval_multiplier(glm_fit_bs, level)

  ci_eta <- cbind(eta - z * se_eta, eta + z * se_eta)
  if (object$family == "binomial") {
    ci_resp <- cbind(
      pmax(0, response - z * se_resp),
      pmin(1, response + z * se_resp)
    )
  } else {
    ci_resp <- ci_eta
  }
  
  list(
    source = "bootstrap_se",
    fit_link = as.numeric(eta),
    fit_response = as.numeric(response),
    se_link = as.numeric(se_eta),
    se_response = as.numeric(se_resp),
    ci_link = ci_eta,
    ci_response = ci_resp,
    residual_sd = NA_real_
  )
}


#' @keywords internal
.prediction_uncertainty_from_validation_glm <- function(object, X_quantile, cov_matrix, newdata, level) {
  refit_info <- .get_active_refit_info(object)
  if (is.null(refit_info)) {
    return(NULL)
  }
  
  glm_fit <- if (!is.null(refit_info$refit_fit)) refit_info$refit_fit else refit_info$glm_fit
  if (is.null(glm_fit)) {
    return(NULL)
  }
  pred_data <- .build_validation_prediction_data(
    object = object,
    X_quantile = X_quantile,
    cov_matrix = cov_matrix
  )
  
  pred_link_raw <- tryCatch(
    stats::predict(glm_fit, newdata = pred_data, type = "link", se.fit = TRUE),
    error = function(e) NULL
  )
  pred_link <- .extract_prediction_with_se(pred_link_raw)
  if (is.null(pred_link) || is.null(pred_link$fit) || is.null(pred_link$se.fit)) {
    return(NULL)
  }
  
  fit_link <- as.numeric(pred_link$fit)
  se_link <- as.numeric(pred_link$se.fit)
  
  if (object$family == "binomial") {
    pred_response_raw <- tryCatch(
      stats::predict(glm_fit, newdata = pred_data, type = "response", se.fit = TRUE),
      error = function(e) NULL
    )
    pred_response <- .extract_prediction_with_se(pred_response_raw)
    if (!is.null(pred_response)) {
      fit_response <- as.numeric(pred_response$fit)
      se_response <- as.numeric(pred_response$se.fit)
    } else {
      fit_response <- stats::plogis(fit_link)
      se_response <- abs(fit_response * (1 - fit_response) * se_link)
    }
  } else {
    fit_response <- fit_link
    se_response <- se_link
  }
  
  crit <- .prediction_interval_multiplier(glm_fit, level)
  
  ci_link <- cbind(fit_link - crit * se_link, fit_link + crit * se_link)
  if (object$family == "binomial") {
    ci_response <- cbind(
      pmax(0, fit_response - crit * se_response),
      pmin(1, fit_response + crit * se_response)
    )
  } else {
    ci_response <- ci_link
  }
  
  resid_sd <- if (object$family == "gaussian") {
    sqrt(tryCatch(summary(glm_fit)$dispersion, error = function(e) NA_real_))
  } else {
    NA_real_
  }
  source <- if (inherits(glm_fit, "svyglm")) "svyglm_refit" else "glm_refit"
  
  list(
    source = source,
    fit_link = fit_link,
    fit_response = fit_response,
    se_link = se_link,
    se_response = se_response,
    ci_link = ci_link,
    ci_response = ci_response,
    residual_sd = resid_sd
  )
}


# .prediction_interval_multiplier() is defined in utils.R
# (handles svyglm, gaussian GLM with estimated dispersion, and binomial/poisson)


#' Build prediction data frame for validation GLM
#'
#' Constructs the design matrix used by \code{predict.glm()} on the
#' validation refit, including group-level WQS indices when applicable
#' and covariate columns.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param cov_matrix Covariate matrix, or \code{NULL}.
#' @return A data frame suitable for \code{predict()}.
#' @keywords internal
.build_validation_prediction_data <- function(object, X_quantile, cov_matrix) {
  refit_info <- .get_active_refit_info(object)
  if (!is.null(refit_info) && isTRUE(refit_info$group_inference) && !is.null(object$groups)) {
    pred_data <- .create_group_wqs_index(
      X_quantile = X_quantile,
      pos_weights = object$pos_weights,
      neg_weights = object$neg_weights,
      groups = object$groups
    )
    names(pred_data) <- make.names(names(pred_data), unique = TRUE)
  } else {
    pred_data <- create_wqs_index(X_quantile, object$pos_weights, object$neg_weights)
  }
  
  if (!is.null(cov_matrix)) {
    cov_df <- as.data.frame(cov_matrix)
    pred_data <- cbind(pred_data, cov_df)
  }
  
  pred_data
}


#' @keywords internal
.extract_prediction_with_se <- function(pred_obj) {
  if (is.null(pred_obj)) {
    return(NULL)
  }
  if (is.list(pred_obj) && !is.null(pred_obj$fit) && !is.null(pred_obj$se.fit)) {
    return(list(
      fit = as.numeric(pred_obj$fit),
      se.fit = as.numeric(pred_obj$se.fit)
    ))
  }
  
  pred_var <- attr(pred_obj, "var")
  if (!is.null(pred_var)) {
    se <- if (is.matrix(pred_var)) {
      sqrt(pmax(diag(pred_var), 0))
    } else {
      sqrt(pmax(as.numeric(pred_var), 0))
    }
    return(list(
      fit = as.numeric(pred_obj),
      se.fit = as.numeric(se)
    ))
  }
  
  NULL
}


#' @keywords internal
.combine_directional_vectors <- function(pos, neg, direction,
                                         pos_prefix = "pos.",
                                         neg_prefix = "neg.") {
  pos_names <- names(pos)
  neg_names <- names(neg)
  pos <- as.numeric(pos)
  neg <- as.numeric(neg)
  names(pos) <- pos_names
  names(neg) <- neg_names
  
  if (direction == "positive") {
    return(pos)
  }
  if (direction == "negative") {
    return(neg)
  }
  
  out <- c(pos, neg)
  names(out) <- c(paste0(pos_prefix, names(pos)), paste0(neg_prefix, names(neg)))
  out
}


#' @keywords internal
.confint_colnames <- function(level) {
  alpha <- (1 - level) / 2
  lower <- sprintf("%.1f%%", 100 * alpha)
  upper <- sprintf("%.1f%%", 100 * (1 - alpha))
  c(lower, upper)
}


#' @keywords internal
.subset_confint_matrix <- function(ci, parm) {
  if (is.null(parm)) {
    return(ci)
  }
  
  if (is.numeric(parm)) {
    return(ci[parm, , drop = FALSE])
  }
  
  if (is.character(parm)) {
    idx <- unique(unlist(lapply(parm, function(p) {
      if (p %in% rownames(ci)) {
        return(match(p, rownames(ci)))
      }
      grep(p, rownames(ci), fixed = TRUE)
    })))
    
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) {
      stop("No matching rows found for `parm`.")
    }
    return(ci[idx, , drop = FALSE])
  }
  
  stop("`parm` must be NULL, numeric, or character.")
}
