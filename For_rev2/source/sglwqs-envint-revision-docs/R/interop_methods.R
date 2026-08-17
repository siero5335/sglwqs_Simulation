#' Tidy Summary for sglwqs Objects
#'
#' Returns model results in a tidy data-frame format compatible with
#' broom-style workflows.
#'
#' @param x A fitted \code{sglwqs} object.
#' @param what Which component to tidy: \code{"weights"} (default),
#'   \code{"coefficients"}, \code{"validation"}, or \code{"bootstrap"}.
#' @param direction Direction filter: \code{"both"} (default),
#'   \code{"positive"}, or \code{"negative"}.
#' @param conf.int Logical. Include confidence intervals when available.
#' @param conf.level Confidence level used when \code{conf.int = TRUE}.
#' @param ... Unused.
#'
#' @return A data frame in tidy format.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:180, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:180],
#'   bootstrap = TRUE,
#'   n_boot = 8,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' generics::tidy(fit, what = "weights")
#' generics::tidy(fit, what = "bootstrap", conf.int = TRUE)
#' }
#'
#' @export
tidy.sglwqs <- function(x,
                        what = c("weights", "coefficients", "validation", "bootstrap"),
                        direction = c("both", "positive", "negative"),
                        conf.int = FALSE,
                        conf.level = 0.95,
                        ...) {
  if (!inherits(x, "sglwqs")) {
    stop("`x` must be a sglwqs object.")
  }
  
  what <- match.arg(what)
  direction <- match.arg(direction)
  
  if (what == "weights") {
    out <- extract_weights(x, direction = direction, by_group = TRUE)
    names(out)[names(out) == "variable"] <- "term"
    return(out)
  }
  
  if (what == "coefficients") {
    coef_df <- .tidy_coefficients_sglwqs(x, direction = direction)
    
    if (conf.int) {
      ci <- tryCatch(
        stats::confint(
          x,
          level = conf.level,
          type = if (x$bootstrap) "bootstrap" else if (x$validation) "validation" else NULL,
          direction = direction
        ),
        error = function(e) NULL
      )
      if (!is.null(ci)) {
        coef_df <- .merge_confint(coef_df, ci)
      }
    }
    
    return(coef_df)
  }
  
  if (what == "validation") {
    if (is.null(.get_active_refit_info(x))) {
      stop("No inference results available. Fit with `validation = TRUE` or `refit = \"full\"` to enable inference.")
    }
    val <- summary_validation(x, conf_level = conf.level)
    out <- data.frame(
      term = val$term,
      estimate = val$estimate,
      std.error = val$std_error,
      statistic = val$estimate / val$std_error,
      p.value = val$p_value,
      direction = val$direction,
      group = val$group,
      stringsAsFactors = FALSE
    )
    if (conf.int) {
      out$conf.low <- val$ci_lower
      out$conf.high <- val$ci_upper
    }
    return(out)
  }
  
  if (!isTRUE(x$bootstrap)) {
    stop("Bootstrap summary requires `bootstrap = TRUE`.")
  }
  boot <- summary_bootstrap(x, conf_level = conf.level, min_freq = 0)
  if (direction != "both") {
    boot <- boot[boot$direction == direction, , drop = FALSE]
  }
  out <- data.frame(
    term = boot$variable,
    direction = boot$direction,
    estimate = boot$mean_coef,
    std.error = boot$se_coef,
    statistic = boot$mean_coef / ifelse(boot$se_coef == 0, NA_real_, boot$se_coef),
    selection_freq = boot$selection_freq,
    weight = boot$weight,
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    out$conf.low <- boot$ci_lower
    out$conf.high <- boot$ci_upper
  }
  out
}


#' Glance Summary for sglwqs Objects
#'
#' Returns one-row model metadata useful for model comparison tables.
#'
#' @param x A fitted \code{sglwqs} object.
#' @param ... Unused.
#'
#' @return A one-row data frame.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:150, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:150],
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' generics::glance(fit)
#' }
#'
#' @export
glance.sglwqs <- function(x, ...) {
  if (!inherits(x, "sglwqs")) {
    stop("`x` must be a sglwqs object.")
  }
  refit_info <- .get_active_refit_info(x)
  source_used <- .get_active_inference_source(x)
  
  data.frame(
    nobs = x$n,
    n_exposures = x$p,
    n_covariates = if (is.null(x$cov_names)) 0 else length(x$cov_names),
    family = x$family,
    lambda = as.character(x$lambda),
    lambda_min = if (!is.null(x$fit$lambda.min)) x$fit$lambda.min else NA_real_,
    lambda_1se = if (!is.null(x$fit$lambda.1se)) x$fit$lambda.1se else NA_real_,
    bootstrap = isTRUE(x$bootstrap),
    n_boot = if (isTRUE(x$bootstrap)) x$boot_info$n_successful else NA_real_,
    validation = identical(source_used, "validation_info"),
    n_train = refit_info$n_train %||% NA_real_,
    n_val = refit_info$n_val %||% NA_real_,
    refit = x$refit,
    refit_engine = if (!is.null(refit_info)) refit_info$engine else NA_character_,
    n_pos_selected = sum(x$pos_weights > 0),
    n_neg_selected = sum(x$neg_weights > 0),
    pos_index_sum = x$pos_index_sum,
    neg_index_sum = x$neg_index_sum,
    stringsAsFactors = FALSE
  )
}


#' Augment Data with sglwqs Predictions
#'
#' Adds fitted values and WQS indices to a dataset.
#'
#' @param x A fitted \code{sglwqs} object.
#' @param data Optional data frame used both as predictors and output base.
#' @param newdata Optional data frame of predictors.
#' @param covariates Optional covariate data for prediction.
#' @param outcome Optional outcome column name (character) or numeric vector
#'   for residual calculation.
#' @param type.predict Prediction scale: \code{"response"} (default) or \code{"link"}.
#' @param ... Unused.
#'
#' @return A data frame with added columns:
#'   \code{.fitted}, \code{.wqs_pos}, \code{.wqs_neg}, and optionally \code{.resid}.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:160, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:160],
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' generics::augment(
#'   fit,
#'   data = sglwqs_example[1:10, exp_vars]
#' )
#' }
#'
#' @export
augment.sglwqs <- function(x, data = NULL, newdata = NULL, covariates = NULL,
                           outcome = NULL,
                           type.predict = c("response", "link"),
                           ...) {
  if (!inherits(x, "sglwqs")) {
    stop("`x` must be a sglwqs object.")
  }
  
  type.predict <- match.arg(type.predict)
  
  if (!is.null(data) && !is.null(newdata)) {
    stop("Provide either `data` or `newdata`, not both.")
  }
  
  base_data <- if (!is.null(data)) {
    data
  } else if (!is.null(newdata)) {
    newdata
  } else {
    data.frame(row_id = seq_len(x$n))
  }
  
  idx <- predict(
    x,
    newdata = if (!is.null(data)) data else newdata,
    covariates = covariates,
    type = "index"
  )
  fitted <- predict(
    x,
    newdata = if (!is.null(data)) data else newdata,
    covariates = covariates,
    type = type.predict
  )
  
  out <- as.data.frame(base_data)
  out$.fitted <- as.numeric(fitted)
  out$.wqs_pos <- idx$wqs_pos
  out$.wqs_neg <- idx$wqs_neg
  
  if (!is.null(outcome)) {
    obs <- if (is.character(outcome) && length(outcome) == 1L) {
      if (!outcome %in% names(out)) {
        stop("`outcome` column not found in `data`/`newdata`.")
      }
      out[[outcome]]
    } else {
      as.numeric(outcome)
    }
    if (length(obs) != nrow(out)) {
      stop("Length of `outcome` does not match rows in augmented data.")
    }
    out$.resid <- as.numeric(obs) - out$.fitted
  }
  
  out
}


#' Validation Metrics for sglwqs Predictions
#'
#' Computes predictive performance metrics on validation or external data.
#'
#' @param object A fitted \code{sglwqs} object.
#' @param newdata Optional new predictor data.
#' @param outcome Optional observed outcome vector, or a column name when
#'   \code{newdata} is a data frame.
#' @param covariates Optional covariate data for prediction.
#' @param threshold Classification threshold for binomial metrics (default: 0.5).
#'
#' @details
#' Survey-aware predictive metrics are not currently implemented. Objects fitted
#' with survey mode (\code{refit_engine = "svyglm"}) will error rather than
#' returning ordinary unweighted metrics.
#'
#' @return A one-row data frame with performance metrics.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit_bin <- sglwqs(
#'   X = sglwqs_example[1:220, exp_vars],
#'   y = sglwqs_example$outcome_bin[1:220],
#'   family = "binomial",
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' validation_metrics(fit_bin)
#' }
#'
#' @export
validation_metrics <- function(object, newdata = NULL, outcome = NULL,
                               covariates = NULL, threshold = 0.5) {
  if (!inherits(object, "sglwqs")) {
    stop("`object` must be a sglwqs object.")
  }
  if (isTRUE(object$survey_mode)) {
    stop(
      "Survey-aware predictive metrics are not yet implemented for survey mode.",
      call. = FALSE
    )
  }
  
  active_fit <- .get_active_glm_fit(object)
  if (is.null(newdata) && is.null(outcome) && !is.null(active_fit)) {
    obs <- stats::model.response(stats::model.frame(active_fit))
    pred <- as.numeric(stats::fitted(active_fit))
  } else {
    pred <- predict(
      object,
      newdata = newdata,
      covariates = covariates,
      type = "response"
    )
    obs <- .resolve_outcome(object, newdata = newdata, outcome = outcome)
  }
  
  if (object$family == "gaussian") {
    rmse <- sqrt(mean((obs - pred)^2))
    mae <- mean(abs(obs - pred))
    sse <- sum((obs - pred)^2)
    sst <- sum((obs - mean(obs))^2)
    r2 <- if (sst > 0) 1 - sse / sst else NA_real_
    
    out <- data.frame(
      family = "gaussian",
      n = length(obs),
      rmse = rmse,
      mae = mae,
      r_squared = r2,
      stringsAsFactors = FALSE
    )
  } else {
    pred <- pmin(pmax(pred, 1e-12), 1 - 1e-12)
    cls <- as.integer(pred >= threshold)
    tp <- sum(cls == 1 & obs == 1)
    tn <- sum(cls == 0 & obs == 0)
    fp <- sum(cls == 1 & obs == 0)
    fn <- sum(cls == 0 & obs == 1)
    
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    accuracy <- mean(cls == obs)
    brier <- mean((obs - pred)^2)
    logloss <- -mean(obs * log(pred) + (1 - obs) * log(1 - pred))
    
    out <- data.frame(
      family = "binomial",
      n = length(obs),
      threshold = threshold,
      auc = .calc_auc(obs, pred),
      brier = brier,
      logloss = logloss,
      accuracy = accuracy,
      sensitivity = sensitivity,
      specificity = specificity,
      precision = precision,
      stringsAsFactors = FALSE
    )
  }
  
  class(out) <- c("sglwqs_validation_metrics", "data.frame")
  out
}


#' Print Validation Metrics
#'
#' @param x An object of class \code{"sglwqs_validation_metrics"}.
#' @param ... Unused.
#'
#' @export
print.sglwqs_validation_metrics <- function(x, ...) {
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}


#' Calibration Table for Binomial sglwqs Models
#'
#' Creates a calibration table by grouping predictions into quantile bins.
#'
#' @param object A fitted \code{sglwqs} object with \code{family = "binomial"}.
#' @param newdata Optional new predictor data.
#' @param outcome Optional observed outcome vector, or a column name when
#'   \code{newdata} is a data frame.
#' @param covariates Optional covariate data for prediction.
#' @param n_bins Number of prediction bins (default: 10).
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Survey-aware calibration is not currently implemented. Objects fitted with
#' survey mode (\code{refit_engine = "svyglm"}) will error rather than
#' returning ordinary unweighted calibration summaries.
#'
#' @return A data frame with bin-level calibration summaries.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit_bin <- sglwqs(
#'   X = sglwqs_example[1:220, exp_vars],
#'   y = sglwqs_example$outcome_bin[1:220],
#'   family = "binomial",
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' cal <- calibrate(fit_bin, n_bins = 6)
#' plot(cal)
#' }
#'
#' @export
calibrate <- function(object, ...) {
  UseMethod("calibrate")
}


#' @rdname calibrate
#' @export
calibrate.sglwqs <- function(object, newdata = NULL, outcome = NULL,
                             covariates = NULL, n_bins = 10, ...) {
  if (!inherits(object, "sglwqs")) {
    stop("`object` must be a sglwqs object.")
  }
  if (isTRUE(object$survey_mode)) {
    stop(
      "Survey-aware calibration is not yet implemented for survey mode.",
      call. = FALSE
    )
  }
  if (object$family != "binomial") {
    stop("Calibration is only defined for binomial models.")
  }

  active_fit <- .get_active_glm_fit(object)
  if (is.null(newdata) && is.null(outcome) && !is.null(active_fit)) {
    obs <- stats::model.response(stats::model.frame(active_fit))
    pred <- as.numeric(stats::fitted(active_fit))
  } else {
    pred <- predict(
      object,
      newdata = newdata,
      covariates = covariates,
      type = "response"
    )
    if (!is.null(newdata) && is.null(outcome)) {
      stop("When `newdata` is supplied, provide the observed `outcome` as well.")
    }
    obs <- .resolve_outcome(object, newdata = newdata, outcome = outcome)
  }

  out <- .build_calibration_table(pred = pred, obs = obs, n_bins = n_bins)
  class(out) <- c("sglwqs_calibration", "data.frame")
  out
}


#' @rdname calibrate
#' @export
calibrate.sglwqs_mids <- function(object, newdata = NULL, outcome = NULL,
                                  covariates = NULL, n_bins = 10, ...) {
  if (!inherits(object, "sglwqs_mids")) {
    stop("`object` must be a sglwqs_mids object.")
  }
  if (object$family != "binomial") {
    stop("Calibration is only defined for binomial models.")
  }

  fits <- Filter(function(f) inherits(f, "sglwqs"), object$fits)
  if (length(fits) == 0) {
    stop("No successful sglwqs fits available for calibration.")
  }

  if (!is.null(newdata) || !is.null(outcome)) {
    if (is.null(newdata) || is.null(outcome)) {
      stop("For sglwqs_mids, provide both `newdata` and `outcome` when calibrating on external data.")
    }

    pred_list <- lapply(fits, function(f) {
      tryCatch(
        as.numeric(predict(f, newdata = newdata, covariates = covariates, type = "response")),
        error = function(e) NULL
      )
    })
    pred_list <- Filter(Negate(is.null), pred_list)
    if (length(pred_list) == 0) {
      stop("Unable to compute pooled predictions from successful imputations.")
    }

    pred_len <- length(pred_list[[1]])
    pred_list <- Filter(function(x) length(x) == pred_len, pred_list)
    if (length(pred_list) == 0) {
      stop("No successful imputations returned prediction vectors of a common length.")
    }

    pred <- rowMeans(do.call(cbind, pred_list), na.rm = TRUE)
    fit_ref <- fits[[1]]
    obs <- .resolve_outcome(fit_ref, newdata = newdata, outcome = outcome)

    out <- .build_calibration_table(pred = pred, obs = obs, n_bins = n_bins)
    class(out) <- c("sglwqs_mids_calibration", "sglwqs_calibration", "data.frame")
    attr(out, "subtitle") <- paste0("Mean predicted probability across ", length(pred_list), " successful imputations")
    return(out)
  }

  val_tables <- lapply(fits, function(f) {
    active_fit <- .get_active_glm_fit(f)
    if (is.null(active_fit)) {
      return(NULL)
    }
    obs <- stats::model.response(stats::model.frame(active_fit))
    pred <- as.numeric(stats::fitted(active_fit))
    .build_calibration_table(pred = pred, obs = obs, n_bins = n_bins)
  })
  val_tables <- Filter(Negate(is.null), val_tables)

  if (length(val_tables) == 0) {
    stop("Provide `newdata` and `outcome`, or fit with `validation = TRUE` or `refit = \"full\"` so calibration can be pooled from in-sample predictions.")
  }

  max_bins <- max(vapply(val_tables, nrow, integer(1)))
  pooled <- do.call(rbind, lapply(seq_len(max_bins), function(i) {
    rows <- lapply(val_tables, function(df) if (nrow(df) >= i) df[i, , drop = FALSE] else NULL)
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0) return(NULL)
    data.frame(
      bin = i,
      n = mean(vapply(rows, function(x) x$n[[1]], numeric(1))),
      pred_mean = mean(vapply(rows, function(x) x$pred_mean[[1]], numeric(1))),
      obs_rate = mean(vapply(rows, function(x) x$obs_rate[[1]], numeric(1))),
      pred_min = mean(vapply(rows, function(x) x$pred_min[[1]], numeric(1))),
      pred_max = mean(vapply(rows, function(x) x$pred_max[[1]], numeric(1))),
      stringsAsFactors = FALSE
    )
  }))
  rownames(pooled) <- NULL
  class(pooled) <- c("sglwqs_mids_calibration", "sglwqs_calibration", "data.frame")
  attr(pooled, "subtitle") <- paste0("Average bin-level calibration across ", length(val_tables), " imputation-specific validation samples")
  attr(pooled, "caption") <- "Each imputation is calibrated on its own validation predictions, then bin summaries are averaged."
  pooled
}


#' Plot Calibration Table
#'
#' @param x An object of class \code{"sglwqs_calibration"}.
#' @param ... Unused.
#'
#' @return A ggplot object.
#'
#' @export
plot.sglwqs_calibration <- function(x, ...) {
  title <- attr(x, "title") %||% "Calibration Plot"
  subtitle <- attr(x, "subtitle")
  caption <- attr(x, "caption")

  ggplot2::ggplot(x, ggplot2::aes(x = .data$pred_mean, y = .data$obs_rate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_point(size = 2, color = "#1F78B4") +
    ggplot2::geom_line(color = "#1F78B4", alpha = 0.8) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = "Mean predicted probability",
      y = "Observed event rate"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.caption = ggplot2::element_text(hjust = 0, size = 9,
                                           color = "gray40", face = "italic")
    )
}


#' Autoplot Method for sglwqs Objects
#'
#' Unified plotting interface for common sglwqs diagnostics.
#'
#' @param object A fitted \code{sglwqs} or \code{sglwqs_mids} object.
#' @param type Plot type: \code{"cv"}, \code{"weights"}, \code{"validation"},
#'   \code{"selection"}, \code{"bootstrap"}, \code{"combined"}, or
#'   \code{"calibration"}.
#' @param ... Additional arguments passed to downstream plotting functions.
#'
#' @return A ggplot object (or a combined plot object when \code{type = "combined"}).
#'
#' @export
autoplot.sglwqs <- function(object,
                            type = c("cv", "weights", "validation", "selection", "bootstrap", "combined", "calibration"),
                            ...) {
  type <- match.arg(type)

  if (type == "cv") return(plot_cv(object, ...))
  if (type == "weights") return(plot_weights(object, ...))
  if (type == "validation") return(plot_validation_results(object, ...))
  if (type == "selection") return(plot_selection_frequency(object, ...))
  if (type == "bootstrap") return(plot_bootstrap_ci(object, ...))
  if (type == "calibration") return(plot(calibrate(object, ...)))
  plot_combined_results(object, ...)
}


#' @rdname autoplot.sglwqs
#' @export
autoplot.sglwqs_mids <- function(object,
                                 type = c("cv", "weights", "validation", "selection", "bootstrap", "combined", "calibration"),
                                 ...) {
  type <- match.arg(type)

  if (type == "cv") return(plot_cv(object, ...))
  if (type == "weights") return(plot_weights(object, ...))
  if (type == "validation") return(plot_validation_results(object, ...))
  if (type == "selection") return(plot_selection_frequency(object, ...))
  if (type == "bootstrap") return(plot_bootstrap_ci(object, ...))
  if (type == "calibration") return(plot(calibrate(object, ...)))
  plot_combined_results(object, ...)
}


#' Convert sglwqs to Data Frame
#'
#' Coerces a fitted object to a data frame for downstream tabular workflows.
#'
#' @param x A fitted \code{sglwqs} object.
#' @param row.names Ignored.
#' @param optional Ignored.
#' @param what Which table to extract: \code{"weights"} (default),
#'   \code{"coefficients"}, \code{"validation"}, or \code{"bootstrap"}.
#' @param direction Direction filter when applicable.
#' @param ... Additional arguments passed to \code{tidy.sglwqs()}.
#'
#' @return A data frame.
#'
#' @export
as.data.frame.sglwqs <- function(x, row.names = NULL, optional = FALSE,
                                 what = c("weights", "coefficients", "validation", "bootstrap"),
                                 direction = c("both", "positive", "negative"),
                                 ...) {
  what <- match.arg(what)
  direction <- match.arg(direction)
  
  tidy.sglwqs(
    x = x,
    what = what,
    direction = direction,
    ...
  )
}


#' Convert sglwqs Output to Tibble
#'
#' Convenience helper for users who prefer tibble output.
#'
#' @param x A fitted \code{sglwqs} object.
#' @param ... Additional arguments passed to \code{as.data.frame.sglwqs()}.
#'
#' @return A tibble.
#'
#' @export
as_tibble_sglwqs <- function(x, ...) {
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required. Install with install.packages('tibble').")
  }
  tibble::as_tibble(as.data.frame.sglwqs(x, ...))
}


# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @keywords internal
.build_calibration_table <- function(pred, obs, n_bins = 10) {
  pred <- as.numeric(pred)
  obs <- as.numeric(obs)

  keep <- is.finite(pred) & is.finite(obs)
  pred <- pred[keep]
  obs <- obs[keep]

  if (length(pred) == 0) {
    stop("No finite predictions and outcomes available for calibration.")
  }
  if (length(pred) != length(obs)) {
    stop("Prediction and outcome vectors must have the same length.")
  }

  n_bins <- max(2L, as.integer(n_bins))
  probs <- seq(0, 1, length.out = n_bins + 1)
  cuts <- unique(stats::quantile(pred, probs = probs, na.rm = TRUE, type = 7))
  if (length(cuts) < 2) {
    cuts <- c(min(pred), max(pred) + 1e-12)
  }

  bin <- cut(pred, breaks = cuts, include.lowest = TRUE, labels = FALSE)
  bins <- split(seq_along(pred), bin)

  out <- do.call(rbind, lapply(seq_along(bins), function(i) {
    idx <- bins[[i]]
    data.frame(
      bin = i,
      n = length(idx),
      pred_mean = mean(pred[idx]),
      obs_rate = mean(obs[idx]),
      pred_min = min(pred[idx]),
      pred_max = max(pred[idx]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}


#' @keywords internal
.tidy_coefficients_sglwqs <- function(x, direction = "both") {
  pos <- data.frame(
    term = names(x$pos_coef),
    estimate = as.numeric(x$pos_coef),
    direction = "positive",
    stringsAsFactors = FALSE
  )
  neg <- data.frame(
    term = names(x$neg_coef),
    estimate = as.numeric(x$neg_coef),
    direction = "negative",
    stringsAsFactors = FALSE
  )
  
  out <- rbind(pos, neg)
  if (direction != "both") {
    out <- out[out$direction == direction, , drop = FALSE]
  }
  
  if (!is.null(x$cov_coef)) {
    cov_df <- data.frame(
      term = names(x$cov_coef),
      estimate = as.numeric(x$cov_coef),
      direction = "covariate",
      stringsAsFactors = FALSE
    )
    out <- rbind(
      data.frame(
        term = "(Intercept)",
        estimate = as.numeric(x$intercept)[1],
        direction = "intercept",
        stringsAsFactors = FALSE
      ),
      out,
      cov_df
    )
  } else {
    out <- rbind(
      data.frame(
        term = "(Intercept)",
        estimate = as.numeric(x$intercept)[1],
        direction = "intercept",
        stringsAsFactors = FALSE
      ),
      out
    )
  }
  
  out$std.error <- NA_real_
  out$statistic <- NA_real_
  out$p.value <- NA_real_
  out
}


#' @keywords internal
.merge_confint <- function(df, ci_matrix) {
  ci_df <- data.frame(
    ci_term = rownames(ci_matrix),
    conf.low = ci_matrix[, 1],
    conf.high = ci_matrix[, 2],
    stringsAsFactors = FALSE
  )
  
  if ("direction" %in% names(df)) {
    df$ci_key <- ifelse(
      df$direction %in% c("positive", "negative"),
      paste0(df$term, " [", df$direction, "]"),
      df$term
    )
  } else {
    df$ci_key <- df$term
  }
  
  idx <- match(df$ci_key, ci_df$ci_term)
  df$conf.low <- ci_df$conf.low[idx]
  df$conf.high <- ci_df$conf.high[idx]
  df$ci_key <- NULL
  df
}


#' @keywords internal
.resolve_outcome <- function(object, newdata = NULL, outcome = NULL) {
  if (!is.null(outcome)) {
    if (is.character(outcome) && length(outcome) == 1L) {
      if (is.null(newdata) || !is.data.frame(newdata) || !outcome %in% names(newdata)) {
        stop("When `outcome` is a column name, it must exist in `newdata`.")
      }
      return(as.numeric(newdata[[outcome]]))
    }
    return(as.numeric(outcome))
  }
  
  active_fit <- .get_active_glm_fit(object)
  if (!is.null(active_fit)) {
    return(stats::model.response(stats::model.frame(active_fit)))
  }
  
  stop(
    "Observed outcome is required. Provide `outcome`, or fit with `validation = TRUE` ",
    "or `refit = \"full\"` and call without `newdata`."
  )
}


#' @keywords internal
.get_active_glm_fit <- function(object) {
  info <- .get_active_refit_info(object)
  if (is.null(info)) {
    return(NULL)
  }
  info$glm_fit %||% info$refit_fit
}


#' @keywords internal
.calc_auc <- function(y_true, p_hat) {
  y_true <- as.integer(y_true)
  pos <- y_true == 1
  neg <- y_true == 0
  
  n_pos <- sum(pos)
  n_neg <- sum(neg)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  
  ranks <- rank(p_hat, ties.method = "average")
  (sum(ranks[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
