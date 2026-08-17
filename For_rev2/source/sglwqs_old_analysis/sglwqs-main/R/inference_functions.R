#' Summary of Bootstrap Results
#'
#' Summarizes bootstrap aggregation results including selection frequencies,
#' coefficient estimates, and confidence intervals.
#' For covariates, bootstrap-aggregated summaries are stored in
#' \code{object$boot_info$mean_cov_coef} and related fields, while
#' \code{object$cov_coef} remains the coefficient vector from the single
#' full-data fit used for backward compatibility.
#'
#' @param object A sglwqs object fitted with bootstrap = TRUE.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @param min_freq Minimum selection frequency to display (default: 0).
#'
#' @return A data frame with bootstrap summary statistics.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:200, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:200],
#'   bootstrap = TRUE,
#'   n_boot = 10,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' summary_bootstrap(fit)
#' summary_bootstrap(fit, min_freq = 0.5)  # Only variables selected > 50% of time
#' }
#'
#' @export
summary_bootstrap <- function(object, conf_level = 0.95, min_freq = 0) {
  
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }
  
  if (!object$bootstrap) {
    stop("Object was not fitted with bootstrap = TRUE")
  }
  
  boot_info <- object$boot_info
  alpha <- 1 - conf_level
  z_val <- qnorm(1 - alpha / 2)
  
  # Check if bootstrap matrices are stored
  has_boot_matrices <- !is.null(boot_info$boot_pos_coef) && !is.null(boot_info$boot_neg_coef)
  
  if (has_boot_matrices) {
    # Compute statistics from bootstrap distribution
    boot_pos <- boot_info$boot_pos_coef
    boot_neg <- boot_info$boot_neg_coef
    
    # Use only valid bootstraps (prefer boot_success if available)
    if (!is.null(boot_info$boot_success)) {
      valid_boots <- boot_info$boot_success
    } else {
      valid_boots <- rowSums(abs(boot_pos) + abs(boot_neg)) > 0
    }
    boot_pos_valid <- boot_pos[valid_boots, , drop = FALSE]
    boot_neg_valid <- boot_neg[valid_boots, , drop = FALSE]

    n_valid <- sum(valid_boots)
    
    # Positive direction
    pos_results <- data.frame(
      variable = object$var_names,
      direction = "positive",
      weight = as.numeric(object$pos_weights),
      mean_coef = colMeans(boot_pos_valid),
      se_coef = apply(boot_pos_valid, 2, sd),
      ci_lower = apply(boot_pos_valid, 2, quantile, probs = alpha / 2),
      ci_upper = apply(boot_pos_valid, 2, quantile, probs = 1 - alpha / 2),
      selection_freq = boot_info$selection_freq_pos,
      stringsAsFactors = FALSE
    )
    
    # Negative direction
    neg_results <- data.frame(
      variable = object$var_names,
      direction = "negative",
      weight = as.numeric(object$neg_weights),
      mean_coef = colMeans(boot_neg_valid),
      se_coef = apply(boot_neg_valid, 2, sd),
      ci_lower = apply(boot_neg_valid, 2, quantile, probs = alpha / 2),
      ci_upper = apply(boot_neg_valid, 2, quantile, probs = 1 - alpha / 2),
      selection_freq = boot_info$selection_freq_neg,
      stringsAsFactors = FALSE
    )
  } else {
    # When bootstrap matrices are unavailable (keep_boot_matrices=FALSE)
    # Compute CI via normal approximation using SE
    message("Note: Bootstrap matrices not stored. Using normal approximation for CI.")
    
    # Positive direction
    pos_results <- data.frame(
      variable = object$var_names,
      direction = "positive",
      weight = as.numeric(object$pos_weights),
      mean_coef = as.numeric(object$pos_coef),
      se_coef = boot_info$se_pos_coef,
      ci_lower = as.numeric(object$pos_coef) - z_val * boot_info$se_pos_coef,
      ci_upper = as.numeric(object$pos_coef) + z_val * boot_info$se_pos_coef,
      selection_freq = boot_info$selection_freq_pos,
      stringsAsFactors = FALSE
    )
    
    # Negative direction
    neg_results <- data.frame(
      variable = object$var_names,
      direction = "negative",
      weight = as.numeric(object$neg_weights),
      mean_coef = as.numeric(object$neg_coef),
      se_coef = boot_info$se_neg_coef,
      ci_lower = as.numeric(object$neg_coef) - z_val * boot_info$se_neg_coef,
      ci_upper = as.numeric(object$neg_coef) + z_val * boot_info$se_neg_coef,
      selection_freq = boot_info$selection_freq_neg,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine
  results <- rbind(pos_results, neg_results)
  
  # Add group info
  if (!is.null(object$groups)) {
    results$group <- .assign_group_labels(results$variable, object$groups)
  }

  # Filtering
  results <- results[results$selection_freq >= min_freq, ]
  
  # Sort (by selection frequency, descending)
 results <- results[order(-results$selection_freq, -abs(results$mean_coef)), ]
  
  # Add meta info as attributes
  attr(results, "n_boot") <- object$boot_info$n_successful
  attr(results, "conf_level") <- conf_level
  if (!is.null(boot_info$mean_cov_coef)) {
    cov_summary <- data.frame(
      term = names(boot_info$mean_cov_coef),
      mean_coef = as.numeric(boot_info$mean_cov_coef),
      se_coef = as.numeric(boot_info$se_cov_coef),
      ci_lower = as.numeric(boot_info$ci_lower_cov),
      ci_upper = as.numeric(boot_info$ci_upper_cov),
      stringsAsFactors = FALSE
    )
    attr(results, "covariate_summary") <- cov_summary
  } else {
    attr(results, "covariate_summary") <- NULL
  }
  
  class(results) <- c("sglwqs_bootstrap_summary", "data.frame")
  
  return(results)
}


#' Print Bootstrap Summary
#'
#' @param x An object of class \code{"sglwqs_bootstrap_summary"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.sglwqs_bootstrap_summary <- function(x, ...) {
  cat("Bootstrap Summary (", attr(x, "n_boot"), " iterations, ", 
      attr(x, "conf_level") * 100, "% CI)\n", sep = "")
  cat(strrep("=", 60), "\n\n")
  
  # Positive direction
  pos_data <- x[x$direction == "positive" & x$selection_freq > 0, ]
  if (nrow(pos_data) > 0) {
    cat("POSITIVE Direction:\n")
    cat(strrep("-", 40), "\n")
    for (i in seq_len(nrow(pos_data))) {
      row <- pos_data[i, ]
      cat(sprintf("  %-20s: weight=%.3f, sel.freq=%.1f%%, coef=%.4f [%.4f, %.4f]\n",
                  row$variable, row$weight, row$selection_freq * 100,
                  row$mean_coef, row$ci_lower, row$ci_upper))
    }
    cat("\n")
  }
  
  # Negative direction
  neg_data <- x[x$direction == "negative" & x$selection_freq > 0, ]
  if (nrow(neg_data) > 0) {
    cat("NEGATIVE Direction:\n")
    cat(strrep("-", 40), "\n")
    for (i in seq_len(nrow(neg_data))) {
      row <- neg_data[i, ]
      cat(sprintf("  %-20s: weight=%.3f, sel.freq=%.1f%%, coef=%.4f [%.4f, %.4f]\n",
                  row$variable, row$weight, row$selection_freq * 100,
                  row$mean_coef, row$ci_lower, row$ci_upper))
    }
  }

  cov_summary <- attr(x, "covariate_summary")
  if (!is.null(cov_summary) && nrow(cov_summary) > 0) {
    cat("\nCOVARIATES (bootstrap-aggregated coefficients):\n")
    cat(strrep("-", 40), "\n")
    for (i in seq_len(nrow(cov_summary))) {
      row <- cov_summary[i, ]
      cat(sprintf("  %-20s: coef=%.4f [%.4f, %.4f], SE=%.4f\n",
                  row$term, row$mean_coef, row$ci_lower, row$ci_upper, row$se_coef))
    }
  }

  if (!is.null(cov_summary) && !is.null(attr(x, "cov_level_note"))) {
    cat("\n", attr(x, "cov_level_note"), "\n", sep = "")
  } else if (!is.null(cov_summary)) {
    cat("\nNote: top-level `cov_coef` remains the single full-data fit coefficient vector; bootstrap-aggregated covariate summaries are shown above.\n")
  }
  
  invisible(x)
}


#' Summary of Validation Results
#'
#' Summarizes validation inference results including estimates, standard errors,
#' confidence intervals, and p-values. When groups are specified, provides
#' group-specific inference.
#'
#' @param object A \code{sglwqs} object with downstream inference enabled via
#'   \code{validation = TRUE} or \code{refit = "full"}.
#' @param conf_level Confidence level for intervals (default: 0.95).
#'
#' @return A data frame with downstream inference summary statistics.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' groups <- list(
#'   Metals = paste0("metal", 1:4),
#'   Pesticides = paste0("pesticide", 1:4)
#' )
#' fit <- sglwqs(
#'   X = sglwqs_example[1:240, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:240],
#'   groups = groups,
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' summary_validation(fit)
#' }
#'
#' @export
summary_validation <- function(object, conf_level = 0.95) {
  
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }
  
  val_info <- .get_active_refit_info(object)
  if (is.null(val_info)) {
    stop("No inference results available. Fit with `validation = TRUE` or `refit = \"full\"` to enable inference.")
  }

  # Use model-aware critical value (t for gaussian/survey, z for binomial)
  refit_fit <- if (!is.null(val_info$refit_fit)) val_info$refit_fit else val_info$glm_fit
  z_val <- if (!is.null(refit_fit)) {
    .prediction_interval_multiplier(refit_fit, conf_level)
  } else {
    stats::qnorm(1 - (1 - conf_level) / 2)
  }

  results_list <- list()
  
  # If group-level results are available
  if (!is.null(val_info$group_results) && length(val_info$group_results) > 0) {
    for (grp_name in names(val_info$group_results)) {
      grp_res <- val_info$group_results[[grp_name]]
      
      # Positive
      results_list[[length(results_list) + 1]] <- data.frame(
        term = paste0(grp_name, " (positive)"),
        group = grp_name,
        direction = "positive",
        estimate = grp_res$pos_estimate,
        std_error = grp_res$pos_se,
        ci_lower = grp_res$pos_estimate - z_val * grp_res$pos_se,
        ci_upper = grp_res$pos_estimate + z_val * grp_res$pos_se,
        p_value = grp_res$pos_pvalue,
        stringsAsFactors = FALSE
      )
      
      # Negative
      results_list[[length(results_list) + 1]] <- data.frame(
        term = paste0(grp_name, " (negative)"),
        group = grp_name,
        direction = "negative",
        estimate = grp_res$neg_estimate,
        std_error = grp_res$neg_se,
        ci_lower = grp_res$neg_estimate - z_val * grp_res$neg_se,
        ci_upper = grp_res$neg_estimate + z_val * grp_res$neg_se,
        p_value = grp_res$neg_pvalue,
        stringsAsFactors = FALSE
      )
    }
  } else {
    # No groups: overall WQS results
    results_list[[1]] <- data.frame(
      term = "WQS Positive",
      group = "Overall",
      direction = "positive",
      estimate = val_info$wqs_pos_estimate,
      std_error = val_info$wqs_pos_se,
      ci_lower = val_info$wqs_pos_estimate - z_val * val_info$wqs_pos_se,
      ci_upper = val_info$wqs_pos_estimate + z_val * val_info$wqs_pos_se,
      p_value = val_info$wqs_pos_pvalue,
      stringsAsFactors = FALSE
    )
    
    results_list[[2]] <- data.frame(
      term = "WQS Negative",
      group = "Overall",
      direction = "negative",
      estimate = val_info$wqs_neg_estimate,
      std_error = val_info$wqs_neg_se,
      ci_lower = val_info$wqs_neg_estimate - z_val * val_info$wqs_neg_se,
      ci_upper = val_info$wqs_neg_estimate + z_val * val_info$wqs_neg_se,
      p_value = val_info$wqs_neg_pvalue,
      stringsAsFactors = FALSE
    )
  }
  
  results <- do.call(rbind, results_list)
  
  # Significance markers
  results$signif <- ""
  results$signif[results$p_value < 0.1] <- "."
  results$signif[results$p_value < 0.05] <- "*"
  results$signif[results$p_value < 0.01] <- "**"
  results$signif[results$p_value < 0.001] <- "***"
  
  # Meta info
  attr(results, "n_train") <- val_info$n_train %||% NA_integer_
  attr(results, "n_val") <- val_info$n_val %||% NA_integer_
  attr(results, "conf_level") <- conf_level
  attr(results, "family") <- object$family
  attr(results, "group_inference") <- !is.null(val_info$group_results) && length(val_info$group_results) > 0
  attr(results, "formula") <- val_info$formula
  attr(results, "source") <- .get_active_inference_source(object)
  
  class(results) <- c("sglwqs_validation_summary", "data.frame")
  
  return(results)
}


#' Print Validation Summary
#'
#' @param x An object of class \code{"sglwqs_validation_summary"}.
#' @param ... Additional arguments passed to \code{print()}.
#' @export
print.sglwqs_validation_summary <- function(x, ...) {
  source_label <- if (identical(attr(x, "source"), "refit_info")) "Refit" else "Validation"
  if (identical(attr(x, "source"), "boot_info")) {
    source_label <- "Bootstrap"
  }
  cat(source_label, " Inference Summary\n", sep = "")
  cat(strrep("=", 60), "\n")
  if (!is.na(attr(x, "n_train")) || !is.na(attr(x, "n_val"))) {
    cat("Training set: ", attr(x, "n_train"), " | Validation set: ",
        attr(x, "n_val"), "\n", sep = "")
  }
  cat("Family: ", attr(x, "family"), " | CI level: ", 
      attr(x, "conf_level") * 100, "%\n", sep = "")
  
  if (!is.null(attr(x, "formula"))) {
    cat("Model: ", attr(x, "formula"), "\n", sep = "")
  }
  cat("\n")
  
  # For group-specific inference
  if (isTRUE(attr(x, "group_inference"))) {
    groups <- unique(x$group)
    for (grp in groups) {
      cat("--- ", grp, " ---\n", sep = "")
      grp_data <- x[x$group == grp, ]
      for (i in seq_len(nrow(grp_data))) {
        row <- grp_data[i, ]
        dir_label <- ifelse(row$direction == "positive", "Positive", "Negative")
        cat(sprintf("  %s: estimate=%8.4f (SE=%.4f), p=%s %s\n", 
                    dir_label, row$estimate, row$std_error,
                    format.pval(row$p_value, digits = 3), row$signif))
        cat(sprintf("           %d%% CI: [%.4f, %.4f]\n",
                    attr(x, "conf_level") * 100, row$ci_lower, row$ci_upper))
      }
      cat("\n")
    }
  } else {
    # Overall results
    for (i in seq_len(nrow(x))) {
      row <- x[i, ]
      cat(sprintf("%s:\n", row$term))
      cat(sprintf("  Estimate: %8.4f (SE: %.4f)\n", row$estimate, row$std_error))
      cat(sprintf("  %d%% CI:   [%.4f, %.4f]\n", 
                  attr(x, "conf_level") * 100, row$ci_lower, row$ci_upper))
      cat(sprintf("  P-value:  %s %s\n\n", format.pval(row$p_value, digits = 4), row$signif))
    }
  }
  
  cat("Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n")
  
  invisible(x)
}


#' Summary of Available Inference Results
#'
#' Information-source agnostic summary of available inference results.
#'
#' @param object A fitted \code{sglwqs} or \code{sglwqs_mids} object.
#' @param conf_level Confidence level for intervals (default: 0.95).
#' @export
summary_inference <- function(object, conf_level = 0.95) {
  if (inherits(object, "sglwqs_mids")) {
    return(summary_inference.sglwqs_mids(object, conf_level = conf_level))
  }
  if (inherits(object, "sglwqs")) {
    return(summary_inference.sglwqs(object, conf_level = conf_level))
  }
  stop("Object must be of class 'sglwqs' or 'sglwqs_mids'.")
}


#' @rdname summary_inference
#' @export
summary_inference.sglwqs <- function(object, conf_level = 0.95) {
  if (!is.null(.get_active_refit_info(object))) {
    return(summary_validation(object, conf_level = conf_level))
  }

  .summary_bootstrap_inference_sglwqs(object, conf_level = conf_level)
}


#' @keywords internal
.bootstrap_p_value <- function(estimate, se) {
  if (is.na(estimate) || is.na(se) || se < 0) {
    return(NA_real_)
  }
  if (se == 0) {
    return(if (estimate == 0) NA_real_ else 0)
  }
  2 * stats::pnorm(-abs(estimate / se))
}


#' @keywords internal
.build_bootstrap_inference_summary <- function(group_results = NULL,
                                               wqs_pos = NULL,
                                               wqs_neg = NULL,
                                               covariates = NULL,
                                               conf_level = 0.95,
                                               family = NA_character_,
                                               source = "boot_info") {
  alpha <- 1 - conf_level
  z_val <- stats::qnorm(1 - alpha / 2)
  rows <- list()
  summarize_boot_term <- function(term, group, direction, type, res) {
    est <- res$estimate %||% NA_real_
    se <- res$se %||% NA_real_
    data.frame(
      term = term,
      group = group,
      direction = direction,
      estimate = est,
      std_error = se,
      ci_lower = res$ci_lower %||% if (is.finite(est) && is.finite(se)) est - z_val * se else NA_real_,
      ci_upper = res$ci_upper %||% if (is.finite(est) && is.finite(se)) est + z_val * se else NA_real_,
      p_value = res$p_value %||% .bootstrap_p_value(est, se),
      df = res$df %||% NA_real_,
      fmi = res$fmi %||% NA_real_,
      t_stat = res$t_stat %||% NA_real_,
      type = type,
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(group_results) && length(group_results) > 0) {
    for (grp_name in names(group_results)) {
      grp_res <- group_results[[grp_name]]
      for (dir_name in c("positive", "negative")) {
        rows[[length(rows) + 1]] <- summarize_boot_term(
          term = paste0(grp_name, " (", dir_name, ")"),
          group = grp_name,
          direction = dir_name,
          type = "WQS Index",
          res = grp_res[[dir_name]]
        )
      }
    }
  } else {
    if (!is.null(wqs_pos)) {
      rows[[length(rows) + 1]] <- summarize_boot_term(
        term = "WQS Positive",
        group = "Overall",
        direction = "positive",
        type = "WQS Index",
        res = wqs_pos
      )
    }
    if (!is.null(wqs_neg)) {
      rows[[length(rows) + 1]] <- summarize_boot_term(
        term = "WQS Negative",
        group = "Overall",
        direction = "negative",
        type = "WQS Index",
        res = wqs_neg
      )
    }
  }

  if (!is.null(covariates) && length(covariates) > 0) {
    for (cov_name in names(covariates)) {
      rows[[length(rows) + 1]] <- summarize_boot_term(
        term = cov_name,
        group = "Covariate",
        direction = NA_character_,
        type = "Covariate",
        res = covariates[[cov_name]]
      )
    }
  }

  if (length(rows) == 0) {
    stop("No inference results available. Fit with `validation = TRUE`, `refit = \"full\"`, or `bootstrap = TRUE` to enable inference.")
  }

  out <- do.call(rbind, rows)
  out$signif <- ""
  out$signif[!is.na(out$p_value) & out$p_value < 0.1] <- "."
  out$signif[!is.na(out$p_value) & out$p_value < 0.05] <- "*"
  out$signif[!is.na(out$p_value) & out$p_value < 0.01] <- "**"
  out$signif[!is.na(out$p_value) & out$p_value < 0.001] <- "***"
  attr(out, "n_train") <- NA_integer_
  attr(out, "n_val") <- NA_integer_
  attr(out, "conf_level") <- conf_level
  attr(out, "family") <- family
  attr(out, "group_inference") <- !is.null(group_results) && length(group_results) > 0
  attr(out, "formula") <- NULL
  attr(out, "has_df") <- any(is.finite(out$df))
  attr(out, "has_fmi") <- any(is.finite(out$fmi))
  attr(out, "source") <- source
  class(out) <- c("sglwqs_validation_summary", "data.frame")
  out
}


#' @keywords internal
.summary_bootstrap_inference_sglwqs <- function(object, conf_level = 0.95) {
  bi <- object$boot_info
  if (is.null(bi)) {
    stop("No inference results available. Fit with `validation = TRUE`, `refit = \"full\"`, or `bootstrap = TRUE` to enable inference.")
  }

  if (!is.null(object$groups) && !is.null(bi$mean_index_sum_by_group_pos)) {
    group_results <- lapply(names(object$groups), function(grp) {
      list(
        positive = list(
          estimate = bi$mean_index_sum_by_group_pos[[grp]] %||% NA_real_,
          se = bi$se_index_sum_by_group_pos[[grp]] %||% NA_real_
        ),
        negative = list(
          estimate = bi$mean_index_sum_by_group_neg[[grp]] %||% NA_real_,
          se = bi$se_index_sum_by_group_neg[[grp]] %||% NA_real_
        )
      )
    })
    names(group_results) <- names(object$groups)
    covariates <- NULL
    if (!is.null(bi$mean_cov_coef)) {
      covariates <- lapply(names(bi$mean_cov_coef), function(nm) {
        list(
          estimate = bi$mean_cov_coef[[nm]],
          se = bi$se_cov_coef[[nm]] %||% NA_real_
        )
      })
      names(covariates) <- names(bi$mean_cov_coef)
    }
    return(.build_bootstrap_inference_summary(
      group_results = group_results,
      covariates = covariates,
      conf_level = conf_level,
      family = object$family,
      source = "boot_info"
    ))
  }

  covariates <- NULL
  if (!is.null(bi$mean_cov_coef)) {
    covariates <- lapply(names(bi$mean_cov_coef), function(nm) {
      list(
        estimate = bi$mean_cov_coef[[nm]],
        se = bi$se_cov_coef[[nm]] %||% NA_real_
      )
    })
    names(covariates) <- names(bi$mean_cov_coef)
  }

  .build_bootstrap_inference_summary(
    wqs_pos = list(estimate = bi$mean_index_sum_pos %||% object$pos_index_sum,
                   se = bi$se_index_sum_pos %||% NA_real_),
    wqs_neg = list(estimate = bi$mean_index_sum_neg %||% object$neg_index_sum,
                   se = bi$se_index_sum_neg %||% NA_real_),
    covariates = covariates,
    conf_level = conf_level,
    family = object$family,
    source = "boot_info"
  )
}


#' @rdname summary_inference
#' @export
summary_inference.sglwqs_mids <- function(object, conf_level = 0.95) {
  inf <- .get_mi_pooled_inference(object)
  if (!is.null(inf)) {
    return(.build_mi_validation_summary_table(object, conf_level = conf_level))
  }

  boot_inf <- object$pooled$bootstrap_inference
  if (is.null(boot_inf)) {
    stop("No inference results available. Fit with `validation = TRUE`, `refit = \"full\"`, or `bootstrap = TRUE` to enable inference.")
  }

  .build_bootstrap_inference_summary(
    group_results = boot_inf$group_results,
    wqs_pos = boot_inf$wqs_pos,
    wqs_neg = boot_inf$wqs_neg,
    covariates = boot_inf$covariates,
    conf_level = conf_level,
    family = object$family,
    source = boot_inf$source_used %||% "boot_info"
  )
}


#' Plot Selection Frequency
#'
#' Generic function to visualize selection frequencies for fitted
#' \code{sglwqs} and \code{sglwqs_mids} objects.
#'
#' For \code{sglwqs} objects, this reproduces the bootstrap-frequency bar plot.
#' For \code{sglwqs_mids} objects, a pooled butterfly plot of MI or pooled
#' bootstrap selection frequency is returned.
#'
#' @param object A fitted \code{sglwqs} or \code{sglwqs_mids} object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' groups <- list(
#'   Metals = paste0("metal", 1:4),
#'   Pesticides = paste0("pesticide", 1:4)
#' )
#' fit <- sglwqs(
#'   X = sglwqs_example[1:200, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:200],
#'   groups = groups,
#'   bootstrap = TRUE,
#'   n_boot = 10,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' plot_selection_frequency(fit)
#' plot_selection_frequency(fit, direction = "positive", top_n = 10)
#' plot_selection_frequency(fit, facet_by_group = TRUE, sort_by = "frequency")
#' }
#'
#' @export
plot_selection_frequency <- function(object, ...) {
  UseMethod("plot_selection_frequency")
}


#' @rdname plot_selection_frequency
#' @param top_n Number of top variables to show per group/direction (default: 15).
#' @param min_freq Minimum selection frequency to display (default: 0.1).
#' @param direction Which direction to plot: "both", "positive", or "negative".
#' @param sort_by How to sort variables: "frequency" (default), "group", or "name".
#' @param facet_by_group Whether to create separate panels for each chemical group (default: TRUE if groups exist).
#' @param show_threshold Whether to show the 50% threshold line (default: TRUE).
#' @export
plot_selection_frequency.sglwqs <- function(object, top_n = 15, min_freq = 0.1,
                                            direction = c("both", "positive", "negative"),
                                            sort_by = c("frequency", "group", "name"),
                                            facet_by_group = NULL,
                                            show_threshold = TRUE,
                                            ...) {

  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }

  if (!object$bootstrap) {
    stop("Object was not fitted with bootstrap = TRUE. Selection frequency requires bootstrap.")
  }

  direction <- match.arg(direction)
  sort_by <- match.arg(sort_by)

  if (is.null(facet_by_group)) {
    facet_by_group <- !is.null(object$groups)
  }

  pos_df <- data.frame(
    variable = object$var_names,
    selection_freq = object$boot_info$selection_freq_pos,
    direction = "Positive",
    stringsAsFactors = FALSE
  )

  neg_df <- data.frame(
    variable = object$var_names,
    selection_freq = object$boot_info$selection_freq_neg,
    direction = "Negative",
    stringsAsFactors = FALSE
  )

  if (!is.null(object$groups)) {
    pos_df$group <- .assign_group_labels(pos_df$variable, object$groups)
    neg_df$group <- .assign_group_labels(neg_df$variable, object$groups)
    group_levels <- .group_display_order(object$groups, unique(c(pos_df$group, neg_df$group)))
    pos_df$group <- factor(pos_df$group, levels = group_levels)
    neg_df$group <- factor(neg_df$group, levels = group_levels)
  } else {
    pos_df$group <- "All"
    neg_df$group <- "All"
  }

  if (direction == "positive") {
    plot_data <- pos_df
  } else if (direction == "negative") {
    plot_data <- neg_df
  } else {
    plot_data <- rbind(pos_df, neg_df)
  }

  plot_data <- plot_data[plot_data$selection_freq >= min_freq, ]

  if (nrow(plot_data) == 0) {
    message("No variables with selection frequency >= ", min_freq * 100, "%")
    return(invisible(NULL))
  }

  if (facet_by_group && !is.null(object$groups)) {
    selected_data <- do.call(rbind, lapply(split(plot_data, list(plot_data$group, plot_data$direction)), function(df) {
      if (nrow(df) == 0) return(NULL)
      df <- df[order(-df$selection_freq), ]
      head(df, top_n)
    }))
    plot_data <- selected_data
  } else {
    if (direction == "both") {
      pos_top <- head(pos_df[order(-pos_df$selection_freq), ], top_n)
      neg_top <- head(neg_df[order(-neg_df$selection_freq), ], top_n)
      top_vars <- unique(c(pos_top$variable, neg_top$variable))
      plot_data <- plot_data[plot_data$variable %in% top_vars, ]
    } else {
      plot_data <- head(plot_data[order(-plot_data$selection_freq), ], top_n)
    }
  }

  if (nrow(plot_data) == 0) {
    message("No variables to plot after filtering")
    return(invisible(NULL))
  }

  if (sort_by == "frequency") {
    if (facet_by_group && !is.null(object$groups)) {
      plot_data <- plot_data[order(plot_data$group, -plot_data$selection_freq), ]
      plot_data$var_label <- paste0(plot_data$group, "_", plot_data$variable)
      var_order <- unique(plot_data$var_label)
      plot_data$var_label <- factor(plot_data$var_label, levels = rev(var_order))
    } else {
      plot_data <- plot_data[order(-plot_data$selection_freq), ]
      var_order <- unique(plot_data$variable)
      plot_data$variable <- factor(plot_data$variable, levels = rev(var_order))
    }
  } else if (sort_by == "group") {
    plot_data <- plot_data[order(plot_data$group, plot_data$variable), ]
    var_order <- unique(plot_data$variable)
    plot_data$variable <- factor(plot_data$variable, levels = rev(var_order))
  } else {
    plot_data <- plot_data[order(plot_data$variable), ]
    var_order <- unique(plot_data$variable)
    plot_data$variable <- factor(plot_data$variable, levels = rev(var_order))
  }

  if (facet_by_group && !is.null(object$groups) && sort_by == "frequency") {
    plot_data$display_name <- gsub("^[^_]+_", "", as.character(plot_data$var_label))

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = stats::reorder(display_name, selection_freq), 
                                                 y = selection_freq * 100))
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = selection_freq * 100))
  }

  p <- p +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Bootstrap Selection Frequency",
      subtitle = paste0("Based on ", object$boot_info$n_successful, " bootstrap iterations",
                        " (min freq: ", min_freq * 100, "%)"),
      x = NULL,
      y = "Selection Frequency (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))

  if (!is.null(object$groups) && !facet_by_group) {
    p <- p + 
      ggplot2::geom_bar(ggplot2::aes(fill = group), stat = "identity", width = 0.7) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Chemical Group")
  } else if (direction == "both" && !facet_by_group) {
    p <- p + 
      ggplot2::geom_bar(ggplot2::aes(fill = direction), stat = "identity", width = 0.7) +
      ggplot2::scale_fill_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B"),
                                 name = "Direction")
  } else {
    p <- p + ggplot2::geom_bar(stat = "identity", fill = "steelblue", width = 0.7)
  }

  if (facet_by_group && !is.null(object$groups)) {
    if (direction == "both") {
      p <- p + ggplot2::facet_grid(group ~ direction, scales = "free_y", space = "free_y")
    } else {
      p <- p + ggplot2::facet_wrap(~ group, scales = "free_y", ncol = 1)
    }
  } else if (direction == "both") {
    p <- p + ggplot2::facet_wrap(~ direction, scales = "free_y")
  }

  if (show_threshold) {
    p <- p + ggplot2::geom_hline(yintercept = 50, linetype = "dashed", color = "red", alpha = 0.5)
  }

  return(p)
}


#' Plot Bootstrap Confidence Intervals
#'
#' Generic function to visualize coefficient estimates with confidence intervals
#' for fitted \code{sglwqs} and \code{sglwqs_mids} objects.
#'
#' @param object A fitted \code{sglwqs} or \code{sglwqs_mids} object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A ggplot object.
#'
#' @export
plot_bootstrap_ci <- function(object, ...) {
  UseMethod("plot_bootstrap_ci")
}


#' @rdname plot_bootstrap_ci
#' @param top_n Number of top variables to show (default: 15).
#' @param conf_level Confidence level (default: 0.95).
#' @param direction Which direction to plot: "both", "positive", or "negative".
#' @param pos_color Color for positive direction.
#' @param neg_color Color for negative direction.
#' @param base_size Base font size.
#' @param title Optional plot title. When \code{NULL}, a default title is used.
#' @export
plot_bootstrap_ci.sglwqs <- function(object, top_n = 15, conf_level = 0.95,
                                     direction = c("both", "positive", "negative"),
                                     pos_color = "#2166AC",
                                     neg_color = "#B2182B",
                                     base_size = 11,
                                     title = NULL,
                                     ...) {

  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }

  if (!object$bootstrap) {
    stop("Object was not fitted with bootstrap = TRUE")
  }

  direction <- match.arg(direction)
  boot_summary <- summary_bootstrap(object, conf_level = conf_level, min_freq = 0)

  if (direction == "positive") {
    plot_data <- boot_summary[boot_summary$direction == "positive", ]
  } else if (direction == "negative") {
    plot_data <- boot_summary[boot_summary$direction == "negative", ]
  } else {
    plot_data <- boot_summary
  }

  plot_data <- plot_data[abs(plot_data$mean_coef) > 1e-6, ]

  if (nrow(plot_data) == 0) {
    message("No non-zero coefficients to plot")
    return(invisible(NULL))
  }

  plot_data <- head(plot_data[order(-abs(plot_data$mean_coef)), ], top_n)
  plot_data$variable <- factor(plot_data$variable, 
                               levels = rev(unique(plot_data$variable[order(plot_data$mean_coef)])))

  if (is.null(title)) {
    title <- paste0("Bootstrap Coefficient Estimates (", conf_level * 100, "% CI)")
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = variable, y = mean_coef)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lower, ymax = ci_upper, color = direction),
                             size = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Based on ", object$boot_info$n_successful, " bootstrap iterations"),
      x = NULL,
      y = "Coefficient"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::scale_color_manual(values = c("positive" = pos_color, "negative" = neg_color),
                                name = "Direction") +
    ggplot2::theme(legend.position = "bottom")

  if (direction == "both") {
    p <- p + ggplot2::facet_wrap(~ direction, scales = "free")
  }

  return(p)
}


#' @rdname plot_validation_results
#' @param conf_level Confidence level (default: 0.95).
#' @param include_covariates Whether to include covariate coefficients (default: FALSE).
#' @export
plot_validation_results.sglwqs <- function(object, conf_level = 0.95, include_covariates = FALSE, ...) {
  
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }
  
  val_info <- .get_active_refit_info(object)
  if (is.null(val_info)) {
    stop("No inference results available. Fit with `validation = TRUE` or `refit = \"full\"` to enable inference.")
  }

  # Model-aware critical value
  refit_fit_pv <- if (!is.null(val_info$refit_fit)) val_info$refit_fit else val_info$glm_fit
  z_val <- if (!is.null(refit_fit_pv)) {
    .prediction_interval_multiplier(refit_fit_pv, conf_level)
  } else {
    stats::qnorm(1 - (1 - conf_level) / 2)
  }

  # Get data using summary_validation
  val_summary <- summary_validation(object, conf_level = conf_level)
  
  # Create plot data
  plot_data <- data.frame(
    term = val_summary$term,
    estimate = val_summary$estimate,
    se = val_summary$std_error,
    ci_lower = val_summary$ci_lower,
    ci_upper = val_summary$ci_upper,
    p_value = val_summary$p_value,
    group = val_summary$group,
    direction = val_summary$direction,
    type = "WQS Index",
    stringsAsFactors = FALSE
  )
  
  # Identify NA effects and generate caption
  na_effects <- plot_data[is.na(plot_data$estimate), ]
  caption_text <- NULL
  
  if (nrow(na_effects) > 0) {
    # Group by direction
    pos_na_groups <- na_effects$group[na_effects$direction == "positive"]
    neg_na_groups <- na_effects$group[na_effects$direction == "negative"]
    
    caption_parts <- c()
    
    if (length(pos_na_groups) > 0) {
      if (length(pos_na_groups) == 1) {
        caption_parts <- c(caption_parts, 
          paste0("Positive effects not shown in ", pos_na_groups))
      } else {
        caption_parts <- c(caption_parts,
          paste0("Positive effects not shown in ", 
                 paste(pos_na_groups, collapse = " and ")))
      }
    }
    
    if (length(neg_na_groups) > 0) {
      if (length(neg_na_groups) == 1) {
        caption_parts <- c(caption_parts,
          paste0("Negative effects not shown in ", neg_na_groups))
      } else {
        caption_parts <- c(caption_parts,
          paste0("Negative effects not shown in ",
                 paste(neg_na_groups, collapse = " and ")))
      }
    }
    
    if (length(caption_parts) > 0) {
      caption_text <- paste0(
        paste(caption_parts, collapse = "; "),
        "\n(coefficients = 0 due to regularization)"
      )
    }
  }
  
  # When including covariates
  if (include_covariates && !is.null(object$cov_coef)) {
    coef_table <- val_info$coef_table
    cov_names <- object$cov_names
    p_col <- val_info$p_col
    
    for (cov_name in cov_names) {
      if (cov_name %in% rownames(coef_table)) {
        est <- coef_table[cov_name, "Estimate"]
        se <- coef_table[cov_name, "Std. Error"]
        cov_row <- data.frame(
          term = cov_name,
          estimate = est,
          se = se,
          ci_lower = est - z_val * se,
          ci_upper = est + z_val * se,
          p_value = coef_table[cov_name, p_col],
          group = "Covariate",
          direction = NA,
          type = "Covariate",
          stringsAsFactors = FALSE
        )
        plot_data <- rbind(plot_data, cov_row)
      }
    }
  }
  
  # Significance
  plot_data$significant <- plot_data$p_value < 0.05
  
  # Plot data excluding NA
  plot_data_valid <- plot_data[!is.na(plot_data$estimate), ]

  # Set variable order
  plot_data_valid$term <- factor(plot_data_valid$term, levels = rev(unique(plot_data_valid$term)))
  
  # For group-specific inference
  if (isTRUE(attr(val_summary, "group_inference"))) {
    # Color-code by group
    p <- ggplot2::ggplot(plot_data_valid, ggplot2::aes(x = term, y = estimate)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lower, ymax = ci_upper, 
                                             color = group, shape = direction),
                               size = 1, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Validation Inference Results (by Group)",
        subtitle = paste0("Training: ", val_info$n_train, " | Validation: ", val_info$n_val,
                          " | ", conf_level * 100, "% CI"),
        x = NULL,
        y = "Estimate",
        caption = caption_text
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_brewer(palette = "Set1", name = "Group") +
      ggplot2::scale_shape_manual(values = c("positive" = 16, "negative" = 17),
                                   name = "Direction",
                                   na.value = 15) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.caption = ggplot2::element_text(hjust = 0, size = 9, 
                                              color = "gray40", face = "italic")
      )
  } else {
    # Overall WQS
    p <- ggplot2::ggplot(plot_data_valid, ggplot2::aes(x = term, y = estimate)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_pointrange(ggplot2::aes(ymin = ci_lower, ymax = ci_upper, 
                                             color = significant, shape = type),
                               size = 1, na.rm = TRUE) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = "Validation Inference Results",
        subtitle = paste0("Training: ", val_info$n_train, " | Validation: ", val_info$n_val,
                          " | ", conf_level * 100, "% CI"),
        x = NULL,
        y = "Estimate",
        caption = caption_text
      ) +
      ggplot2::theme_minimal() +
      ggplot2::scale_color_manual(values = c("TRUE" = "#D7191C", "FALSE" = "gray50"),
                                   name = "p < 0.05",
                                   labels = c("TRUE" = "Yes", "FALSE" = "No")) +
      ggplot2::scale_shape_manual(values = c("WQS Index" = 16, "Covariate" = 17),
                                   name = "Type") +
      ggplot2::theme(
        legend.position = "bottom",
        plot.caption = ggplot2::element_text(hjust = 0, size = 9, 
                                              color = "gray40", face = "italic")
      )
  }
  
  # Annotate p-values
  p <- p + ggplot2::geom_text(
    ggplot2::aes(label = sprintf("p=%s", format.pval(p_value, digits = 2))),
    hjust = -0.2, size = 3, na.rm = TRUE
  )
  
  return(p)
}

#' Plot Available Inference Results
#'
#' Information-source agnostic plot of available inference results.
#'
#' @inheritParams plot_validation_results
#' @param object A fitted \code{sglwqs} or \code{sglwqs_mids} object with
#'   downstream inference (\code{validation = TRUE} or \code{refit = "full"})
#'   or bootstrap-only inference (\code{bootstrap = TRUE}).
#' @export
plot_inference_results <- function(object, ...) {
  if (inherits(object, "sglwqs_mids")) {
    return(plot_inference_results.sglwqs_mids(object, ...))
  }
  if (inherits(object, "sglwqs")) {
    return(plot_inference_results.sglwqs(object, ...))
  }
  stop("Object must be of class 'sglwqs' or 'sglwqs_mids'.")
}


#' @keywords internal
.plot_inference_summary_df <- function(summary_df,
                                       title = "Inference Results",
                                       pos_color = "#0072B2",
                                       neg_color = "#E69F00",
                                       show_significance = c("none", "asterisk", "pvalue"),
                                       base_size = 11) {
  show_significance <- match.arg(show_significance)
  plot_data <- summary_df
  plot_data_valid <- plot_data[!is.na(plot_data$estimate), , drop = FALSE]
  if (nrow(plot_data_valid) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = "No inference effects available", size = 4) +
        ggplot2::theme_void() +
        ggplot2::labs(title = title)
    )
  }

  plot_data_valid$term <- factor(plot_data_valid$term, levels = rev(unique(plot_data_valid$term)))
  plot_data_valid$significant <- plot_data_valid$p_value < 0.05
  plot_data_valid$color_group <- ifelse(
    plot_data_valid$type == "Covariate", "covariate",
    ifelse(plot_data_valid$direction == "positive", "positive", "negative")
  )

  if (show_significance == "asterisk") {
    plot_data_valid$sig_label <- plot_data_valid$signif
  } else if (show_significance == "pvalue") {
    plot_data_valid$sig_label <- ifelse(is.na(plot_data_valid$p_value), "", format.pval(plot_data_valid$p_value, digits = 2))
  } else {
    plot_data_valid$sig_label <- ""
  }

  p <- ggplot2::ggplot(plot_data_valid, ggplot2::aes(x = .data$term, y = .data$estimate)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_pointrange(
      ggplot2::aes(
        ymin = .data$ci_lower,
        ymax = .data$ci_upper,
        color = .data$color_group,
        shape = .data$type
      ),
      size = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = title,
      subtitle = paste0(attr(summary_df, "source") %||% "inference", " | ", attr(summary_df, "conf_level") * 100, "% CI"),
      x = NULL,
      y = "Estimate"
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::scale_color_manual(
      values = c(positive = pos_color, negative = neg_color, covariate = "gray40"),
      name = "Direction"
    ) +
    ggplot2::scale_shape_manual(values = c("WQS Index" = 16, "Covariate" = 17), name = "Type") +
    ggplot2::theme(legend.position = "bottom")

  if (show_significance != "none") {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$sig_label),
      hjust = -0.3,
      size = 3,
      na.rm = TRUE
    )
  }

  p
}


#' @rdname plot_inference_results
#' @export
plot_inference_results.sglwqs <- function(object,
                                          conf_level = 0.95,
                                          include_covariates = FALSE,
                                          pos_color = "#0072B2",
                                          neg_color = "#E69F00",
                                          show_significance = c("none", "asterisk", "pvalue"),
                                          base_size = 11,
                                          title = NULL,
                                          ...) {
  if (!inherits(object, "sglwqs")) {
    stop("Object must be of class 'sglwqs'")
  }

  if (!is.null(.get_active_refit_info(object))) {
    return(plot_validation_results(
      object,
      conf_level = conf_level,
      include_covariates = include_covariates,
      ...
    ))
  }

  summary_df <- summary_inference(object, conf_level = conf_level)
  if (!isTRUE(include_covariates)) {
    summary_df <- summary_df[summary_df$type != "Covariate", , drop = FALSE]
  }
  .plot_inference_summary_df(
    summary_df,
    title = title %||% "Bootstrap Inference Results",
    pos_color = pos_color,
    neg_color = neg_color,
    show_significance = show_significance,
    base_size = base_size
  )
}


#' @rdname plot_inference_results
#' @export
plot_inference_results.sglwqs_mids <- function(object,
                                               conf_level = 0.95,
                                               include_covariates = FALSE,
                                               pos_color = "#0072B2",
                                               neg_color = "#E69F00",
                                               show_significance = c("none", "asterisk", "pvalue"),
                                               base_size = 11,
                                               title = NULL,
                                               ...) {
  if (!inherits(object, "sglwqs_mids")) {
    stop("Object must be of class 'sglwqs_mids'")
  }

  if (!is.null(.get_mi_pooled_inference(object))) {
    return(plot_validation_results(
      object,
      conf_level = conf_level,
      include_covariates = include_covariates,
      pos_color = pos_color,
      neg_color = neg_color,
      show_significance = show_significance,
      base_size = base_size,
      title = title,
      ...
    ))
  }

  summary_df <- summary_inference(object, conf_level = conf_level)
  if (!isTRUE(include_covariates)) {
    summary_df <- summary_df[summary_df$type != "Covariate", , drop = FALSE]
  }
  .plot_inference_summary_df(
    summary_df,
    title = title %||% "Pooled Bootstrap Inference Results",
    pos_color = pos_color,
    neg_color = neg_color,
    show_significance = show_significance,
    base_size = base_size
  )
}


#' @rdname plot_combined_results
#' @export
plot_combined_results.sglwqs <- function(object,
                                   top_n = 10,
                                   sort_by = c("combined", "selection_freq", "weight"),
                                   pos_color = "#0072B2",
                                   neg_color = "#E69F00",
                                   layout = c("auto", "vertical", "horizontal"),
                                   show_significance = c("none", "asterisk", "pvalue"),
                                   conf_level = 0.95,
                                   base_size = 11,
                                   include_covariates = FALSE,
                                   ...) {
  
  sort_by <- match.arg(sort_by)
  layout <- match.arg(layout)
  show_significance <- match.arg(show_significance)
  
  # Determine common variable order and group information
  var_info <- .compute_unified_var_order(object, sort_by = sort_by, top_n = top_n)
  
  plots <- list()
  
  # A) Selection frequency plot (butterfly format)
  if (object$bootstrap) {
    plots$selection <- .plot_selection_butterfly(
      object, var_info = var_info,
      pos_color = pos_color, neg_color = neg_color,
      base_size = base_size
    ) + ggplot2::ggtitle("A) Selection Frequency")
  }
  
  # B) Weight plot (butterfly format)
  plots$weights <- .plot_weights_butterfly(
    object, var_info = var_info,
    pos_color = pos_color, neg_color = neg_color,
    base_size = base_size
  ) + ggplot2::ggtitle("B) Variable Weights")
  
  # C) Inference result plot
  if (!is.null(.get_active_refit_info(object)) || !is.null(object$boot_info)) {
    plots$validation <- plot_inference_results(
      object,
      conf_level = conf_level,
      include_covariates = include_covariates,
      pos_color = pos_color,
      neg_color = neg_color,
      show_significance = show_significance,
      base_size = base_size,
      title = NULL
    ) + ggplot2::ggtitle("C) Mixture Effect Estimates")
  }
  
  # Determine layout
  if (layout == "auto") {
    if (length(plots) == 3) {
      layout <- "vertical"  # A|B on top, C on bottom
    } else {
      layout <- "horizontal"
    }
  }
  
  # Combine if patchwork available
  if (requireNamespace("patchwork", quietly = TRUE)) {
    if (length(plots) >= 2) {
      if (layout == "vertical" && length(plots) == 3) {
        # A | B
        # ---C---
        top_row <- patchwork::wrap_plots(plots$selection, plots$weights, ncol = 2)
        combined <- top_row / plots$validation + 
          patchwork::plot_layout(heights = c(2, 1))
      } else {
        combined <- patchwork::wrap_plots(plots, ncol = 2)
      }
      return(combined)
    }
  }
  
  # Return list if patchwork unavailable
  if (length(plots) > 0) {
    message("Install 'patchwork' package for combined plot layout.")
    message("Returning list of individual plots.")
    return(plots)
  }
  
  return(invisible(NULL))
}


#' Compute Unified Variable Order
#' @keywords internal
.compute_unified_var_order <- function(object, sort_by = "combined", top_n = 10) {
  
  df <- data.frame(
    variable = object$var_names,
    pos_weight = as.numeric(object$pos_weights),
    neg_weight = as.numeric(object$neg_weights),
    stringsAsFactors = FALSE
  )
  
  # Add bootstrap information if available
  if (object$bootstrap && !is.null(object$boot_info)) {
    df$selection_freq_pos <- object$boot_info$selection_freq_pos[df$variable]
    df$selection_freq_neg <- object$boot_info$selection_freq_neg[df$variable]
    df$selection_freq <- pmax(df$selection_freq_pos, df$selection_freq_neg, na.rm = TRUE)
  } else {
    df$selection_freq <- 1
  }
  
  # Calculate sort criteria
  df$weight_max <- pmax(df$pos_weight, df$neg_weight)
  
  if (sort_by == "combined") {
    df$sort_val <- df$selection_freq * df$weight_max
  } else if (sort_by == "selection_freq") {
    df$sort_val <- df$selection_freq
  } else {
    df$sort_val <- df$weight_max
  }
  
  # If group information exists, sort within groups
  if (!is.null(object$groups)) {
    df$group <- .assign_group_labels(df$variable, object$groups)
    group_order <- .group_display_order(object$groups, unique(df$group))

    # Select top_n per group
    df_list <- lapply(group_order, function(grp) {
      g <- df[df$group == grp, ]
      g <- g[order(-g$sort_val), ]
      head(g, top_n)
    })
    df_top <- do.call(rbind, df_list)
  } else {
    df_top <- df[order(-df$sort_val), ]
    df_top <- head(df_top, top_n)
    df_top$group <- "All"
  }
  
  # Finalize variable order (group order, then sort_val descending)
  var_order <- df_top$variable
  
  # Calculate group boundaries (variable index positions)
  group_boundaries <- NULL
  if (!is.null(object$groups) && length(unique(df_top$group)) > 1) {
    # Get index of last variable in each group
    n_vars <- nrow(df_top)
    current_group <- df_top$group[1]
    boundaries <- c()
    for (i in 1:(n_vars - 1)) {
      if (df_top$group[i] != df_top$group[i + 1]) {
        # For coord_flip, boundary is between variables (i + 0.5)
        boundaries <- c(boundaries, i + 0.5)
      }
    }
    group_boundaries <- boundaries
  }
  
  # Return as list
  return(list(
    var_order = var_order,
    group_boundaries = group_boundaries,
    df = df_top
  ))
}


#' Plot Selection Frequency as Butterfly Chart
#' @keywords internal
.plot_selection_butterfly <- function(object, var_info = NULL,
                                       pos_color = "#3498db", neg_color = "#e74c3c",
                                       base_size = 11) {
  
  if (!object$bootstrap || is.null(object$boot_info)) {
    stop("Bootstrap information not available")
  }
  
  df <- data.frame(
    variable = names(object$boot_info$selection_freq_pos),
    pos_freq = as.numeric(object$boot_info$selection_freq_pos),
    neg_freq = as.numeric(object$boot_info$selection_freq_neg),
    stringsAsFactors = FALSE
  )
  
  # Get filter, sort, group information from var_info
  has_groups <- FALSE
  if (!is.null(var_info) && is.list(var_info) && !is.null(var_info$df)) {
    var_order <- var_info$var_order
    df <- df[df$variable %in% var_order, ]
    
    # Merge group information
    group_info <- var_info$df[, c("variable", "group")]
    df <- merge(df, group_info, by = "variable", all.x = TRUE)
    
    # Check if multiple groups exist
    if (length(unique(df$group)) > 1) {
      has_groups <- TRUE
      # Maintain group order
      group_order <- unique(var_info$df$group)
      df$group <- factor(df$group, levels = group_order)
    }
    
    # Sort within groups
    if (has_groups) {
      df_list <- lapply(split(df, df$group), function(g) {
        g_order <- var_order[var_order %in% g$variable]
        g$variable <- factor(g$variable, levels = rev(g_order))
        g
      })
      df <- do.call(rbind, df_list)
    } else {
      df$variable <- factor(df$variable, levels = rev(var_order))
    }
  } else if (!is.null(var_info)) {
    var_order <- if (is.list(var_info)) var_info$var_order else var_info
    df <- df[df$variable %in% var_order, ]
    df$variable <- factor(df$variable, levels = rev(var_order))
  }
  
  # Create plot
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = -.data$neg_freq),
      stat = "identity", fill = neg_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = .data$pos_freq),
      stat = "identity", fill = pos_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", 
                        color = "gray60", linewidth = 0.3) +
    ggplot2::scale_y_continuous(
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      labels = function(x) sprintf("%.0f%%", abs(x) * 100)
    ) +
    ggplot2::labs(x = NULL, y = "Negative                    Positive") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  
  # Add facet for multiple groups
  if (has_groups) {
    p <- p + ggplot2::facet_wrap(
      ~ group, 
      scales = "free_y", 
      ncol = 1,
      strip.position = "right"
    ) +
    ggplot2::theme(
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )
  }
  
  return(p)
}


#' Plot Weights as Butterfly Chart
#' @keywords internal
.plot_weights_butterfly <- function(object, var_info = NULL,
                                     pos_color = "#3498db", neg_color = "#e74c3c",
                                     base_size = 11) {
  
  df <- data.frame(
    variable = object$var_names,
    pos_weight = as.numeric(object$pos_weights),
    neg_weight = as.numeric(object$neg_weights),
    stringsAsFactors = FALSE
  )
  
  # Get filter, sort, group information from var_info
  has_groups <- FALSE
  if (!is.null(var_info) && is.list(var_info) && !is.null(var_info$df)) {
    var_order <- var_info$var_order
    df <- df[df$variable %in% var_order, ]
    
    # Merge group information
    group_info <- var_info$df[, c("variable", "group")]
    df <- merge(df, group_info, by = "variable", all.x = TRUE)
    
    # Check if multiple groups exist
    if (length(unique(df$group)) > 1) {
      has_groups <- TRUE
      # Maintain group order
      group_order <- unique(var_info$df$group)
      df$group <- factor(df$group, levels = group_order)
    }
    
    # Sort within groups
    if (has_groups) {
      df_list <- lapply(split(df, df$group), function(g) {
        g_order <- var_order[var_order %in% g$variable]
        g$variable <- factor(g$variable, levels = rev(g_order))
        g
      })
      df <- do.call(rbind, df_list)
    } else {
      df$variable <- factor(df$variable, levels = rev(var_order))
    }
  } else if (!is.null(var_info)) {
    var_order <- if (is.list(var_info)) var_info$var_order else var_info
    df <- df[df$variable %in% var_order, ]
    df$variable <- factor(df$variable, levels = rev(var_order))
  }
  
  # Determine axis range
  max_weight <- max(c(df$pos_weight, df$neg_weight), na.rm = TRUE)
  axis_limit <- ceiling(max_weight * 10) / 10
  axis_limit <- max(axis_limit, 0.1)
  
  # Create plot
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = -.data$neg_weight),
      stat = "identity", fill = neg_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = .data$pos_weight),
      stat = "identity", fill = pos_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5) +
    ggplot2::scale_y_continuous(
      limits = c(-axis_limit, axis_limit),
      labels = function(x) sprintf("%.2f", abs(x))
    ) +
    ggplot2::labs(x = NULL, y = "Negative                    Positive") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  
  # Add facet for multiple groups
  if (has_groups) {
    p <- p + ggplot2::facet_wrap(
      ~ group, 
      scales = "free_y", 
      ncol = 1,
      strip.position = "right"
    ) +
    ggplot2::theme(
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )
  }
  
  return(p)
}


#' Improved Forest Plot with Asterisk Significance
#' @keywords internal
.plot_forest_improved <- function(object,
                                   pos_color = "#3498db", neg_color = "#e74c3c",
                                   show_significance = "asterisk",
                                   conf_level = 0.95,
                                   base_size = 11,
                                   include_covariates = FALSE) {
  active_info <- .get_active_refit_info(object)
  if (is.null(active_info)) {
    stop("No inference results available. Fit with `validation = TRUE` or `refit = \"full\"` to enable inference.")
  }
  
  # Get results by group
  if (!is.null(active_info$group_results) &&
      length(active_info$group_results) > 0) {
    results <- active_info$group_results
    
    plot_data <- do.call(rbind, lapply(names(results), function(grp) {
      r <- results[[grp]]
      data.frame(
        group = grp,
        direction = c("Positive", "Negative"),
        estimate = c(r$pos_estimate, r$neg_estimate),
        se = c(r$pos_se, r$neg_se),
        p_value = c(r$pos_pvalue, r$neg_pvalue),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    # Only overall results
    vi <- active_info
    plot_data <- data.frame(
      group = "Overall",
      direction = c("Positive", "Negative"),
      estimate = c(vi$wqs_pos_estimate, vi$wqs_neg_estimate),
      se = c(vi$wqs_pos_se, vi$wqs_neg_se),
      p_value = c(vi$wqs_pos_pvalue, vi$wqs_neg_pvalue),
      stringsAsFactors = FALSE
    )
  }

  # Append covariate rows when requested
  if (isTRUE(include_covariates) && !is.null(object$cov_names)) {
    coef_table <- active_info$coef_table
    p_col <- active_info$p_col
    for (cov_nm in object$cov_names) {
      if (cov_nm %in% rownames(coef_table)) {
        cov_row <- data.frame(
          group = cov_nm,
          direction = "Covariate",
          estimate = coef_table[cov_nm, "Estimate"],
          se = coef_table[cov_nm, "Std. Error"],
          p_value = coef_table[cov_nm, p_col],
          stringsAsFactors = FALSE
        )
        plot_data <- rbind(plot_data, cov_row)
      }
    }
  }

  if (!is.data.frame(plot_data) || nrow(plot_data) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5, y = 0.5,
          label = "No validation effects available",
          size = 4
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Validation Inference")
    )
  }
  
  # Identify NA effects and generate caption
  na_effects <- plot_data[is.na(plot_data$estimate), , drop = FALSE]
  caption_text <- NULL
  
  if (nrow(na_effects) > 0) {
    # Group by direction
    pos_na_groups <- na_effects$group[na_effects$direction == "Positive"]
    neg_na_groups <- na_effects$group[na_effects$direction == "Negative"]
    
    caption_parts <- c()
    
    if (length(pos_na_groups) > 0) {
      if (length(pos_na_groups) == 1) {
        caption_parts <- c(caption_parts, 
          paste0("Positive effects not shown in ", pos_na_groups))
      } else {
        caption_parts <- c(caption_parts,
          paste0("Positive effects not shown in ", 
                 paste(pos_na_groups, collapse = " and ")))
      }
    }
    
    if (length(neg_na_groups) > 0) {
      if (length(neg_na_groups) == 1) {
        caption_parts <- c(caption_parts,
          paste0("Negative effects not shown in ", neg_na_groups))
      } else {
        caption_parts <- c(caption_parts,
          paste0("Negative effects not shown in ",
                 paste(neg_na_groups, collapse = " and ")))
      }
    }
    
    if (length(caption_parts) > 0) {
      caption_text <- paste0(
        paste(caption_parts, collapse = "; "),
        "\n(coefficients = 0 due to regularization)"
      )
    }
  }
  
  # Calculate CI (model-aware: t for gaussian/survey, z for binomial)
  refit_fit <- if (!is.null(object$refit_info$refit_fit)) {
    object$refit_info$refit_fit
  } else if (!is.null(active_info$glm_fit)) {
    active_info$glm_fit
  } else {
    NULL
  }
  z_val <- if (!is.null(refit_fit)) {
    .prediction_interval_multiplier(refit_fit, conf_level)
  } else {
    stats::qnorm(1 - (1 - conf_level) / 2)
  }
  plot_data$ci_lower <- plot_data$estimate - z_val * plot_data$se
  plot_data$ci_upper <- plot_data$estimate + z_val * plot_data$se
  
  # Create significance labels
  if (show_significance == "asterisk") {
    plot_data$sig_label <- sapply(plot_data$p_value, function(p) {
      if (is.na(p)) return("")
      if (p < 0.001) return("***")
      if (p < 0.01) return("**")
      if (p < 0.05) return("*")
      return("")
    })
  } else if (show_significance == "pvalue") {
    plot_data$sig_label <- sapply(plot_data$p_value, function(p) {
      if (is.na(p)) return("")
      format.pval(p, digits = 2)
    })
  } else {
    plot_data$sig_label <- ""
  }
  
  # Create labels (group name + direction)
  plot_data$label <- ifelse(
    plot_data$direction == "Covariate",
    plot_data$group,
    paste(plot_data$group, plot_data$direction, sep = "\n")
  )
  plot_data$label <- factor(plot_data$label, levels = rev(unique(plot_data$label)))

  # Set colors
  plot_data$color <- ifelse(
    plot_data$direction == "Positive", pos_color,
    ifelse(plot_data$direction == "Covariate", "gray40", neg_color)
  )
  
  # Plot data excluding NA (prevent warnings)
  plot_data_valid <- plot_data[!is.na(plot_data$estimate), , drop = FALSE]
  if (nrow(plot_data_valid) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5, y = 0.5,
          label = "No non-missing validation effects available",
          size = 4
        ) +
        ggplot2::theme_void() +
        ggplot2::labs(title = "Validation Inference", caption = caption_text)
    )
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data_valid, ggplot2::aes(x = .data$estimate, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper),
      orientation = "y",
      width = 0.2, linewidth = 0.8, color = "gray30",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(size = 3, color = plot_data_valid$color, na.rm = TRUE)
  
  # Position asterisks away from CI right end
  if (any(plot_data_valid$sig_label != "", na.rm = TRUE)) {
    # Calculate x-axis range and determine offset
    x_range <- max(plot_data_valid$ci_upper, na.rm = TRUE) - min(plot_data_valid$ci_lower, na.rm = TRUE)
    sig_offset <- x_range * 0.08  # 8% offset
    
    p <- p + ggplot2::geom_text(
      ggplot2::aes(x = .data$ci_upper + sig_offset, label = .data$sig_label),
      hjust = 0, vjust = 0.5, size = 4, fontface = "bold",
      na.rm = TRUE
    )
  }
  
  p <- p +
    ggplot2::labs(
      x = "Coefficient (95% CI)", 
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
