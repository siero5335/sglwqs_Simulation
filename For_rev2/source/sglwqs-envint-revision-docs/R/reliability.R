#' Compute diagnostic summaries for an sglwqs fit
#'
#' Returns a list of observed model and data characteristics. Each item reports
#' an observed value and, where applicable, a commonly cited reference. This
#' function does not label fits as reliable or unreliable.
#'
#' @param fit A fitted \code{sglwqs} or \code{sglwqs_mids} object.
#' @return A list of class \code{"sglwqs_diagnostics"} or
#'   \code{"sglwqs_mids_diagnostics"}.
#' @export
compute_diagnostics <- function(fit) {
  UseMethod("compute_diagnostics")
}

#' @export
compute_diagnostics.sglwqs <- function(fit) {
  items <- list()

  n <- fit$n %||% NA_integer_
  p_exp <- fit$p %||% length(fit$var_names %||% character(0))
  p_cov <- length(fit$cov_names %||% character(0))
  p_total <- p_exp + p_cov

  if (is.finite(n) && p_total > 0) {
    items$sample_to_predictor_ratio <- .mk_diagnostic_item(
      value = n / p_total,
      text = sprintf(
        "n = %d, predictors (exposures + covariates) = %d, ratio = %.2f",
        n, p_total, n / p_total
      ),
      reference = "Harrell (2015) Section 4.4 discusses n/p guidelines."
    )
  }

  if (identical(fit$family, "binomial") && !is.null(fit$y) && p_total > 0) {
    event_counts <- table(fit$y)
    n_events <- if (length(event_counts) >= 2) min(event_counts) else 0L
    items$events_to_predictor_ratio <- .mk_diagnostic_item(
      value = n_events / p_total,
      text = sprintf(
        "Minority class count = %d, predictors = %d, EPV = %.2f",
        n_events, p_total, n_events / p_total
      ),
      reference = "Peduzzi et al. (1996) suggests EPV >= 10."
    )
  }

  n_nz_pos <- sum(fit$pos_weights > 0, na.rm = TRUE)
  n_nz_neg <- sum(fit$neg_weights > 0, na.rm = TRUE)
  n_nz_union <- sum((fit$pos_weights > 0) | (fit$neg_weights > 0), na.rm = TRUE)
  items$weights_nonzero_count <- .mk_diagnostic_item(
    value = n_nz_union,
    text = sprintf(
      "Non-zero weights: positive = %d, negative = %d, union = %d / %d exposures",
      n_nz_pos, n_nz_neg, n_nz_union, p_exp
    )
  )

  if (is.finite(p_exp) && p_exp > 0) {
    items$weights_nonzero_fraction <- .mk_diagnostic_item(
      value = n_nz_union / p_exp,
      text = sprintf(
        "Non-zero weight fraction: %d / %d exposures (%.1f%%)",
        n_nz_union, p_exp, 100 * n_nz_union / p_exp
      )
    )
  }

  if (n_nz_union == 0) {
    items$weights_all_zero <- .mk_diagnostic_item(
      value = 0,
      text = "All positive and negative weights are zero (SGL shrank all coefficients)."
    )
  }

  sel_diag <- fit$selection_diagnostics
  if (!is.null(sel_diag)) {
    items$selection_lambda_path <- .mk_diagnostic_item(
      value = sel_diag$lambda_path_length %||% NA_integer_,
      text = sprintf(
        "Selection lambda path: length = %s, range = [%s, %s], selected = %s%s",
        sel_diag$lambda_path_length %||% NA_integer_,
        signif(sel_diag$lambda_path_min %||% NA_real_, 4),
        signif(sel_diag$lambda_path_max %||% NA_real_, 4),
        signif(sel_diag$selected_lambda %||% NA_real_, 4),
        if (isTRUE(sel_diag$selected_lambda_at_path_boundary)) " (path boundary)" else ""
      )
    )
    items$selection_exposure_nonzero <- .mk_diagnostic_item(
      value = (sel_diag$nonzero_positive_coef %||% 0L) +
        (sel_diag$nonzero_negative_coef %||% 0L),
      text = sprintf(
        "Selection non-zero exposure coefficients: positive = %d, negative = %d; weight sums: positive = %.3f, negative = %.3f",
        sel_diag$nonzero_positive_coef %||% 0L,
        sel_diag$nonzero_negative_coef %||% 0L,
        sel_diag$positive_weight_sum %||% NA_real_,
        sel_diag$negative_weight_sum %||% NA_real_
      )
    )
    if (isTRUE(sel_diag$all_zero_exposure)) {
      items$selection_all_zero_exposure <- .mk_diagnostic_item(
        value = TRUE,
        text = "Selection returned all-zero exposure coefficients; downstream WQS exposure indices may be absent."
      )
    }
  }

  overlap_n <- sum((fit$pos_coef > 0) & (fit$neg_coef > 0), na.rm = TRUE)
  if (overlap_n > 0) {
    items$exposure_directional_overlap <- .mk_diagnostic_item(
      value = overlap_n,
      text = sprintf("%d exposures have non-zero coefficients in both directions.", overlap_n)
    )
  }

  if (!is.null(fit$X_quantile) && !is.null(fit$y) && ncol(fit$X_quantile) > 0) {
    cor_vals <- tryCatch(
      apply(fit$X_quantile, 2, function(x) stats::cor(x, fit$y, use = "complete.obs")),
      error = function(e) NULL
    )
    if (!is.null(cor_vals)) {
      max_abs_r <- max(abs(cor_vals), na.rm = TRUE)
      if (is.finite(max_abs_r)) {
        items$max_abs_exposure_outcome_r <- .mk_diagnostic_item(
          value = max_abs_r,
          text = sprintf("Maximum absolute exposure-outcome correlation: %.3f", max_abs_r)
        )
      }
    }
  }

  if (!is.null(fit$cov_matrix) && !is.null(fit$y) && identical(fit$family, "gaussian")) {
    cov_lm <- tryCatch(
      suppressWarnings(stats::lm(fit$y ~ fit$cov_matrix)),
      error = function(e) NULL
    )
    if (!is.null(cov_lm)) {
      cov_lm_s <- summary(cov_lm)
      cov_r2 <- cov_lm_s$r.squared
      f_p <- if (!is.null(cov_lm_s$fstatistic)) {
        stats::pf(
          cov_lm_s$fstatistic[1],
          cov_lm_s$fstatistic[2],
          cov_lm_s$fstatistic[3],
          lower.tail = FALSE
        )
      } else {
        NA_real_
      }
      if (is.finite(cov_r2)) {
        items$covariates_marginal_r2 <- .mk_diagnostic_item(
          value = cov_r2,
          text = sprintf(
            "Covariate-only linear model: R^2 = %.3f, F-test p = %.3f",
            cov_r2, f_p
          )
        )
      }
    }
  }

  if (!is.null(fit$fit$cvm)) {
    cvm <- fit$fit$cvm
    rel_range <- (max(cvm) - min(cvm)) / max(cvm)
    items$cv_curve_relative_range <- .mk_diagnostic_item(
      value = rel_range,
      text = sprintf(
        "sparsegl CV curve: range = [%.4f, %.4f], relative range = %.3f",
        min(cvm), max(cvm), rel_range
      )
    )
  }

  if (!is.null(fit$boot_info)) {
    n_total <- length(fit$boot_info$boot_success %||% logical(0))
    n_success <- fit$boot_info$n_successful %||% NA_integer_
    if (n_total > 0 && is.finite(n_success)) {
      items$bootstrap_success_rate <- .mk_diagnostic_item(
        value = n_success / n_total,
        text = sprintf(
          "Bootstrap iterations succeeded: %d / %d (%.1f%%)",
          n_success, n_total, 100 * n_success / n_total
        )
      )
    }
  }

  active_info <- .get_active_refit_info(fit)
  if (!is.null(active_info) && identical(attr(active_info, "source"), "validation_info")) {
    n_val <- active_info$n_val %||% NA_integer_
    n_params <- NA_integer_
    fit_obj <- active_info$refit_fit %||% active_info$glm_fit
    if (!is.null(fit_obj)) {
      n_params <- tryCatch(length(stats::coef(fit_obj)), error = function(e) NA_integer_)
    }
    if (is.finite(n_val)) {
      items$validation_set_size <- .mk_diagnostic_item(
        value = n_val,
        text = if (is.finite(n_params)) {
          sprintf("Validation set size = %d with %d coefficients in the downstream model.", n_val, n_params)
        } else {
          sprintf("Validation set size = %d.", n_val)
        }
      )
    }
  }
  fit_obj_active <- NULL
  if (!is.null(active_info)) {
    fit_obj_active <- active_info$refit_fit %||% active_info$glm_fit
  }
  if (!is.null(fit_obj_active)) {
    aliased <- tryCatch(
      summary(fit_obj_active)$aliased,
      error = function(e) NULL
    )
    if (!is.null(aliased) && any(aliased)) {
      items$rank_deficient_coefficients <- .mk_diagnostic_item(
        value = sum(aliased),
        text = sprintf(
          "Refit GLM: %d of %d coefficients not defined due to rank deficiency. Aliased terms: %s",
          sum(aliased), length(aliased),
          paste(names(aliased)[aliased], collapse = ", ")
        )
      )
    }
  }

  has_inference <- !is.null(active_info)
  has_boot <- !is.null(fit$boot_info)
  if (!has_inference && !has_boot) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = "point_estimates_only",
      text = "Configured path: point estimates only."
    )
  } else if (has_boot && !has_inference) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = "bootstrap_only",
      text = "Configured path: bootstrap-only summaries (no downstream GLM refit)."
    )
  } else if (has_inference && identical(attr(active_info, "source"), "validation_info")) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = "validation_info",
      text = "Configured path: training/validation downstream GLM inference."
    )
  } else if (has_inference) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = if (has_boot) "bootstrap_plus_refit_full" else "refit_full",
      text = if (has_boot) {
        "Configured path: bootstrap-stabilized weights with full-data downstream GLM inference."
      } else {
        "Configured path: full-data downstream GLM inference."
      }
    )
  }

  if (!has_inference && !has_boot) {
    items$inference_mode <- .mk_diagnostic_item(
      value = "none",
      text = "No inference mode enabled (bootstrap = FALSE, refit = 'none'). Only point estimates available."
    )
  } else if (has_boot && !has_inference) {
    items$inference_mode <- .mk_diagnostic_item(
      value = "bootstrap_only",
      text = "Inference source: bootstrap summaries only (no downstream GLM refit available)."
    )
  } else if (has_inference) {
    items$inference_mode <- .mk_diagnostic_item(
      value = attr(active_info, "source") %||% "unknown",
      text = sprintf("Inference source: %s.", attr(active_info, "source") %||% "unknown")
    )
  }

  structure(items, class = c("sglwqs_diagnostics", "list"))
}

#' @export
compute_diagnostics.sglwqs_mids <- function(fit) {
  items <- list()

  items$imputation_count <- .mk_diagnostic_item(
    value = fit$m,
    text = sprintf("Number of imputations: m = %d", fit$m),
    reference = "White et al. (2011) suggests m >= ceiling(FMI * 100)."
  )

  if (fit$n_successful < fit$m) {
    items$imputation_success_rate <- .mk_diagnostic_item(
      value = fit$n_successful / fit$m,
      text = sprintf(
        "Imputations that fit successfully: %d / %d",
        fit$n_successful, fit$m
      )
    )
  }

  inf <- fit$pooled$inference %||% fit$pooled$validation
  has_boot_inf <- !is.null(fit$pooled$bootstrap_inference)
  if (is.null(inf) && !has_boot_inf) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = "point_estimates_only",
      text = "Configured path: pooled point estimates only."
    )
  } else if (is.null(inf) && has_boot_inf) {
    items$two_stage_path <- .mk_diagnostic_item(
      value = "bootstrap_only",
      text = "Configured path: Rubin-pooled bootstrap-only summaries."
    )
  } else if (!is.null(inf)) {
    source_used <- inf$source_used %||% "unknown"
    items$two_stage_path <- .mk_diagnostic_item(
      value = source_used,
      text = if (identical(source_used, "validation_info")) {
        "Configured path: Rubin-pooled validation-split downstream GLM inference."
      } else if (identical(source_used, "refit_info")) {
        "Configured path: Rubin-pooled full-data downstream GLM inference."
      } else {
        sprintf("Configured path: Rubin-pooled downstream inference from %s.", source_used)
      }
    )
    items$pooled_inference_source <- .mk_diagnostic_item(
      value = source_used,
      text = sprintf("Pooled downstream inference source: %s.", source_used)
    )
  }

  if (!is.null(inf)) {
    for (dir_name in c("wqs_pos", "wqs_neg")) {
      fmi <- inf[[dir_name]]$fmi
      if (!is.null(fmi) && !is.na(fmi)) {
        items[[paste0("fmi_", dir_name)]] <- .mk_diagnostic_item(
          value = fmi,
          text = sprintf("Fraction of missing information (%s): %.2f", dir_name, fmi)
        )
      }
    }
  }

  if (has_boot_inf) {
    items$pooled_bootstrap_inference <- .mk_diagnostic_item(
      value = TRUE,
      text = "Pooled bootstrap-derived inference is available."
    )
  }

  per_imp <- lapply(Filter(Negate(is.null), fit$fits), compute_diagnostics)
  n_all_zero <- sum(vapply(per_imp, function(di) {
    wz <- di$weights_nonzero_count
    !is.null(wz) && identical(wz$value, 0)
  }, logical(1)))
  if (n_all_zero > 0) {
    items$imputations_with_all_zero_weights <- .mk_diagnostic_item(
      value = n_all_zero,
      text = sprintf(
        "Imputations with all weights zero: %d / %d",
        n_all_zero, length(per_imp)
      )
    )
  }

  items$per_imputation <- per_imp
  structure(items, class = c("sglwqs_mids_diagnostics", "list"))
}

#' @export
print.sglwqs_diagnostics <- function(x, ...) {
  items <- x[names(x) != "per_imputation"]

  if (length(items) == 0) {
    cat("Diagnostics: no items reported.\n")
    return(invisible(x))
  }

  cat("Diagnostics:\n")
  for (nm in names(items)) {
    cat(sprintf("  %s\n", items[[nm]]$text))
    if (!is.null(items[[nm]]$reference)) {
      cat(sprintf("    Reference: %s\n", items[[nm]]$reference))
    }
  }
  invisible(x)
}

#' @export
print.sglwqs_mids_diagnostics <- print.sglwqs_diagnostics

#' Emit diagnostics as fit-time messages
#' @keywords internal
.emit_diagnostics <- function(diagnostics, verbose = TRUE) {
  if (!isTRUE(verbose) || is.null(diagnostics)) {
    return(invisible(NULL))
  }

  for (nm in names(diagnostics)) {
    if (identical(nm, "per_imputation")) next
    item <- diagnostics[[nm]]
    if (identical(nm, "rank_deficient_coefficients")) {
      warning(sprintf("[%s] %s", nm, item$text), call. = FALSE, immediate. = TRUE)
    } else {
      message(sprintf("[%s] %s", nm, item$text))
    }
  }
  invisible(NULL)
}

#' Settings advice printed at fit start
#' @keywords internal
.advise_on_settings <- function(bootstrap, validation, refit, is_mids, verbose = TRUE) {
  if (!isTRUE(verbose)) {
    return(invisible(NULL))
  }

  if (!bootstrap && !validation && identical(refit, "none")) {
    message(
      "[settings] bootstrap = FALSE, validation = FALSE, refit = 'none'. ",
      "No standard errors or confidence intervals will be computed."
    )
  }

  if (isTRUE(is_mids) && !bootstrap && identical(refit, "none")) {
    message(
      "[settings] Multiple imputation without bootstrap or refit. ",
      "Rubin's rules require a within-imputation variance source; pooled outputs contain point estimates only."
    )
  }

  if (bootstrap && identical(refit, "full")) {
    message(
      "[settings] bootstrap = TRUE, refit = 'full'. WQS composite indices will be estimated from bootstrap-averaged weights and refit on full data."
    )
  }

  invisible(NULL)
}

#' @keywords internal
.mk_diagnostic_item <- function(value, text, reference = NULL) {
  list(
    value = value,
    text = text,
    reference = reference
  )
}
