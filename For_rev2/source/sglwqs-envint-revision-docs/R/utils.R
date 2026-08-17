#' Null-coalescing helper
#'
#' @param a Primary value.
#' @param b Fallback value used when \code{a} is \code{NULL}.
#' @return \code{a} if non-\code{NULL}, otherwise \code{b}.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Validate and Complete Group Definitions
#'
#' Ensures group definitions are well-formed: named, disjoint, and reference
#' only known exposure variables. Optionally adds an \code{"Other"} group for
#' unassigned exposures so downstream code always operates on a complete
#' partition.
#'
#' @param var_names Character vector of exposure variable names.
#' @param groups Named list of character vectors, each listing exposure names
#'   belonging to that group.
#' @param add_other Logical. If \code{TRUE} (default), unassigned exposures
#'   are collected into an \code{"Other"} group.
#' @param other_name Character. Name of the auto-generated group for unassigned
#'   exposures (default: \code{"Other"}).
#'
#' @return A validated (and possibly augmented) named list of groups, or
#'   \code{NULL} when \code{groups} is \code{NULL}.
#'
#' @keywords internal
.validate_and_complete_groups <- function(var_names, groups,
                                          add_other = TRUE,
                                          other_name = "Other") {
  if (is.null(groups)) return(NULL)

  if (!is.list(groups) || is.null(names(groups)) || any(!nzchar(names(groups)))) {
    stop("`groups` must be a named list.", call. = FALSE)
  }
  if (anyDuplicated(names(groups))) {
    dup <- unique(names(groups)[duplicated(names(groups))])
    stop("Duplicated group names: ", paste(dup, collapse = ", "), call. = FALSE)
  }

  groups <- lapply(groups, function(x) unique(as.character(x)))

  flat <- unlist(groups, use.names = FALSE)
  unknown <- setdiff(flat, var_names)
  if (length(unknown) > 0) {
    stop("Unknown exposure(s) in `groups`: ",
         paste(unique(unknown), collapse = ", "), call. = FALSE)
  }

  dup_vars <- unique(flat[duplicated(flat)])
  if (length(dup_vars) > 0) {
    stop(
      "Each exposure may belong to at most one group. Overlapping exposure(s): ",
      paste(dup_vars, collapse = ", "),
      call. = FALSE
    )
  }

  assigned <- flat
  missing_vars <- setdiff(var_names, assigned)
  if (add_other && length(missing_vars) > 0) {
    if (other_name %in% names(groups)) {
      stop(
        "`groups` includes an explicit group named '", other_name,
        "', but some exposures are still unassigned: ",
        paste(missing_vars, collapse = ", "), ". ",
        "When '", other_name, "' is supplied explicitly, it must contain all intended members. ",
        "Either add these exposures to `groups[['", other_name, "']]` (or another group), ",
        "or remove the explicit '", other_name, "' group and let it be created automatically.",
           call. = FALSE)
    }
    groups[[other_name]] <- missing_vars
  }

  groups
}


#' Assign group labels to a variable column
#'
#' Given a character vector of variable names and a canonical groups list,
#' returns a character vector of group labels.  Variables not found in any
#' group are labelled \code{"Other"}.
#'
#' @param variables Character vector of variable names.
#' @param groups A named list mapping group names to variable name vectors,
#'   or \code{NULL}.
#' @return Character vector the same length as \code{variables}.
#' @keywords internal
.assign_group_labels <- function(variables, groups) {
  out <- rep(NA_character_, length(variables))
  if (!is.null(groups)) {
    for (grp_name in names(groups)) {
      out[variables %in% groups[[grp_name]]] <- grp_name
    }
  }
  out[is.na(out)] <- "Other"
  out
}


#' Canonical group ordering
#'
#' Returns the group display order, de-duplicating \code{"Other"} if it is
#' already present in the group names.  Only groups that actually appear in
#' \code{present} are kept.
#'
#' @param groups A named list of groups, or \code{NULL}.
#' @param present Character vector of group labels actually observed in data.
#' @return Character vector of ordered group labels.
#' @keywords internal
.group_display_order <- function(groups, present) {
  if (is.null(groups)) return("All")
  order <- unique(c(names(groups), "Other"))
  order[order %in% present]
}


#' Aggregate Bootstrap Weights (Mean of Per-Draw Weights)
#'
#' For each successful bootstrap draw, converts raw coefficients into
#' normalised weights using \code{calculate_weights()}, then averages those
#' weight vectors across draws.  This is the estimator described in the
#' manuscript (mean of weights, not weights-of-mean-coefficients).
#'
#' @param boot_pos_coef Bootstrap positive coefficient matrix (n_boot x p).
#' @param boot_neg_coef Bootstrap negative coefficient matrix (n_boot x p).
#' @param boot_success Logical vector of length n_boot.
#' @param var_names Character vector of variable names.
#' @param groups Optional named list of groups.
#' @param handle_collinearity Passed to \code{calculate_weights()}.
#'
#' @return A list with \code{pos_weights} and \code{neg_weights}.
#'
#' @keywords internal
.aggregate_bootstrap_weights <- function(boot_pos_coef, boot_neg_coef,
                                         boot_success, var_names,
                                         groups = NULL,
                                         handle_collinearity = "net") {
  if (is.null(boot_success)) {
    boot_success <- rowSums(abs(boot_pos_coef) + abs(boot_neg_coef)) > 0
  }
  idx <- which(boot_success)
  if (length(idx) == 0) {
    stop("No successful bootstrap iterations available for weight aggregation.")
  }

  p <- length(var_names)
  pos_accum <- stats::setNames(rep(0, p), var_names)
  neg_accum <- stats::setNames(rep(0, p), var_names)
  pos_index_sum_accum <- 0
  neg_index_sum_accum <- 0
  pos_index_sum_by_group_accum <- NULL
  neg_index_sum_by_group_accum <- NULL

  for (b in idx) {
    # Bootstrap aggregation should not emit one collinearity warning per draw.
    # The final fitted object still reports the high-level warning when present.
    w <- suppressWarnings(calculate_weights(
      pos_coef = boot_pos_coef[b, ],
      neg_coef = boot_neg_coef[b, ],
      var_names = var_names,
      groups = groups,
      handle_collinearity = handle_collinearity
    ))
    pos_accum <- pos_accum + w$pos_weights
    neg_accum <- neg_accum + w$neg_weights
    pos_index_sum_accum <- pos_index_sum_accum + w$pos_index_sum
    neg_index_sum_accum <- neg_index_sum_accum + w$neg_index_sum

    if (!is.null(w$pos_index_sum_by_group)) {
      if (is.null(pos_index_sum_by_group_accum)) {
        pos_index_sum_by_group_accum <- w$pos_index_sum_by_group
        neg_index_sum_by_group_accum <- w$neg_index_sum_by_group
      } else {
        for (nm in names(w$pos_index_sum_by_group)) {
          pos_index_sum_by_group_accum[[nm]] <- (if (is.null(pos_index_sum_by_group_accum[[nm]])) 0 else pos_index_sum_by_group_accum[[nm]]) + w$pos_index_sum_by_group[[nm]]
          neg_index_sum_by_group_accum[[nm]] <- (if (is.null(neg_index_sum_by_group_accum[[nm]])) 0 else neg_index_sum_by_group_accum[[nm]]) + w$neg_index_sum_by_group[[nm]]
        }
      }
    }
  }

  n <- length(idx)

  pos_index_sum_by_group_mean <- if (!is.null(pos_index_sum_by_group_accum)) {
    lapply(pos_index_sum_by_group_accum, function(x) x / n)
  }
  neg_index_sum_by_group_mean <- if (!is.null(neg_index_sum_by_group_accum)) {
    lapply(neg_index_sum_by_group_accum, function(x) x / n)
  }

  list(
    pos_weights = pos_accum / n,
    neg_weights = neg_accum / n,
    pos_index_sum = pos_index_sum_accum / n,
    neg_index_sum = neg_index_sum_accum / n,
    pos_index_sum_by_group = pos_index_sum_by_group_mean,
    neg_index_sum_by_group = neg_index_sum_by_group_mean
  )
}


#' Model-Aware CI Multiplier
#'
#' Returns the appropriate critical value for confidence intervals, using
#' the t-distribution when residual degrees of freedom are available
#' (survey GLMs, gaussian GLMs with estimated dispersion), and the
#' normal distribution otherwise (binomial, poisson).
#'
#' @param fit A fitted model object (glm, svyglm, or similar).
#' @param level Confidence level in (0, 1).
#'
#' @return A single numeric critical value.
#'
#' @keywords internal
.prediction_interval_multiplier <- function(fit, level) {
  alpha <- 1 - level

  if (inherits(fit, "svyglm")) {
    df <- tryCatch(stats::df.residual(fit), error = function(e) NA_real_)
    if (is.finite(df) && df > 0) {
      return(stats::qt(1 - alpha / 2, df = df))
    }
  }

  fam <- tryCatch(stats::family(fit)$family, error = function(e) NA_character_)
  if (!is.na(fam) && !fam %in% c("binomial", "poisson", "quasibinomial", "quasipoisson")) {
    df <- tryCatch(stats::df.residual(fit), error = function(e) NA_real_)
    if (is.finite(df) && df > 0) {
      return(stats::qt(1 - alpha / 2, df = df))
    }
  }

  stats::qnorm(1 - alpha / 2)
}


#' Sanitize Names for Safe Use in R Formulas
#'
#' Creates formula-safe versions of variable/group names while preserving
#' a mapping to the original names for display purposes. Uses \code{make.names()}
#' to ensure names are syntactically valid R identifiers.
#'
#' @param names Character vector of names to sanitize.
#'
#' @return A list containing:
#' \itemize{
#'   \item safe: Character vector of sanitized names (valid R identifiers)
#'   \item original: Character vector of original names
#'   \item to_original: Named character vector mapping safe -> original
#'   \item to_safe: Named character vector mapping original -> safe
#'   \item was_modified: Logical indicating whether any names were changed
#' }
#'
#' @details
#' This function is critical for supporting metabolomics-style variable names
#' that may contain parentheses, colons, slashes, and other special characters
#' (e.g., "PC(16:0/18:1)", "2-Aminobutyric acid", "alpha-Tocopherol").
#' Such names break R formula parsing if used directly.
#'
#' @examples
#' \dontrun{
#' result <- sanitize_for_formula(c("PC(16:0/18:1)", "LysoPC 18:0", "normal_name"))
#' result$safe       # c("PC.16.0.18.1.", "LysoPC.18.0", "normal_name")
#' result$to_original["PC.16.0.18.1."]  # "PC(16:0/18:1)"
#' }
#'
#' @keywords internal
sanitize_for_formula <- function(names) {
  if (is.null(names) || length(names) == 0) {
    return(list(
      safe = names,
      original = names,
      to_original = stats::setNames(character(0), character(0)),
      to_safe = stats::setNames(character(0), character(0)),
      was_modified = FALSE
    ))
  }
  
  safe <- make.names(names, unique = TRUE)
  was_modified <- !identical(safe, names)
  
  to_original <- stats::setNames(names, safe)
  to_safe <- stats::setNames(safe, names)
  
  list(
    safe = safe,
    original = names,
    to_original = to_original,
    to_safe = to_safe,
    was_modified = was_modified
  )
}


#' Validate Variable Names and Warn About Special Characters
#'
#' Checks exposure variable names and group names for characters that require
#' internal sanitization. Issues an informative message if names are modified.
#'
#' @param var_names Character vector of variable names.
#' @param groups Optional named list of groups.
#' @param verbose Logical. Whether to display messages.
#'
#' @return Invisible NULL. Called for side effects (messages/warnings).
#'
#' @keywords internal
validate_names <- function(var_names, groups = NULL, verbose = TRUE) {
  safe_vars <- make.names(var_names)
  modified_vars <- var_names[safe_vars != var_names]
  
  if (length(modified_vars) > 0 && verbose) {
    n_show <- min(5, length(modified_vars))
    examples <- paste0("  '", modified_vars[seq_len(n_show)], "'", collapse = "\n")
    if (length(modified_vars) > n_show) {
      examples <- paste0(examples, "\n  ... and ", length(modified_vars) - n_show, " more")
    }
    message(
      "Note: ", length(modified_vars), " variable name(s) contain special characters:\n",
      examples, "\n",
      "These are handled safely. Original names are preserved in all outputs."
    )
  }
  
  if (!is.null(groups)) {
    safe_grps <- make.names(names(groups))
    modified_grps <- names(groups)[safe_grps != names(groups)]
    
    if (length(modified_grps) > 0 && verbose) {
      message(
        "Note: Group name(s) with special characters: ",
        paste0("'", modified_grps, "'", collapse = ", "), "\n",
        "These are handled safely. Original names are preserved in all outputs."
      )
    }
  }
  
  invisible(NULL)
}


#' Process Covariates with Automatic Dummy Coding
#'
#' Converts covariates to a numeric matrix with proper handling of factors.
#' Column names are sanitized to ensure they are valid R identifiers for
#' safe use in formula construction (e.g., in validation GLM).
#'
#' @param covariates A data frame or matrix of covariates.
#'
#' @return A list containing:
#' \itemize{
#'   \item matrix: Numeric matrix of processed covariates (with sanitized colnames)
#'   \item names: Column names of the processed matrix (sanitized)
#'   \item original_names: Original covariate names (before sanitization)
#' }
#'
#' @keywords internal
process_covariates <- function(covariates) {
  if (is.null(covariates)) {
    return(list(matrix = NULL, names = NULL, original_names = NULL))
  }
  
  # Convert to data frame
  if (!is.data.frame(covariates)) {
    covariates <- as.data.frame(covariates)
  }
  
  original_names <- names(covariates)
  
  # Check for factor/character columns
  has_factor <- any(sapply(covariates, function(x) is.factor(x) || is.character(x)))
  
  if (has_factor) {
    # Dummy coding via model.matrix
    # Wrap column names in backticks for safe formula parsing
    # (original_names may contain spaces/special characters)
    safe_terms <- paste0("`", original_names, "`")
    formula_str <- paste("~ -1 +", paste(safe_terms, collapse = " + "))
    cov_matrix <- model.matrix(as.formula(formula_str), data = covariates)
    # Sanitize output column names for downstream formula safety
    # (factor levels may introduce additional special chars)
    colnames(cov_matrix) <- make.names(colnames(cov_matrix), unique = TRUE)
    cov_names <- colnames(cov_matrix)

    message("Note: Factor/character covariates were automatically dummy-coded.")
  } else {
    cov_matrix <- as.matrix(covariates)
    if (is.null(colnames(cov_matrix))) {
      cov_names <- paste0("Z", seq_len(ncol(cov_matrix)))
    } else {
      # Sanitize column names for downstream formula safety in validation GLM
      cov_names <- make.names(colnames(cov_matrix), unique = TRUE)
      colnames(cov_matrix) <- cov_names
    }
  }
  
  return(list(
    matrix = cov_matrix,
    names = cov_names,
    original_names = original_names
  ))
}


#' Check for Collinearity Warning
#'
#' Checks if any variable has both positive and negative coefficients non-zero.
#'
#' @param pos_coef Positive coefficients vector.
#' @param neg_coef Negative coefficients vector.
#' @param var_names Variable names.
#' @param tol Tolerance for considering a coefficient as non-zero.
#'
#' @return A list with warning information.
#'
#' @keywords internal
check_collinearity <- function(pos_coef, neg_coef, var_names, tol = 1e-6) {
  both_nonzero <- (abs(pos_coef) > tol) & (abs(neg_coef) > tol)
  
  if (any(both_nonzero)) {
    problem_vars <- var_names[both_nonzero]
    warning(
      "The following variables have both positive and negative coefficients non-zero:\n",
      paste("  -", problem_vars, collapse = "\n"), "\n",
      "This may indicate numerical instability. Net coefficients will be used for weights."
    )
    
    return(list(
      has_warning = TRUE,
      problem_vars = problem_vars,
      indices = which(both_nonzero)
    ))
  }
  
  return(list(has_warning = FALSE, problem_vars = NULL, indices = NULL))
}


#' Calculate Weights from Coefficients
#'
#' Calculates normalized weights from coefficients, with optional group-wise normalization.
#'
#' @param pos_coef Positive coefficients.
#' @param neg_coef Negative coefficients.
#' @param var_names Variable names.
#' @param groups Optional group list for group-wise normalization.
#' @param handle_collinearity How to handle variables with both pos and neg coefficients.
#'   "net" uses net coefficient, "warning" just warns.
#'
#' @return A list with weights and index sums.
#'
#' @keywords internal
calculate_weights <- function(pos_coef, neg_coef, var_names, groups = NULL,
                               handle_collinearity = "net") {
  
  p <- length(pos_coef)
  
  # Collinearity check
  collin_check <- check_collinearity(pos_coef, neg_coef, var_names)
  
  # Handle collinearity if detected
  if (collin_check$has_warning && handle_collinearity == "net") {
    for (idx in collin_check$indices) {
      net <- pos_coef[idx] - neg_coef[idx]
      if (net >= 0) {
        pos_coef[idx] <- net
        neg_coef[idx] <- 0
      } else {
        pos_coef[idx] <- 0
        neg_coef[idx] <- -net
      }
    }
  }
  
  # Compute weights by group
  if (!is.null(groups)) {
    pos_weights <- rep(0, p)
    neg_weights <- rep(0, p)
    names(pos_weights) <- var_names
    names(neg_weights) <- var_names
    
    pos_index_sum_by_group <- list()
    neg_index_sum_by_group <- list()
    
    for (grp_name in names(groups)) {
      grp_vars <- groups[[grp_name]]
      grp_idx <- which(var_names %in% grp_vars)
      
      grp_pos_sum <- sum(pos_coef[grp_idx])
      grp_neg_sum <- sum(neg_coef[grp_idx])
      
      pos_index_sum_by_group[[grp_name]] <- grp_pos_sum
      neg_index_sum_by_group[[grp_name]] <- grp_neg_sum
      
      if (grp_pos_sum > 0) {
        pos_weights[grp_idx] <- pos_coef[grp_idx] / grp_pos_sum
      }
      if (grp_neg_sum > 0) {
        neg_weights[grp_idx] <- neg_coef[grp_idx] / grp_neg_sum
      }
    }
    
    # Variables not belonging to any group
    all_grp_vars <- unlist(groups)
    other_idx <- which(!var_names %in% all_grp_vars)
    if (length(other_idx) > 0) {
      other_pos_sum <- sum(pos_coef[other_idx])
      other_neg_sum <- sum(neg_coef[other_idx])
      pos_index_sum_by_group[["Other"]] <- other_pos_sum
      neg_index_sum_by_group[["Other"]] <- other_neg_sum
      
      if (other_pos_sum > 0) {
        pos_weights[other_idx] <- pos_coef[other_idx] / other_pos_sum
      }
      if (other_neg_sum > 0) {
        neg_weights[other_idx] <- neg_coef[other_idx] / other_neg_sum
      }
    }
    
    pos_sum <- sum(pos_coef)
    neg_sum <- sum(neg_coef)
    
  } else {
    # Compute weights across all variables
    pos_sum <- sum(pos_coef)
    neg_sum <- sum(neg_coef)
    
    if (pos_sum > 0) {
      pos_weights <- pos_coef / pos_sum
    } else {
      pos_weights <- rep(0, p)
    }
    names(pos_weights) <- var_names
    
    if (neg_sum > 0) {
      neg_weights <- neg_coef / neg_sum
    } else {
      neg_weights <- rep(0, p)
    }
    names(neg_weights) <- var_names
    
    pos_index_sum_by_group <- NULL
    neg_index_sum_by_group <- NULL
  }
  
  return(list(
    pos_weights = pos_weights,
    neg_weights = neg_weights,
    pos_coef = pos_coef,
    neg_coef = neg_coef,
    pos_index_sum = pos_sum,
    neg_index_sum = neg_sum,
    pos_index_sum_by_group = pos_index_sum_by_group,
    neg_index_sum_by_group = neg_index_sum_by_group,
    collinearity_warning = collin_check$has_warning
  ))
}


#' Create WQS Index from Weights and Data
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param pos_weights Positive weights.
#' @param neg_weights Negative weights.
#'
#' @return A data frame with positive and negative WQS indices.
#'
#' @keywords internal
create_wqs_index <- function(X_quantile, pos_weights, neg_weights) {
  pos_index <- as.numeric(X_quantile %*% pos_weights)
  neg_index <- as.numeric(X_quantile %*% neg_weights)
  
  data.frame(
    wqs_pos = pos_index,
    wqs_neg = neg_index
  )
}
