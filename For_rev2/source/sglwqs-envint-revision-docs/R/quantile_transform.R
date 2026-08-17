#' Quantile Transformation of Data
#'
#' Transforms continuous variables into quantile-based categorical variables.
#' Handles NA values and variables with many tied values gracefully.
#'
#' @param data A numeric data frame or numeric matrix containing the variables
#'   to transform.
#' @param n_quantiles Integer. Number of quantiles to use for transformation (default: 4 for quartiles).
#' @param var_names Optional character vector of variable names. If NULL, uses column names of data.
#' @param breaks_list Optional list of pre-computed breaks for each variable.
#'   Used when applying training-set breaks to validation data.
#' @param weights Optional numeric vector of observation weights used to compute
#'   weighted quantile breaks.
#' @param method Character. Quantile break method. Defaults to
#'   \code{"empirical"} when \code{weights = NULL} and \code{"weighted_ecdf"}
#'   otherwise.
#'
#' @return A list containing:
#' \itemize{
#'   \item q: Matrix with quantile-transformed values (0 to n_quantiles-1)
#'   \item breaks: List of breaks used for each variable (for applying to new data)
#' }
#'
#' @details
#' The function handles several edge cases:
#' \itemize{
#'   \item NA values in input are preserved as NA in output
#'   \item Variables with many tied values may have fewer effective quantile bins
#'   \item When breaks have duplicates, the function automatically adjusts
#' }
#'
#' @examples
#' \dontrun{
#' data("metals", package = "qgcomp")
#' Xnm <- c("arsenic", "barium", "cadmium", "calcium", "chromium")
#' result <- quantile_transform(metals[, Xnm], n_quantiles = 4)
#' x_quantile <- result$q
#' 
#' # Apply same breaks to new data
#' result_new <- quantile_transform(new_data, breaks_list = result$breaks)
#' }
#'
#' @export
quantile_transform <- function(data, n_quantiles = 4, var_names = NULL,
                               breaks_list = NULL, weights = NULL,
                               method = c("empirical", "weighted_ecdf")) {
  if (missing(method)) {
    method <- if (is.null(weights)) "empirical" else "weighted_ecdf"
  } else {
    method <- match.arg(method)
  }
  
  data_colnames <- colnames(data)
  if (is.data.frame(data)) {
    is_numeric_col <- vapply(data, is.numeric, logical(1))
    if (!all(is_numeric_col)) {
      bad_cols <- names(data)[!is_numeric_col]
      stop(
        "`data` must contain only numeric columns for quantile transformation. ",
        "Non-numeric columns: ", paste(bad_cols, collapse = ", "), ".",
        call. = FALSE
      )
    }
    data <- as.matrix(data)
  } else {
    data <- as.matrix(data)
    if (!is.numeric(data)) {
      stop("`data` must be a numeric matrix or a data frame with only numeric columns.", call. = FALSE)
    }
  }
  
  n <- nrow(data)
  p <- ncol(data)
  
  if (!is.null(weights)) {
    weights <- .validate_quantile_weights(weights, n)
  }
  
  # Matrix to store results
  q <- matrix(NA_real_, n, p)
  
  # Store breaks (for validation)
  if (is.null(breaks_list)) {
    breaks_list <- vector("list", p)
    compute_breaks <- TRUE
  } else {
    compute_breaks <- FALSE
  }
  
  # Perform quantile transformation for each variable
  for (i in seq_len(p)) {
    x <- data[, i]
    not_na <- !is.na(x)
    
    if (sum(not_na) == 0) {
      # All NA case
      q[, i] <- NA_real_
      if (compute_breaks) breaks_list[[i]] <- NA
      next
    }
    
    if (compute_breaks) {
      if (identical(method, "weighted_ecdf") && !is.null(weights)) {
        break_info <- .compute_weighted_breaks(
          x = x[not_na],
          weights = weights[not_na],
          n_quantiles = n_quantiles
        )
        breaks <- break_info$breaks
      } else {
        breaks <- stats::quantile(
          x[not_na],
          probs = seq(0, 1, length.out = n_quantiles + 1),
          na.rm = TRUE,
          type = 7
        )
        break_info <- list(
          breaks = breaks,
          fallback = if (length(unique(breaks)) < length(breaks)) "rank" else "cut"
        )
      }
      
      # Handle duplicate breaks
      breaks_unique <- unique(breaks)
      
      if (length(breaks_unique) < 2) {
        # All identical values: assign 0
        q[not_na, i] <- 0
        breaks_list[[i]] <- list(breaks = breaks_unique, fallback = "constant")
        next
      }
      
      if (length(breaks_unique) < length(breaks)) {
        var_label <- if (!is.null(var_names) && length(var_names) >= i) {
          var_names[[i]]
        } else if (!is.null(data_colnames) && length(data_colnames) >= i) {
          data_colnames[[i]]
        } else {
          paste0("column ", i)
        }
        warning("Variable '", var_label, "' has many tied values. Using conservative fallback binning.")
        q[not_na, i] <- .assign_fallback_bins(
          x = x[not_na],
          n_quantiles = n_quantiles,
          weights = if (!is.null(weights)) weights[not_na] else NULL,
          method = break_info$fallback
        )
        breaks_list[[i]] <- list(
          breaks = breaks_unique,
          fallback = break_info$fallback
        )
        next
      }
      
      breaks_list[[i]] <- list(breaks = breaks, fallback = "cut")
    } else {
      break_spec <- .normalize_break_spec(breaks_list[[i]])
      breaks <- break_spec$breaks
      
      if (length(breaks) < 2 || all(is.na(breaks))) {
        q[not_na, i] <- 0
        next
      }
    }
    
    # Apply cut() (clamp values outside breaks range)
    x_to_cut <- x[not_na]
    if (!compute_breaks) {
      # Validation data: extend breaks edges to assign
      # out-of-range values to the nearest bin
      breaks[1] <- min(breaks[1], min(x_to_cut, na.rm = TRUE))
      breaks[length(breaks)] <- max(breaks[length(breaks)], max(x_to_cut, na.rm = TRUE))
    }
    cut_result <- cut(x_to_cut, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    q[not_na, i] <- cut_result - 1  # Convert to 0-based
  }
  
  # Check for NaN and Inf
  if (any(!is.finite(q) & !is.na(q))) {
    warning("Non-finite values detected in quantile transformation. Setting to NA.")
    q[!is.finite(q)] <- NA
  }
  
  # Set column names
  if (!is.null(var_names)) {
    colnames(q) <- var_names
  } else if (!is.null(colnames(data))) {
    colnames(q) <- colnames(data)
  } else {
    colnames(q) <- paste0("X", seq_len(p))
  }
  
  # Set names for breaks as well
  names(breaks_list) <- colnames(q)
  
  return(list(q = q, breaks = breaks_list))
}


#' @keywords internal
.validate_quantile_weights <- function(weights, n) {
  if (!is.numeric(weights)) {
    stop("`weights` must be a numeric vector.", call. = FALSE)
  }
  if (length(weights) != n) {
    stop("`weights` must have length equal to the number of rows in `data`.", call. = FALSE)
  }
  if (any(!is.finite(weights))) {
    stop("`weights` must contain only finite values.", call. = FALSE)
  }
  if (any(weights < 0)) {
    stop("`weights` must be non-negative.", call. = FALSE)
  }
  if (sum(weights) <= 0) {
    stop("`weights` must have positive total weight.", call. = FALSE)
  }
  weights
}


#' @keywords internal
.compute_weighted_breaks <- function(x, weights, n_quantiles) {
  valid <- !is.na(x) & !is.na(weights) & weights > 0
  x <- x[valid]
  weights <- weights[valid]
  
  if (length(x) == 0 || sum(weights) <= 0) {
    return(list(breaks = NA_real_, fallback = "constant"))
  }
  
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- weights[ord]
  probs <- seq(0, 1, length.out = n_quantiles + 1)
  cdf <- cumsum(w_ord) / sum(w_ord)
  idx <- vapply(
    probs,
    function(prob) which(cdf >= prob)[1],
    integer(1)
  )
  idx[is.na(idx)] <- length(x_ord)
  
  list(
    breaks = x_ord[idx],
    fallback = "weighted_rank"
  )
}


#' @keywords internal
.assign_fallback_bins <- function(x, n_quantiles, weights = NULL,
                                  method = c("rank", "weighted_rank")) {
  method <- match.arg(method)
  
  if (identical(method, "weighted_rank") && !is.null(weights)) {
    valid <- !is.na(x) & !is.na(weights)
    x_valid <- x[valid]
    w_valid <- weights[valid]
    probs <- seq(0, 1, length.out = n_quantiles + 1)
    
    level_weights <- tapply(w_valid, x_valid, sum)
    level_weights <- level_weights[order(as.numeric(names(level_weights)))]
    cdf_end <- cumsum(level_weights) / sum(level_weights)
    bins <- cut(cdf_end, breaks = probs, include.lowest = TRUE, labels = FALSE) - 1
    bins[is.na(bins)] <- 0
    
    level_map <- stats::setNames(as.numeric(bins), names(level_weights))
    out <- rep(NA_real_, length(x))
    out[valid] <- unname(level_map[as.character(x_valid)])
    return(out)
  }
  
  ranks <- rank(x, ties.method = "average", na.last = "keep")
  n_valid <- sum(!is.na(ranks))
  out <- floor((ranks - 1) / n_valid * n_quantiles)
  pmin(out, n_quantiles - 1)
}


#' @keywords internal
.normalize_break_spec <- function(break_spec) {
  if (is.list(break_spec) && !is.null(break_spec$breaks)) {
    return(break_spec)
  }
  
  list(
    breaks = break_spec,
    fallback = "cut"
  )
}
