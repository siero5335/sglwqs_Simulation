#' S3 Generic Functions for sglwqs
#'
#' Unified API generics for both \code{sglwqs} and \code{sglwqs_mids}
#' objects.


# =============================================================================
# extract_weights - Unified Weight Extraction
# =============================================================================

#' Extract Weights from sglwqs Objects
#'
#' Generic function to extract weights from sglwqs or sglwqs_mids objects.
#' This is the unified interface for weight extraction - it automatically
#' handles both regular sglwqs objects and MICE-pooled results.
#'
#' @param object A sglwqs or sglwqs_mids object.
#' @param direction Character. "both" (default), "positive", or "negative".
#' @param by_group Logical. If TRUE, include group information in output.
#'   Default is TRUE if groups were specified in the model.
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with variable weights and related information.
#'   For sglwqs objects: variable, weight, coefficient, direction, and optional
#'   columns such as group, selection_freq, and se.
#'   For sglwqs_mids objects: variable, weight, coefficient, mi_selection_freq,
#'   direction, and optional group.
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
#' extract_weights(fit)
#' extract_weights(fit, direction = "positive")
#' extract_weights(fit, by_group = TRUE)
#' }
#'
#' \dontrun{
#' # MICE sglwqs (same interface!)
#' fit_mice <- sglwqs_mice(imp, ...)
#' extract_weights(fit_mice)
#' }
#'
#' @export
extract_weights <- function(object, ...) {
  UseMethod("extract_weights")
}


#' @rdname extract_weights
#' @export
extract_weights.sglwqs <- function(object, 
                                    direction = c("both", "positive", "negative"),
                                    by_group = NULL,
                                    ...) {
  
  direction <- match.arg(direction)
  
  # Default by_group: TRUE if groups exist
  if (is.null(by_group)) {
    by_group <- !is.null(object$groups)
  }
  
  # Get variable names (handle missing names)
  var_names <- object$var_names
  if (is.null(var_names) || length(var_names) == 0) {
    # Try to get from pos_weights names
    var_names <- names(object$pos_weights)
    if (is.null(var_names) || length(var_names) == 0) {
      var_names <- names(object$neg_weights)
    }
  }
  
  # If variable names still unavailable
  if (is.null(var_names) || length(var_names) == 0) {
    warning("No variable names found in object")
    return(data.frame(
      variable = character(0),
      weight = numeric(0),
      coefficient = numeric(0),
      direction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  n_vars <- length(var_names)
  
  # Get weights and coefficients
  # NULL/empty → silent zero (field absent in legacy objects)
  # Present but wrong length → stop (data corruption)
  .safe_vec <- function(x, n, label) {
    v <- as.numeric(x)
    if (is.null(v) || length(v) == 0) return(rep(0, n))
    if (length(v) != n) {
      stop("Length mismatch in extract_weights.sglwqs(): `", label,
           "` has length ", length(v), ", expected ", n, ".")
    }
    v
  }
  pos_weights <- .safe_vec(object$pos_weights, n_vars, "pos_weights")
  neg_weights <- .safe_vec(object$neg_weights, n_vars, "neg_weights")
  pos_coef    <- .safe_vec(object$pos_coef,    n_vars, "pos_coef")
  neg_coef    <- .safe_vec(object$neg_coef,    n_vars, "neg_coef")
  
  # Positive direction data
  pos_df <- data.frame(
    variable = var_names,
    weight = pos_weights,
    coefficient = pos_coef,
    direction = "positive",
    stringsAsFactors = FALSE
  )
  
  # Negative direction data
  neg_df <- data.frame(
    variable = var_names,
    weight = neg_weights,
    coefficient = neg_coef,
    direction = "negative",
    stringsAsFactors = FALSE
  )
  
  # Add bootstrap info if available
  if (object$bootstrap && !is.null(object$boot_info)) {
    sel_pos <- as.numeric(object$boot_info$selection_freq_pos)
    sel_neg <- as.numeric(object$boot_info$selection_freq_neg)
    se_pos <- as.numeric(object$boot_info$se_pos_coef)
    se_neg <- as.numeric(object$boot_info$se_neg_coef)
    
    if (length(sel_pos) == n_vars) pos_df$selection_freq <- sel_pos
    if (length(se_pos) == n_vars) pos_df$se <- se_pos
    if (length(sel_neg) == n_vars) neg_df$selection_freq <- sel_neg
    if (length(se_neg) == n_vars) neg_df$se <- se_neg
  }
  
  # Group info
  if (by_group && !is.null(object$groups)) {
    pos_df$group <- .assign_group_labels(pos_df$variable, object$groups)
    neg_df$group <- .assign_group_labels(neg_df$variable, object$groups)
  }

  # Filter by direction
  if (direction == "positive") {
    return(pos_df)
  } else if (direction == "negative") {
    return(neg_df)
  } else {
    return(rbind(pos_df, neg_df))
  }
}


#' @rdname extract_weights
#' @export
extract_weights.sglwqs_mids <- function(object, 
                                         direction = c("both", "positive", "negative"),
                                         by_group = NULL,
                                         ...) {
  
  direction <- match.arg(direction)
  
  # Default by_group
  if (is.null(by_group)) {
    by_group <- !is.null(.get_effective_groups(object))
  }

  pooled <- object$pooled
  
  # Check pooled object existence
  if (is.null(pooled)) {
    warning("No pooled results found in sglwqs_mids object")
    return(data.frame(
      variable = character(0),
      weight = numeric(0),
      coefficient = numeric(0),
      mi_selection_freq = numeric(0),
      weight_var = numeric(0),
      direction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  # Get variable names
  var_names <- pooled$var_names
  if (is.null(var_names) || length(var_names) == 0) {
    var_names <- names(pooled$pos_weights)
    if (is.null(var_names) || length(var_names) == 0) {
      var_names <- names(pooled$neg_weights)
    }
  }
  
  # If variable names still unavailable
  if (is.null(var_names) || length(var_names) == 0) {
    warning("No variable names found in pooled results")
    return(data.frame(
      variable = character(0),
      weight = numeric(0),
      coefficient = numeric(0),
      mi_selection_freq = numeric(0),
      weight_var = numeric(0),
      direction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  n_vars <- length(var_names)
  
  # Function to safely get vectors
  # All pooled vectors should be present once an sglwqs_mids object is constructed.
  # Missing or malformed vectors indicate an inconsistent pooled object.
  safe_numeric <- function(x, n) {
    result <- as.numeric(x)
    if (is.null(result) || length(result) == 0) {
      stop("Missing pooled vector in extract_weights.sglwqs_mids().")
    }
    if (length(result) == 1L) {
      return(rep(result, n))
    }
    if (length(result) != n) {
      stop("Length mismatch in extract_weights.sglwqs_mids(): expected ",
           n, ", got ", length(result), ".")
    }
    return(result)
  }
  
  # Positive direction data
  pos_df <- data.frame(
    variable = var_names,
    weight = safe_numeric(pooled$pos_weights, n_vars),
    coefficient = safe_numeric(pooled$pos_coef, n_vars),
    mi_selection_freq = safe_numeric(pooled$mi_selection_freq_pos, n_vars),
    weight_var = safe_numeric(pooled$pos_weights_var, n_vars),
    direction = "positive",
    stringsAsFactors = FALSE
  )

  # Negative direction data
  neg_df <- data.frame(
    variable = var_names,
    weight = safe_numeric(pooled$neg_weights, n_vars),
    coefficient = safe_numeric(pooled$neg_coef, n_vars),
    mi_selection_freq = safe_numeric(pooled$mi_selection_freq_neg, n_vars),
    weight_var = safe_numeric(pooled$neg_weights_var, n_vars),
    direction = "negative",
    stringsAsFactors = FALSE
  )
  
  # Add bootstrap information if available
  if (!is.null(pooled$bootstrap)) {
    boot_pos <- safe_numeric(pooled$bootstrap$selection_freq_pos, n_vars)
    boot_neg <- safe_numeric(pooled$bootstrap$selection_freq_neg, n_vars)
    pos_df$boot_selection_freq <- boot_pos
    neg_df$boot_selection_freq <- boot_neg
  }
  
  # Group information (single source of truth)
  groups_to_use <- .get_effective_groups(object)
  if (by_group && !is.null(groups_to_use)) {
    pos_df$group <- .assign_group_labels(pos_df$variable, groups_to_use)
    neg_df$group <- .assign_group_labels(neg_df$variable, groups_to_use)
  }
  
  # Filter by direction
  if (direction == "positive") {
    return(pos_df)
  } else if (direction == "negative") {
    return(neg_df)
  } else {
    return(rbind(pos_df, neg_df))
  }
}


# =============================================================================
# plot_weights - Unified Weight Visualization
# =============================================================================

#' Plot Variable Weights
#'
#' Generic function to create butterfly charts of variable weights.
#' Works with both sglwqs and sglwqs_mids objects - just use the same
#' function name regardless of input type.
#'
#' @param object A sglwqs or sglwqs_mids object.
#' @param top_n Number of top variables to show per group (default: 10).
#' @param sort_by Character. Sorting criterion. For sglwqs: "weight", "selection_freq", 
#'   or "combined". For sglwqs_mids: "combined", "weight", or "mi_selection_freq".
#' @param pos_color Color for positive direction (default: "#0072B2", Okabe-Ito blue).
#' @param neg_color Color for negative direction (default: "#E69F00", Okabe-Ito orange).
#' @param show_freq Whether to show selection frequency labels (default: TRUE for MICE).
#' @param show_group_facet Logical. Facet by group in \code{plot_weights.sglwqs}
#'   when groups are available (default: TRUE).
#' @param base_size Base font size (default: 11).
#' @param title Plot title (default: auto-generated).
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
#' plot_weights(fit)
#' plot_weights(fit, top_n = 5, sort_by = "selection_freq")
#' }
#'
#' \dontrun{
#' # MICE sglwqs (same function!)
#' fit_mice <- sglwqs_mice(imp, ...)
#' plot_weights(fit_mice)
#' }
#'
#' @export
plot_weights <- function(object, ...) {
  UseMethod("plot_weights")
}


#' @rdname plot_weights
#' @export
plot_weights.sglwqs_mids <- function(object, 
                                      top_n = 10, 
                                      sort_by = c("combined", "weight", "mi_selection_freq"),
                                      pos_color = "#0072B2",
                                      neg_color = "#E69F00",
                                      show_freq = TRUE,
                                      base_size = 11,
                                      title = NULL,
                                      ...) {
  
  sort_by <- match.arg(sort_by)
  extra_args <- list(...)
  var_info <- extra_args$var_info
  
  # Default title
  if (is.null(title)) {
    title <- paste0("Pooled Weights from MICE (m=", object$m, ")")
  }
  
  if (is.null(var_info)) {
    var_info <- .compute_mi_unified_var_order(
      object,
      top_n = top_n,
      select_by = sort_by,
      order_by = sort_by
    )
  }
  
  if (length(var_info$var_order) == 0) {
    message("No non-zero weights to plot")
    return(invisible(NULL))
  }
  
  # Get pooled weights
  weights_df <- extract_weights(object, direction = "both", by_group = TRUE)
  
  # Check for empty dataframe
  if (is.null(weights_df) || nrow(weights_df) == 0) {
    message("No weight data available to plot")
    return(invisible(NULL))
  }
  
  # Fill with 1 if mi_selection_freq doesn't exist
  if (!"mi_selection_freq" %in% names(weights_df)) {
    weights_df$mi_selection_freq <- 1
  }
  
  plot_data <- weights_df[
    weights_df$variable %in% var_info$var_order & weights_df$weight > 0,
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_data) == 0) {
    message("No non-zero weights to plot")
    return(invisible(NULL))
  }
  
  # Invert negative direction for butterfly plot
  plot_data$weight_butterfly <- ifelse(
    plot_data$direction == "negative",
    -plot_data$weight,
    plot_data$weight
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
  
  # MI selection frequency label
  if (show_freq) {
    plot_data$freq_label <- sprintf("%.0f%%", plot_data$mi_selection_freq * 100)
    plot_data$freq_label[plot_data$weight == 0] <- ""
  }

  # Set symmetric axes, display labels as absolute values
  max_weight <- max(plot_data$weight, na.rm = TRUE)
  axis_limit <- ceiling(max_weight * 10) / 10
  axis_limit <- max(axis_limit, 0.1)
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$variable, y = .data$weight_butterfly))
  
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
    ggplot2::scale_y_continuous(
      limits = c(-axis_limit, axis_limit),
      labels = function(x) sprintf("%.2f", abs(x)),
      expand = ggplot2::expansion(mult = c(0.12, 0.12))
    ) +
    ggplot2::scale_fill_manual(
      values = c("positive" = pos_color, "negative" = neg_color),
      labels = c("positive" = "Positive", "negative" = "Negative"),
      name = "Direction"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Pooled from ", object$n_successful, " imputations"),
      x = NULL,
      y = "Negative                    Positive"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "gray95", color = NA),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )
  
  # Frequency label
  if (show_freq) {
    pos_label_data <- plot_data[plot_data$direction == "positive" & plot_data$weight > 0, ]
    if (nrow(pos_label_data) > 0) {
      p <- p + ggplot2::geom_text(
        data = pos_label_data,
        ggplot2::aes(label = .data$freq_label),
        hjust = -0.1, size = 3, color = "gray30"
      )
    }
    neg_label_data <- plot_data[plot_data$direction == "negative" & plot_data$weight > 0, ]
    if (nrow(neg_label_data) > 0) {
      p <- p + ggplot2::geom_text(
        data = neg_label_data,
        ggplot2::aes(label = .data$freq_label),
        hjust = 1.1, size = 3, color = "gray30"
      )
    }
  }
  
  return(p)
}


# =============================================================================
# inference_table - Unified Publication-Ready Tables
# =============================================================================

#' Create Inference Table for Publication
#'
#' Generic function to create publication-ready inference tables.
#' Works with both sglwqs and sglwqs_mids objects.
#'
#' @param object A sglwqs or sglwqs_mids object.
#' @param min_freq Minimum selection frequency to include (default: 0).
#' @param top_n Number of top variables to show (default: NULL for all).
#' @param conf_level Confidence level (default: 0.95).
#' @param ... Additional arguments passed to methods.
#'
#' @return A data frame with inference results, formatted for publication.
#'   The output includes weights, selection frequencies, and validation results
#'   (if available) as attributes.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:240, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:240],
#'   bootstrap = TRUE,
#'   n_boot = 10,
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' inference_table(fit)
#' }
#'
#' \dontrun{
#' # For MICE objects
#' fit_mice <- sglwqs_mice(imp, ...)
#' inference_table(fit_mice)
#' }
#'
#' @export
inference_table <- function(object, ...) {
  UseMethod("inference_table")
}


#' @rdname inference_table
#' @export
inference_table.sglwqs <- function(object, 
                                    min_freq = 0,
                                    top_n = NULL,
                                    conf_level = 0.95,
                                    ...) {
  
  var_names <- object$var_names
  
  # Wide format dataframe
  result <- data.frame(
    variable = var_names,
    pos_weight = as.numeric(object$pos_weights),
    neg_weight = as.numeric(object$neg_weights),
    stringsAsFactors = FALSE
  )
  
  # Group information
  if (!is.null(object$groups)) {
    result$group <- .assign_group_labels(result$variable, object$groups)
  }

  # Bootstrap information
  if (object$bootstrap && !is.null(object$boot_info)) {
    result$pos_sel_freq <- as.numeric(object$boot_info$selection_freq_pos[var_names])
    result$neg_sel_freq <- as.numeric(object$boot_info$selection_freq_neg[var_names])
    result$pos_coef_se <- as.numeric(object$boot_info$se_pos_coef[var_names])
    result$neg_coef_se <- as.numeric(object$boot_info$se_neg_coef[var_names])
    
    # Importance score
    result$importance_pos <- result$pos_weight * result$pos_sel_freq
    result$importance_neg <- result$neg_weight * result$neg_sel_freq
  } else {
    result$importance_pos <- result$pos_weight
    result$importance_neg <- result$neg_weight
  }
  
  # Filtering
  if (object$bootstrap && !is.null(object$boot_info)) {
    keep <- (result$pos_sel_freq >= min_freq & result$pos_weight > 0) |
            (result$neg_sel_freq >= min_freq & result$neg_weight > 0)
  } else {
    keep <- result$pos_weight > 0 | result$neg_weight > 0
  }
  result <- result[keep, ]
  
  if (nrow(result) == 0) {
    message("No variables meet the selection criteria")
    return(NULL)
  }
  
  # Sort
  result$max_importance <- pmax(result$importance_pos, result$importance_neg, na.rm = TRUE)
  result <- result[order(-result$max_importance), ]
  result$max_importance <- NULL
  
  # top_n
  if (!is.null(top_n) && nrow(result) > top_n) {
    result <- head(result, top_n)
  }
  
  rownames(result) <- NULL
  
  # Add validation information as attributes
  val <- .get_active_refit_info(object)
  if (!is.null(val)) {
    if (!is.null(val$group_results)) {
      attr(result, "validation_type") <- "group"
      attr(result, "group_results") <- val$group_results
    } else {
      attr(result, "validation_type") <- "overall"
      attr(result, "wqs_pos_estimate") <- val$wqs_pos_estimate
      attr(result, "wqs_pos_se") <- val$wqs_pos_se
      attr(result, "wqs_pos_pvalue") <- val$wqs_pos_pvalue
      attr(result, "wqs_neg_estimate") <- val$wqs_neg_estimate
      attr(result, "wqs_neg_se") <- val$wqs_neg_se
      attr(result, "wqs_neg_pvalue") <- val$wqs_neg_pvalue
    }
    attr(result, "n_train") <- val$n_train
    attr(result, "n_val") <- val$n_val
    attr(result, "inference_source") <- attr(val, "source")
  }
  
  attr(result, "bootstrap") <- object$bootstrap
  attr(result, "validation") <- object$validation
  attr(result, "conf_level") <- conf_level
  
  class(result) <- c("sglwqs_inference_table", "data.frame")
  
  return(result)
}


#' @rdname inference_table
#' @export
inference_table.sglwqs_mids <- function(object, 
                                         min_freq = 0,
                                         top_n = NULL,
                                         conf_level = 0.95,
                                         ...) {
  
  pooled <- object$pooled
  var_names <- pooled$var_names
  
  # Wide format dataframe
  result <- data.frame(
    variable = var_names,
    pos_weight = as.numeric(pooled$pos_weights),
    neg_weight = as.numeric(pooled$neg_weights),
    stringsAsFactors = FALSE
  )
  
  # Group information (single source of truth)
  eff_groups <- .get_effective_groups(object)
  if (!is.null(eff_groups)) {
    result$group <- .assign_group_labels(result$variable, eff_groups)
  }

  # MI selection frequency
  result$pos_mi_freq <- as.numeric(pooled$mi_selection_freq_pos[var_names])
  result$neg_mi_freq <- as.numeric(pooled$mi_selection_freq_neg[var_names])
  
  # Bootstrap selection frequency
  if (!is.null(pooled$bootstrap)) {
    result$pos_boot_freq <- as.numeric(pooled$bootstrap$selection_freq_pos[var_names])
    result$neg_boot_freq <- as.numeric(pooled$bootstrap$selection_freq_neg[var_names])
    result$pos_coef_se <- as.numeric(pooled$bootstrap$se_pos_coef[var_names])
    result$neg_coef_se <- as.numeric(pooled$bootstrap$se_neg_coef[var_names])
  }
  
  # Weight variance
  result$pos_weight_var <- as.numeric(pooled$pos_weights_var[var_names])
  result$neg_weight_var <- as.numeric(pooled$neg_weights_var[var_names])
  
  # Importance score
  result$importance_pos <- result$pos_weight * result$pos_mi_freq
  result$importance_neg <- result$neg_weight * result$neg_mi_freq
  
  # Filtering
  keep <- (result$pos_mi_freq >= min_freq & result$pos_weight > 0) |
          (result$neg_mi_freq >= min_freq & result$neg_weight > 0)
  result <- result[keep, ]
  
  if (nrow(result) == 0) {
    message("No variables meet the selection criteria")
    return(NULL)
  }
  
  # Sort
  result$max_importance <- pmax(result$importance_pos, result$importance_neg, na.rm = TRUE)
  result <- result[order(-result$max_importance), ]
  result$max_importance <- NULL
  
  # top_n
  if (!is.null(top_n) && nrow(result) > top_n) {
    result <- head(result, top_n)
  }
  
  rownames(result) <- NULL
  
  # Add validation information as attributes
  val <- .get_mi_pooled_inference(object)
  if (!is.null(val)) {
    if (isTRUE(val$has_group_results) && !is.null(val$group_results)) {
      attr(result, "validation_type") <- "group"
      attr(result, "group_results") <- val$group_results
    } else {
      attr(result, "validation_type") <- "overall"
      attr(result, "wqs_pos") <- val$wqs_pos
      attr(result, "wqs_neg") <- val$wqs_neg
    }
    attr(result, "n_train") <- val$n_train
    attr(result, "n_val") <- val$n_val
    attr(result, "inference_source") <- val$source_used %||% "validation_info"
  }
  
  attr(result, "m") <- object$m
  attr(result, "n_successful") <- object$n_successful
  attr(result, "conf_level") <- conf_level
  
  class(result) <- c("sglwqs_mids_inference", "data.frame")
  
  return(result)
}


#' Print Method for sglwqs_inference_table
#'
#' @param x An object of class \code{"sglwqs_inference_table"}.
#' @param ... Additional arguments passed to \code{print.data.frame()}.
#' @export
print.sglwqs_inference_table <- function(x, ...) {
  cat("SGL-WQS Inference Table\n")
  cat(strrep("=", 60), "\n\n")
  
  # Validation information
  val_type <- attr(x, "validation_type")
  if (!is.null(val_type) && val_type == "group") {
    group_results <- attr(x, "group_results")
    if (!is.null(group_results)) {
      cat("Mixture Effect Estimates (Validation):\n")
      cat(strrep("-", 60), "\n")
      for (grp in names(group_results)) {
        r <- group_results[[grp]]
        cat(sprintf("  %s:\n", grp))
        if (!is.na(r$pos_estimate) && r$pos_estimate != 0) {
          cat(sprintf("    Positive: %.4f (SE=%.4f, p=%s)\n",
                      r$pos_estimate, r$pos_se,
                      format.pval(r$pos_pvalue, digits = 3)))
        }
        if (!is.na(r$neg_estimate) && r$neg_estimate != 0) {
          cat(sprintf("    Negative: %.4f (SE=%.4f, p=%s)\n",
                      r$neg_estimate, r$neg_se,
                      format.pval(r$neg_pvalue, digits = 3)))
        }
      }
      cat("\n")
    }
  } else if (!is.null(val_type) && val_type == "overall") {
    cat("Overall Mixture Effects (Validation):\n")
    pos_est <- attr(x, "wqs_pos_estimate")
    neg_est <- attr(x, "wqs_neg_estimate")
    if (!is.null(pos_est) && !is.na(pos_est)) {
      cat(sprintf("  Positive Index: %.4f (p=%s)\n",
                  pos_est, format.pval(attr(x, "wqs_pos_pvalue"), digits = 3)))
    }
    if (!is.null(neg_est) && !is.na(neg_est)) {
      cat(sprintf("  Negative Index: %.4f (p=%s)\n",
                  neg_est, format.pval(attr(x, "wqs_neg_pvalue"), digits = 3)))
    }
    cat("\n")
  }
  
  cat("Variable Weights:\n")
  cat(strrep("-", 60), "\n")
  
  # Display the dataframe portion
  print.data.frame(x, row.names = FALSE, ...)
  
  invisible(x)
}


# =============================================================================
# plot_validation_results - Validation / Mixture Effect Forest Plots
# =============================================================================

#' Plot Validation Results
#'
#' Generic function to visualize validation / mixture-effect estimates.
#' Works with both \code{sglwqs} and \code{sglwqs_mids} objects.
#'
#' @param object A \code{sglwqs} or \code{sglwqs_mids} object with downstream
#'   inference available via \code{validation = TRUE} or \code{refit = "full"}.
#' @param ... Additional arguments passed to methods.
#'
#' @return A ggplot object.
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:240, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:240],
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' plot_validation_results(fit)
#' }
#'
#' \dontrun{
#' # Rubin-pooled downstream forest plot for multiply imputed fits
#' fit_mice <- sglwqs_mice(imp, ..., refit = "full")
#' plot_validation_results(fit_mice)
#' }
#'
#' @export
plot_validation_results <- function(object, ...) {
  UseMethod("plot_validation_results")
}


# =============================================================================
# plot_combined_results - Unified Combined Visualization
# =============================================================================

#' Combined Results Plot (Publication-Ready)
#'
#' Generic function to create combined multi-panel plots for publication.
#' Works with both sglwqs and sglwqs_mids objects.
#'
#' @param object A sglwqs or sglwqs_mids object.
#' @param top_n Number of top variables to show per group (default: 10).
#' @param sort_by Character. Sorting criterion.
#' @param pos_color Color for positive direction (default: "#0072B2").
#' @param neg_color Color for negative direction (default: "#E69F00").
#' @param layout Layout option for \code{plot_combined_results.sglwqs}:
#'   "auto", "vertical", or "horizontal".
#' @param show_significance Significance annotation mode for
#'   \code{plot_combined_results.sglwqs}: "none", "asterisk", or "pvalue".
#' @param conf_level Confidence level for confidence-interval based panels
#'   (default: 0.95).
#' @param base_size Base font size (default: 11).
#' @param include_covariates Logical. Whether to include covariate coefficients
#'   in the validation forest plot panel (default: FALSE).
#' @param ... Additional arguments passed to methods.
#'
#' @return A combined ggplot object (requires patchwork package).
#'
#' @examples
#' \donttest{
#' data("sglwqs_example")
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[1:240, exp_vars],
#'   y = sglwqs_example$outcome_cont[1:240],
#'   bootstrap = TRUE,
#'   n_boot = 10,
#'   validation = TRUE,
#'   train_prop = 0.6,
#'   nfolds = 3,
#'   nlambda = 20,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' plot_combined_results(fit)
#' }
#'
#' \dontrun{
#' # For MICE objects
#' fit_mice <- sglwqs_mice(imp, ...)
#' plot_combined_results(fit_mice)
#' }
#'
#' @export
plot_combined_results <- function(object, ...) {
  UseMethod("plot_combined_results")
}


# =============================================================================
# Deprecated Functions (for backward compatibility)
# =============================================================================

#' Deprecated: Use extract_weights() instead
#' @export
#' @keywords internal
extract_weights_grouped <- function(object, ...) {
  .Deprecated("extract_weights(..., by_group = TRUE)", 
              msg = "extract_weights_grouped() is deprecated. Use extract_weights(..., by_group = TRUE) instead.")
  extract_weights(object, by_group = TRUE, ...)
}

#' Deprecated: Use extract_weights() instead
#' @export
#' @keywords internal
extract_pooled_weights <- function(object, ...) {
  .Deprecated("extract_weights", 
              msg = "extract_pooled_weights() is deprecated. Use extract_weights() instead.")
  extract_weights(object, ...)
}

#' Deprecated: Use plot_weights() instead
#' @export
#' @keywords internal  
plot_pooled_weights <- function(object, ...) {
  .Deprecated("plot_weights", 
              msg = "plot_pooled_weights() is deprecated. Use plot_weights() instead.")
  plot_weights(object, ...)
}

#' Deprecated: Use inference_table() instead
#' @export
#' @keywords internal
get_inference_table <- function(object, ...) {
  .Deprecated("inference_table", 
              msg = "get_inference_table() is deprecated. Use inference_table() instead.")
  inference_table(object, ...)
}

#' Deprecated: summary_by_group is no longer needed
#' 
#' Use summary(fit) to see group summaries, or extract_weights(fit, by_group = TRUE)
#' to get a data frame with group information.
#' @export
#' @keywords internal
summary_by_group <- function(object, ...) {
  .Deprecated("summary() or extract_weights(..., by_group = TRUE)",
              msg = "summary_by_group() is deprecated. Use summary() for display or extract_weights() for data extraction.")
  
  if (is.null(object$groups)) {
    message("No groups defined in the model.")
    return(NULL)
  }
  
  # Return aggregation by group
  group_summary <- do.call(rbind, lapply(names(object$groups), function(grp_name) {
    grp_vars <- object$groups[[grp_name]]
    pos_w <- object$pos_weights[grp_vars]
    neg_w <- object$neg_weights[grp_vars]
    
    data.frame(
      group = grp_name,
      n_vars = length(grp_vars),
      pos_index_sum = sum(object$pos_coef[grp_vars]),
      neg_index_sum = sum(object$neg_coef[grp_vars]),
      n_pos_nonzero = sum(pos_w > 0),
      n_neg_nonzero = sum(neg_w > 0),
      stringsAsFactors = FALSE
    )
  }))
  
  return(group_summary)
}
