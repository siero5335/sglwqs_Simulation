#' @rdname plot_weights
#' @export
plot_weights.sglwqs <- function(object,
                         top_n = 10,
                         sort_by = c("weight", "selection_freq", "combined"),
                         pos_color = "#0072B2",
                         neg_color = "#E69F00",
                         title = "Variable Weights",
                         base_size = 11,
                         show_group_facet = TRUE,
                         ...) {
  
  sort_by <- match.arg(sort_by)
  
  # Fall back to weight sorting if bootstrap is unavailable
  if (sort_by %in% c("selection_freq", "combined") && !object$bootstrap) {
    message("Bootstrap not performed. Using sort_by = 'weight'.")
    sort_by <- "weight"
  }
  
  # Compute variable info
  var_info <- .compute_weight_var_order(object, sort_by = sort_by, top_n = top_n)
  
  # Create butterfly chart
  p <- .create_butterfly_chart(
    object = object,
    var_info = var_info,
    pos_color = pos_color,
    neg_color = neg_color,
    title = title,
    base_size = base_size,
    show_group_facet = show_group_facet,
    chart_type = "weights"
  )
  
  return(p)
}


#' Compute Variable Order for Weight Plots
#' @keywords internal
.compute_weight_var_order <- function(object, sort_by = "weight", top_n = 10) {
  
  df <- data.frame(
    variable = object$var_names,
    pos_weight = as.numeric(object$pos_weights),
    neg_weight = as.numeric(object$neg_weights),
    stringsAsFactors = FALSE
  )
  
  # Add bootstrap info if available
  if (object$bootstrap && !is.null(object$boot_info)) {
    df$selection_freq_pos <- object$boot_info$selection_freq_pos[df$variable]
    df$selection_freq_neg <- object$boot_info$selection_freq_neg[df$variable]
    df$selection_freq <- pmax(df$selection_freq_pos, df$selection_freq_neg, na.rm = TRUE)
  } else {
    df$selection_freq <- 1
    df$selection_freq_pos <- 1
    df$selection_freq_neg <- 1
  }
  
  # Compute sort criterion
  df$weight_max <- pmax(df$pos_weight, df$neg_weight)
  
  if (sort_by == "combined") {
    df$sort_val <- df$selection_freq * df$weight_max
  } else if (sort_by == "selection_freq") {
    df$sort_val <- df$selection_freq
  } else {
    df$sort_val <- df$weight_max
  }
  
  # Sort within groups if group info is available
  if (!is.null(object$groups)) {
    df$group <- .assign_group_labels(df$variable, object$groups)
    group_order <- .group_display_order(object$groups, unique(df$group))

    # Select top_n per group
    if (!is.null(top_n)) {
      df_list <- lapply(group_order, function(grp) {
        g <- df[df$group == grp, ]
        g <- g[order(-g$sort_val), ]
        head(g, top_n)
      })
      df_top <- do.call(rbind, df_list)
    } else {
      df_top <- df[order(match(df$group, group_order), -df$sort_val), ]
    }
  } else {
    df_top <- df[order(-df$sort_val), ]
    if (!is.null(top_n)) {
      df_top <- head(df_top, top_n)
    }
    df_top$group <- "All"
  }
  
  # Finalize variable order
  var_order <- df_top$variable
  
  return(list(
    var_order = var_order,
    df = df_top
  ))
}


#' Create Butterfly Chart
#' @keywords internal
.create_butterfly_chart <- function(object, var_info, 
                                     pos_color, neg_color,
                                     title, base_size, 
                                     show_group_facet = TRUE,
                                     chart_type = c("weights", "selection")) {
  
  chart_type <- match.arg(chart_type)
  var_order <- var_info$var_order
  
  if (chart_type == "weights") {
    # Weight data
    df <- data.frame(
      variable = object$var_names,
      pos_val = as.numeric(object$pos_weights),
      neg_val = as.numeric(object$neg_weights),
      stringsAsFactors = FALSE
    )
    y_label <- "Negative                    Positive"
    axis_format <- function(x) sprintf("%.2f", abs(x))
  } else {
    # Selection frequency data
    if (!object$bootstrap || is.null(object$boot_info)) {
      stop("Bootstrap information not available")
    }
    df <- data.frame(
      variable = names(object$boot_info$selection_freq_pos),
      pos_val = as.numeric(object$boot_info$selection_freq_pos),
      neg_val = as.numeric(object$boot_info$selection_freq_neg),
      stringsAsFactors = FALSE
    )
    y_label <- "Negative                    Positive"
    axis_format <- function(x) sprintf("%.0f%%", abs(x) * 100)
  }
  
  # Filter and merge group info
  df <- df[df$variable %in% var_order, ]
  group_info <- var_info$df[, c("variable", "group")]
  df <- merge(df, group_info, by = "variable", all.x = TRUE)
  
  # Check if multiple groups exist
  has_groups <- length(unique(df$group)) > 1 && show_group_facet
  
  # Preserve group order
  if (has_groups) {
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
  
  # Determine axis range
  if (chart_type == "weights") {
    max_val <- max(c(df$pos_val, df$neg_val), na.rm = TRUE)
    axis_limit <- ceiling(max_val * 10) / 10
    axis_limit <- max(axis_limit, 0.1)
  } else {
    axis_limit <- 1
  }
  
  # Create plot
  p <- ggplot2::ggplot(df) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = -.data$neg_val),
      stat = "identity", fill = neg_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::geom_bar(
      ggplot2::aes(x = .data$variable, y = .data$pos_val),
      stat = "identity", fill = pos_color, width = 0.7, alpha = 0.8
    ) +
    ggplot2::coord_flip() +
    ggplot2::geom_hline(yintercept = 0, color = "gray40", linewidth = 0.5)
  
  # Add 50% reference line for selection frequency
  if (chart_type == "selection") {
    p <- p + ggplot2::geom_hline(
      yintercept = c(-0.5, 0.5), 
      linetype = "dashed", 
      color = "gray60", 
      linewidth = 0.3
    )
  }
  
  p <- p +
    ggplot2::scale_y_continuous(
      limits = c(-axis_limit, axis_limit),
      labels = axis_format
    ) +
    ggplot2::labs(x = NULL, y = y_label, title = title) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
  
  # Add facets for multiple groups
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
