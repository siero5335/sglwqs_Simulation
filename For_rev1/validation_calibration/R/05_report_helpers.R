read_validation_table <- function(stem, profile) {
  path <- vc_file("output", "tables", paste0(stem, "_", profile, ".csv"))
  if (!file.exists(path)) {
    return(data.frame())
  }
  readr::read_csv(path, show_col_types = FALSE)
}

empty_note <- function(label) {
  paste0("No rows available for ", label, ". Run the pipeline first or check failed fits.")
}

kable_if_rows <- function(x, caption = NULL, digits = 3) {
  if (nrow(x) == 0L) {
    cat(empty_note(caption %||% "this table"), "\n")
    return(invisible(NULL))
  }
  print(knitr::kable(x, caption = caption, digits = digits))
  invisible(x)
}

compact_columns <- function(x, cols) {
  cols <- intersect(cols, names(x))
  x[, cols, drop = FALSE]
}

safe_rate_plot <- function(x, xvar, yvar, color = "group", title = NULL) {
  if (nrow(x) == 0L) {
    return(NULL)
  }
  ggplot2::ggplot(
    x,
    ggplot2::aes(x = .data[[xvar]], y = .data[[yvar]], color = .data[[color]])
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(.data[[color]], direction)), linewidth = 0.4) +
    ggplot2::facet_wrap(~ direction) +
    ggplot2::labs(title = title, x = xvar, y = yvar, color = color) +
    ggplot2::theme_minimal(base_size = 12)
}
