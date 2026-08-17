test_that("plot_combined_results.sglwqs_mids accepts weight-based ordering when select_by differs", {
  fit_mice <- structure(
    list(
      m = 5,
      n_successful = 5,
      validation = FALSE,
      groups = list(GroupA = c("x1", "x2"), GroupB = "x3"),
      pooled = list(
        var_names = c("x1", "x2", "x3"),
        pos_weights = c(x1 = 0.30, x2 = 0.05, x3 = 0.00),
        neg_weights = c(x1 = 0.00, x2 = 0.15, x3 = 0.25),
        pos_coef = c(x1 = 0.30, x2 = 0.05, x3 = 0.00),
        neg_coef = c(x1 = 0.00, x2 = -0.15, x3 = -0.25),
        mi_selection_freq_pos = c(x1 = 0.45, x2 = 0.90, x3 = 0.00),
        mi_selection_freq_neg = c(x1 = 0.00, x2 = 0.20, x3 = 0.85),
        pos_weights_var = c(x1 = 0, x2 = 0, x3 = 0),
        neg_weights_var = c(x1 = 0, x2 = 0, x3 = 0)
      )
    ),
    class = "sglwqs_mids"
  )

  expect_no_error({
    res <- plot_combined_results(fit_mice, top_n = NULL, sort_by = "combined")
    expect_false(is.null(res))
  })
})
