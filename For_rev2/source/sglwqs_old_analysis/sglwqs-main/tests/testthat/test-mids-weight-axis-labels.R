test_that("sglwqs_mids weight plot uses absolute-value axis labels", {
  fit_mice <- structure(
    list(
      m = 3,
      n_successful = 3,
      family = "gaussian",
      bootstrap = FALSE,
      validation = FALSE,
      groups = list(GroupA = c("x1", "x2"), GroupB = "x3"),
      pooled = list(
        var_names = c("x1", "x2", "x3"),
        pos_weights = c(x1 = 0.40, x2 = 0.10, x3 = 0.00),
        neg_weights = c(x1 = 0.00, x2 = 0.20, x3 = 0.30),
        pos_coef = c(x1 = 0.50, x2 = 0.20, x3 = 0.00),
        neg_coef = c(x1 = 0.00, x2 = -0.25, x3 = -0.40),
        mi_selection_freq_pos = c(x1 = 0.80, x2 = 0.40, x3 = 0.00),
        mi_selection_freq_neg = c(x1 = 0.00, x2 = 0.50, x3 = 0.70),
        pos_weights_var = c(x1 = 0, x2 = 0, x3 = 0),
        neg_weights_var = c(x1 = 0, x2 = 0, x3 = 0)
      )
    ),
    class = "sglwqs_mids"
  )

  p <- plot_weights(fit_mice, show_freq = FALSE)

  expect_s3_class(p, "ggplot")
  expect_true(any(p$data$weight_butterfly < 0))
  expect_identical(p$labels$y, "Negative                    Positive")

  yscale <- p$scales$get_scales("y")
  expect_false(is.null(yscale))
  expect_equal(yscale$labels(c(-0.2, 0, 0.2)), c("0.20", "0.00", "0.20"))
})
