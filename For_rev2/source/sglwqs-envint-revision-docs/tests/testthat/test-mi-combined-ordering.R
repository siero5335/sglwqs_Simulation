test_that("MI helper preserves weight-based order and retains negative-only variables", {
  fit_mice <- structure(
    list(
      m = 5,
      n_successful = 5,
      validation = FALSE,
      groups = list(
        GroupA = c("x1", "x2"),
        GroupB = c("x3")
      ),
      pooled = list(
        var_names = c("x1", "x2", "x3"),
        pos_weights = c(x1 = 0.30, x2 = 0.10, x3 = 0.00),
        neg_weights = c(x1 = 0.00, x2 = 0.05, x3 = 0.40),
        pos_coef = c(x1 = 0.30, x2 = 0.10, x3 = 0.00),
        neg_coef = c(x1 = 0.00, x2 = -0.05, x3 = -0.40),
        mi_selection_freq_pos = c(x1 = 0.50, x2 = 0.90, x3 = 0.00),
        mi_selection_freq_neg = c(x1 = 0.00, x2 = 0.20, x3 = 0.80),
        pos_weights_var = c(x1 = 0, x2 = 0, x3 = 0),
        neg_weights_var = c(x1 = 0, x2 = 0, x3 = 0)
      )
    ),
    class = "sglwqs_mids"
  )
  
  var_info <- sglwqs:::.compute_mi_unified_var_order(
    fit_mice,
    top_n = NULL,
    select_by = "combined",
    order_by = "weight"
  )
  
  expect_equal(var_info$var_order, c("x1", "x2", "x3"))
  expect_true("x3" %in% var_info$var_order)
})


test_that("MI weight plot keeps negative-only variables when ordered by shared var_info", {
  fit_mice <- structure(
    list(
      m = 5,
      n_successful = 5,
      validation = FALSE,
      groups = NULL,
      pooled = list(
        var_names = c("x1", "x2", "x3"),
        pos_weights = c(x1 = 0.35, x2 = 0.00, x3 = 0.00),
        neg_weights = c(x1 = 0.00, x2 = 0.25, x3 = 0.40),
        pos_coef = c(x1 = 0.35, x2 = 0.00, x3 = 0.00),
        neg_coef = c(x1 = 0.00, x2 = -0.25, x3 = -0.40),
        mi_selection_freq_pos = c(x1 = 0.70, x2 = 0.00, x3 = 0.00),
        mi_selection_freq_neg = c(x1 = 0.00, x2 = 0.60, x3 = 0.90),
        pos_weights_var = c(x1 = 0, x2 = 0, x3 = 0),
        neg_weights_var = c(x1 = 0, x2 = 0, x3 = 0)
      )
    ),
    class = "sglwqs_mids"
  )
  
  var_info <- sglwqs:::.compute_mi_unified_var_order(
    fit_mice,
    top_n = NULL,
    select_by = "weight",
    order_by = "weight"
  )
  
  p <- plot_weights(fit_mice, top_n = NULL, sort_by = "weight", var_info = var_info)
  
  expect_s3_class(p, "ggplot")
  expect_true(any(as.character(p$data$variable) == "x3"))
  expect_false(any(is.na(p$data$variable)))
})
