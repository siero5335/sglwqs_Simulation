make_validation_mids_fixture <- function() {
  coef_mat1 <- matrix(
    c(
      0.40, 0.10, 0.02,
     -0.25, 0.12, 0.04,
      0.55, 0.20, 0.01
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("wqs_pos", "wqs_neg", "age"), c("Estimate", "Std. Error", "Pr(>|z|)"))
  )

  coef_mat2 <- matrix(
    c(
      0.36, 0.11, 0.03,
     -0.22, 0.13, 0.06,
      0.48, 0.18, 0.02
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("wqs_pos", "wqs_neg", "age"), c("Estimate", "Std. Error", "Pr(>|z|)"))
  )

  fit1 <- structure(
    list(
      validation = TRUE,
      cov_names = "age",
      validation_info = list(
        coef_table = coef_mat1,
        p_col = "Pr(>|z|)"
      )
    ),
    class = "sglwqs"
  )

  fit2 <- structure(
    list(
      validation = TRUE,
      cov_names = "age",
      validation_info = list(
        coef_table = coef_mat2,
        p_col = "Pr(>|z|)"
      )
    ),
    class = "sglwqs"
  )

  structure(
    list(
      m = 2,
      n_successful = 2,
      validation = TRUE,
      fits = list(fit1, fit2),
      pooled = list(
        var_names = c("x1", "x2"),
        pos_weights = c(x1 = 0.6, x2 = 0.1),
        neg_weights = c(x1 = 0.0, x2 = 0.3),
        pos_coef = c(x1 = 0.6, x2 = 0.1),
        neg_coef = c(x1 = 0.0, x2 = -0.3),
        mi_selection_freq_pos = c(x1 = 1.0, x2 = 0.4),
        mi_selection_freq_neg = c(x1 = 0.0, x2 = 0.8),
        pos_weights_var = c(x1 = 0.0, x2 = 0.0),
        neg_weights_var = c(x1 = 0.0, x2 = 0.0),
        validation = list(
          has_group_results = FALSE,
          group_results = NULL,
          wqs_pos = list(estimate = 0.38, se = 0.10, p_value = 0.02, fmi = 0.10, df = 12),
          wqs_neg = list(estimate = -0.24, se = 0.11, p_value = 0.05, fmi = 0.12, df = 11)
        )
      )
    ),
    class = "sglwqs_mids"
  )
}


test_that("sglwqs_mids validation forest plot can include pooled covariates", {
  fit_mice <- make_validation_mids_fixture()

  p <- plot_validation_results(
    fit_mice,
    include_covariates = TRUE,
    show_significance = "asterisk"
  )

  expect_s3_class(p, "ggplot")
  expect_true("Covariate" %in% p$data$group)
  expect_true(any(grepl("age", as.character(p$data$label))))
  expect_true("sig_label" %in% names(p$data))
})


test_that("sglwqs_mids combined plot accepts validation-panel options", {
  fit_mice <- make_validation_mids_fixture()

  p <- plot_combined_results(
    fit_mice,
    layout = "horizontal",
    conf_level = 0.90,
    include_covariates = TRUE,
    show_significance = "pvalue"
  )

  expect_true(
    inherits(p, "patchwork") ||
      inherits(p, "ggplot") ||
      is.list(p)
  )

  if (is.list(p)) {
    expect_true("validation" %in% names(p))
  }
})
