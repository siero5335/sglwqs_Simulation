make_summary_validation_mids_fixture <- function() {
  coef_mat1 <- matrix(
    c(
      0.80, 0.15, 0.001,
      0.42, 0.10, 0.020,
     -0.21, 0.11, 0.060,
      0.55, 0.12, 0.010,
     -0.18, 0.09, 0.040
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      c("(Intercept)", "wqs_pos", "wqs_neg", "age", "sexM"),
      c("Estimate", "Std. Error", "Pr(>|z|)")
    )
  )

  coef_mat2 <- matrix(
    c(
      0.75, 0.14, 0.002,
      0.36, 0.11, 0.030,
     -0.18, 0.12, 0.080,
      0.48, 0.10, 0.008,
     -0.15, 0.08, 0.050
    ),
    nrow = 5,
    byrow = TRUE,
    dimnames = list(
      c("(Intercept)", "wqs_pos", "wqs_neg", "age", "sexM"),
      c("Estimate", "Std. Error", "Pr(>|z|)")
    )
  )

  fit1 <- structure(
    list(
      validation = TRUE,
      intercept = 0.80,
      cov_coef = c(age = 0.55, sexM = -0.18),
      cov_names = c("age", "sexM"),
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
      intercept = 0.75,
      cov_coef = c(age = 0.48, sexM = -0.15),
      cov_names = c("age", "sexM"),
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
      bootstrap = FALSE,
      family = "gaussian",
      outcome_var = "y",
      exposure_vars = c("x1", "x2"),
      covariate_vars = c("age", "sexM"),
      penalize_covariates = FALSE,
      fits = list(fit1, fit2),
      pooled = list(
        pos_index_sum = 0.78,
        neg_index_sum = 0.31,
        pos_weights = c(x1 = 0.6, x2 = 0.2),
        neg_weights = c(x1 = 0.0, x2 = 0.4),
        mi_selection_freq_pos = c(x1 = 1.0, x2 = 0.5),
        mi_selection_freq_neg = c(x1 = 0.0, x2 = 1.0),
        validation = list(
          has_group_results = FALSE,
          group_results = NULL,
          wqs_pos = rubin_pool(c(0.42, 0.36), c(0.10^2, 0.11^2), 2),
          wqs_neg = rubin_pool(c(-0.21, -0.18), c(0.11^2, 0.12^2), 2)
        )
      )
    ),
    class = "sglwqs_mids"
  )
}

make_summary_training_mids_fixture <- function() {
  fit1 <- structure(
    list(
      validation = FALSE,
      intercept = 1.10,
      cov_coef = c(age = 0.30, sexM = -0.12)
    ),
    class = "sglwqs"
  )

  fit2 <- structure(
    list(
      validation = FALSE,
      intercept = 0.90,
      cov_coef = c(age = 0.26, sexM = -0.10)
    ),
    class = "sglwqs"
  )

  structure(
    list(
      m = 2,
      n_successful = 2,
      validation = FALSE,
      bootstrap = FALSE,
      family = "gaussian",
      outcome_var = "y",
      exposure_vars = c("x1", "x2"),
      covariate_vars = c("age", "sexM"),
      penalize_covariates = FALSE,
      fits = list(fit1, fit2),
      pooled = list(
        pos_index_sum = 0.82,
        neg_index_sum = 0.24,
        pos_weights = c(x1 = 0.7, x2 = 0.1),
        neg_weights = c(x1 = 0.0, x2 = 0.3),
        mi_selection_freq_pos = c(x1 = 1.0, x2 = 0.5),
        mi_selection_freq_neg = c(x1 = 0.0, x2 = 1.0),
        validation = NULL
      )
    ),
    class = "sglwqs_mids"
  )
}


test_that("validation coefficient table includes Rubin-pooled covariates", {
  fit_mi <- make_summary_validation_mids_fixture()

  coef_df <- .build_mi_validation_coef_table(fit_mi, conf_level = 0.95)

  expect_true(all(c("WQS Positive", "WQS Negative", "age", "sexM") %in% coef_df$Term))
  expect_false("(Intercept)" %in% coef_df$Term)

  age_expected <- rubin_pool(c(0.55, 0.48), c(0.12^2, 0.10^2), 2)
  age_row <- coef_df[coef_df$Term == "age", , drop = FALSE]

  expect_equal(age_row$Estimate, age_expected$estimate)
  expect_equal(age_row$`Std.Error`, age_expected$se)
})


test_that("summary.sglwqs_mids prints Rubin-pooled covariates when validation is available", {
  fit_mi <- make_summary_validation_mids_fixture()

  out <- paste(capture.output(summary(fit_mi)), collapse = "\n")

  expect_match(out, "Rubin-Pooled Validation Coefficients")
  expect_match(out, "age")
  expect_match(out, "sexM")
  expect_match(out, "Penalize covariates:\\s+No")
  expect_match(out, "Rubin's rules are applied term-wise")
})


test_that("summary.sglwqs_mids prints training covariates without validation", {
  fit_mi <- make_summary_training_mids_fixture()

  out <- paste(capture.output(summary(fit_mi)), collapse = "\n")

  expect_match(out, "Pooled Training Coefficients")
  expect_match(out, "age")
  expect_match(out, "sexM")
  expect_match(out, "Standard errors for training coefficients are not Rubin-pooled")
})


test_that("mids variable spec resolver accepts names indices and logical selectors", {
  data_names <- c("y", "x1", "x2", "age")

  expect_identical(.resolve_mids_var_spec(c("x1", "x2"), data_names, "X"), c("x1", "x2"))
  expect_identical(.resolve_mids_var_spec(c(2, 4), data_names, "covariates"), c("x1", "age"))
  expect_identical(.resolve_mids_var_spec(c(FALSE, TRUE, FALSE, TRUE), data_names, "covariates"), c("x1", "age"))
  expect_error(.resolve_single_mids_var(c("y", "x1"), data_names, "y"), "exactly one variable")
})
