test_that("sglwqs_mids selection and bootstrap plots return ggplot objects", {
  fit_mice <- structure(
    list(
      m = 3,
      n_successful = 3,
      family = "gaussian",
      bootstrap = TRUE,
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
        neg_weights_var = c(x1 = 0, x2 = 0, x3 = 0),
        bootstrap = list(
          selection_freq_pos = c(x1 = 0.75, x2 = 0.35, x3 = 0.00),
          selection_freq_neg = c(x1 = 0.00, x2 = 0.45, x3 = 0.65),
          se_pos_coef = c(x1 = 0.12, x2 = 0.09, x3 = 0.00),
          se_neg_coef = c(x1 = 0.00, x2 = 0.10, x3 = 0.11)
        )
      )
    ),
    class = "sglwqs_mids"
  )

  p_sel <- plot_selection_frequency(fit_mice, source = "mi")
  p_boot <- plot_bootstrap_ci(fit_mice)

  expect_s3_class(p_sel, "ggplot")
  expect_s3_class(p_boot, "ggplot")
  expect_true(all(c("x1", "x2", "x3") %in% as.character(p_sel$data$variable)))
})


test_that("sglwqs_mids pooled CV plot interpolates successful imputations", {
  fit1 <- structure(
    list(
      fit = list(
        lambda = exp(seq(0, -3, length.out = 5)),
        cvm = c(1.5, 1.2, 1.0, 1.1, 1.3),
        cvsd = c(0.10, 0.09, 0.08, 0.09, 0.11),
        lambda.min = exp(-1.5),
        lambda.1se = exp(-0.75)
      )
    ),
    class = "sglwqs"
  )
  fit2 <- structure(
    list(
      fit = list(
        lambda = exp(seq(-0.1, -3.1, length.out = 6)),
        cvm = c(1.6, 1.25, 1.05, 1.02, 1.12, 1.35),
        cvsd = c(0.11, 0.10, 0.08, 0.08, 0.09, 0.12),
        lambda.min = exp(-1.9),
        lambda.1se = exp(-1.0)
      )
    ),
    class = "sglwqs"
  )

  fit_mice <- structure(
    list(
      n_successful = 2,
      groups = NULL,
      pooled = list(
        var_names = c("x1", "x2", "x3"),
        pos_weights = c(x1 = 0.3, x2 = 0.0, x3 = 0.2),
        neg_weights = c(x1 = 0.0, x2 = 0.4, x3 = 0.0)
      ),
      fits = list(fit1, fit2)
    ),
    class = "sglwqs_mids"
  )

  p <- plot_cv(fit_mice)

  expect_s3_class(p, "ggplot")
  expect_true(all(c("lambda", "cvm", "cvsd") %in% names(p$data)))
  expect_gt(nrow(p$data), 1)
})


test_that("sglwqs_mids calibration and autoplot methods work with pooled validation predictions", {
  set.seed(123)
  dat1 <- data.frame(
    y = rbinom(40, 1, 0.5),
    x = rnorm(40)
  )
  dat2 <- data.frame(
    y = rbinom(40, 1, 0.5),
    x = rnorm(40)
  )

  glm1 <- glm(y ~ x, data = dat1, family = binomial())
  glm2 <- glm(y ~ x, data = dat2, family = binomial())

  fit1 <- structure(
    list(
      family = "binomial",
      validation = TRUE,
      validation_info = list(glm_fit = glm1)
    ),
    class = "sglwqs"
  )
  fit2 <- structure(
    list(
      family = "binomial",
      validation = TRUE,
      validation_info = list(glm_fit = glm2)
    ),
    class = "sglwqs"
  )

  fit_mice <- structure(
    list(
      n_successful = 2,
      family = "binomial",
      validation = TRUE,
      fits = list(fit1, fit2),
      pooled = list(
        pos_weights = c(x1 = 0.4),
        neg_weights = c(x1 = 0.0),
        var_names = c("x1")
      )
    ),
    class = "sglwqs_mids"
  )

  cal <- calibrate(fit_mice, n_bins = 5)
  p <- ggplot2::autoplot(fit_mice, type = "calibration", n_bins = 5)

  expect_true(is.data.frame(cal))
  expect_true(all(c("bin", "pred_mean", "obs_rate") %in% names(cal)))
  expect_s3_class(p, "ggplot")
})
