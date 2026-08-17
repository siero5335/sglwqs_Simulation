library(testthat)
library(sglwqs)

test_that("compute_diagnostics records weight counts and inference mode", {
  dat <- make_simple_data(n = 120, seed = 900)

  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 1,
    verbose = FALSE
  )

  diag <- fit$diagnostics
  expect_s3_class(diag, "sglwqs_diagnostics")
  expect_true(!is.null(diag$weights_nonzero_count))
  expect_equal(
    diag$weights_nonzero_count$value,
    sum(fit$pos_weights > 0 | fit$neg_weights > 0)
  )
  expect_true(!is.null(diag$inference_mode))
  expect_true(!is.null(diag$two_stage_path))
  expect_true(!is.null(diag$weights_nonzero_fraction))
  expect_true(!is.null(diag$max_abs_exposure_outcome_r))
})

test_that("diagnostics print uses fact-only formatting", {
  dat <- make_simple_data(n = 120, seed = 901)

  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 2,
    verbose = FALSE
  )

  out <- capture.output(print(fit$diagnostics))
  expect_true(any(grepl("^Diagnostics:", out)))
  expect_false(any(grepl("CRITICAL|WARNING|⚠|interpret with caution", out)))
})

test_that("summary_validation and plot_validation_results work for refit = 'full'", {
  dat <- make_simple_data(n = 160, seed = 902)
  cov_df <- data.frame(age = rnorm(nrow(dat)))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    nfolds = 3,
    nlambda = 20,
    seed = 3,
    verbose = FALSE
  ))

  val <- summary_validation(fit)
  inf <- summary_inference(fit)
  plt <- plot_validation_results(fit, include_covariates = TRUE)
  plt2 <- plot_inference_results(fit, include_covariates = TRUE)

  expect_true(is.data.frame(val))
  expect_equal(attr(val, "source"), "refit_info")
  expect_equal(val, inf)
  expect_s3_class(plt, "ggplot")
  expect_s3_class(plt2, "ggplot")
})

test_that("validation_metrics and calibrate use active refit info", {
  dat <- make_binomial_data(n = 180, seed = 903)

  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    family = "binomial",
    refit = "full",
    nfolds = 3,
    nlambda = 20,
    seed = 4,
    verbose = FALSE
  )

  vm <- validation_metrics(fit)
  cal <- calibrate(fit, n_bins = 5)

  expect_s3_class(vm, "sglwqs_validation_metrics")
  expect_true(is.data.frame(cal))
  expect_true("accuracy" %in% names(vm))
  expect_true(nrow(cal) >= 2)
})

test_that("summary_inference and plot_inference_results fall back to bootstrap summaries", {
  dat <- make_simple_data(n = 150, seed = 904)
  cov_df <- data.frame(age = rnorm(nrow(dat)))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
    bootstrap = TRUE,
    n_boot = 8,
    nfolds = 3,
    nlambda = 20,
    seed = 5,
    verbose = FALSE
  ))

  inf <- summary_inference(fit)
  plt <- plot_inference_results(fit)
  plt_no_cov <- plot_inference_results(fit, include_covariates = FALSE)

  expect_s3_class(inf, "sglwqs_validation_summary")
  expect_identical(attr(inf, "source"), "boot_info")
  expect_true(all(c("G1 (positive)", "G2 (negative)") %in% inf$term))
  expect_true("age" %in% inf$term)
  expect_s3_class(plt, "ggplot")
  expect_false("age" %in% plt_no_cov$data$term)
  expect_identical(fit$diagnostics$two_stage_path$value, "bootstrap_only")
})

test_that("summary_inference and plot_inference_results work for mids bootstrap inference", {
  fit_mids <- structure(
    list(
      m = 3,
      n_successful = 3,
      family = "gaussian",
      fits = list(NULL, NULL, NULL),
      pooled = list(
        inference = NULL,
        validation = NULL,
        bootstrap_inference = list(
          source_used = "boot_info",
          has_group_results = TRUE,
          group_results = list(
            G1 = list(
              positive = list(estimate = 0.6, se = 0.15, p_value = 0.01, ci_lower = 0.2, ci_upper = 1.1, df = 12, fmi = 0.3),
              negative = list(estimate = 0.1, se = 0.08, p_value = 0.20, ci_lower = -0.1, ci_upper = 0.3, df = 12, fmi = 0.2)
            ),
            G2 = list(
              positive = list(estimate = 0.2, se = 0.10, p_value = 0.08, ci_lower = 0.01, ci_upper = 0.39, df = 12, fmi = 0.1),
              negative = list(estimate = 0.5, se = 0.12, p_value = 0.02, ci_lower = 0.2, ci_upper = 0.8, df = 12, fmi = 0.4)
            )
          ),
          covariates = list(
            age = list(estimate = 0.3, se = 0.09, p_value = 0.01, ci_lower = 0.05, ci_upper = 0.55, df = 11, fmi = 0.25)
          )
        )
      ),
      groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
      group_by_compound = TRUE,
      validation = FALSE,
      bootstrap = TRUE,
      exposure_vars = c("x1", "x2", "x3", "x4"),
      covariate_vars = "age",
      outcome_var = "y"
    ),
    class = "sglwqs_mids"
  )
  fit_mids$diagnostics <- compute_diagnostics(fit_mids)

  inf <- summary_inference(fit_mids)
  plt <- plot_inference_results(fit_mids)
  plt_no_cov <- plot_inference_results(fit_mids, include_covariates = FALSE)
  out <- capture.output(summary(fit_mids))

  expect_s3_class(inf, "sglwqs_validation_summary")
  expect_identical(attr(inf, "source"), "boot_info")
  expect_true("age" %in% inf$term)
  expect_equal(inf$ci_lower[inf$term == "G1 (positive)"], 0.2)
  expect_equal(inf$ci_upper[inf$term == "age"], 0.55)
  expect_equal(inf$df[inf$term == "G1 (positive)"], 12)
  expect_equal(inf$fmi[inf$term == "age"], 0.25)
  expect_s3_class(plt, "ggplot")
  expect_false("age" %in% plt_no_cov$data$term)
  expect_true(any(grepl("Rubin-Pooled Bootstrap Inference", out, fixed = TRUE)))
  expect_identical(fit_mids$diagnostics$two_stage_path$value, "bootstrap_only")
})

test_that("compute_diagnostics reports all-zero weights and validation set size when applicable", {
  fit_zero <- structure(
    list(
      n = 100L,
      p = 4L,
      var_names = c("x1", "x2", "x3", "x4"),
      cov_names = character(0),
      family = "gaussian",
      pos_weights = c(0, 0, 0, 0),
      neg_weights = c(0, 0, 0, 0),
      pos_coef = c(0, 0, 0, 0),
      neg_coef = c(0, 0, 0, 0),
      fit = list(cvm = c(1, 0.9, 0.85)),
      X_quantile = matrix(rnorm(400), ncol = 4),
      y = rnorm(100),
      validation_info = NULL,
      refit_info = NULL,
      boot_info = NULL,
      cov_matrix = NULL
    ),
    class = "sglwqs"
  )
  diag_zero <- compute_diagnostics(fit_zero)
  expect_true(!is.null(diag_zero$weights_all_zero))
  expect_identical(diag_zero$weights_nonzero_fraction$value, 0)

  dat <- make_simple_data(n = 160, seed = 905)
  fit_val <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    validation = TRUE,
    train_prop = 0.7,
    nfolds = 3,
    nlambda = 20,
    seed = 6,
    verbose = FALSE
  ))
  diag_val <- compute_diagnostics(fit_val)
  expect_true(!is.null(diag_val$validation_set_size))
  expect_identical(diag_val$inference_mode$value, "validation_info")
  expect_identical(diag_val$two_stage_path$value, "validation_info")
})

test_that("compute_diagnostics records pooled inference source for mids refit path", {
  fit_mids <- structure(
    list(
      m = 3,
      n_successful = 3,
      pooled = list(
        inference = list(
          source_used = "refit_info",
          wqs_pos = list(fmi = 0.10),
          wqs_neg = list(fmi = 0.20)
        ),
        validation = NULL,
        bootstrap_inference = NULL
      ),
      fits = list(),
      pos_weights = numeric(),
      neg_weights = numeric()
    ),
    class = "sglwqs_mids"
  )

  diag <- compute_diagnostics(fit_mids)
  expect_identical(diag$two_stage_path$value, "refit_info")
  expect_identical(diag$pooled_inference_source$value, "refit_info")
})
