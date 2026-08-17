library(testthat)
library(sglwqs)

make_mice_refit_pool_fixture <- function(seed = 42, n = 120, m = 2) {
  set.seed(seed)
  X <- matrix(runif(n * 6), n, 6)
  colnames(X) <- paste0("x", 1:6)
  Z <- data.frame(
    z1 = rnorm(n),
    z2 = rbinom(n, 1, 0.5)
  )

  X_q <- sglwqs:::quantile_transform(X, n_quantiles = 4)$q
  beta_x <- c(0.4, 0.3, 0, 0, -0.4, 0)
  beta_z <- c(0.5, 0.6)
  y <- as.numeric(X_q %*% beta_x + as.matrix(Z) %*% beta_z + rnorm(n, sd = 0.7))

  dat <- data.frame(y = y, X, Z)
  dat$x1[sample(n, 12)] <- NA
  dat$z1[sample(n, 8)] <- NA

  suppressWarnings(mice::mice(dat, m = m, maxit = 1, printFlag = FALSE, seed = seed))
}

test_that("sglwqs_mice pools refit_info when refit = 'full'", {
  skip_if_not_installed("mice")

  imp <- make_mice_refit_pool_fixture(seed = 42, n = 120, m = 2)
  fit <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = paste0("x", 1:6),
    outcome_var = "y",
    covariate_vars = c("z1", "z2"),
    groups = list(g1 = paste0("x", 1:3), g2 = paste0("x", 4:6)),
    family = "gaussian",
    validation = FALSE,
    refit = "full",
    nfolds = 3,
    nlambda = 20,
    seed = 1,
    verbose = FALSE
  ))

  expect_false(is.null(fit$pooled$inference))
  expect_identical(fit$pooled$inference$source_used, "refit_info")
  expect_false(is.null(fit$pooled$covariates))
  expect_true("z1" %in% names(fit$pooled$covariates))
  expect_true("z2" %in% names(fit$pooled$covariates))
  expect_true(is.finite(fit$pooled$covariates$z1$estimate))

  coef_df <- .build_mi_validation_coef_table(fit, conf_level = 0.95)
  expect_true(all(c("g1 (positive)", "g1 (negative)", "z1", "z2") %in% coef_df$Term))
  expect_false("(Intercept)" %in% coef_df$Term)

  out <- paste(capture.output(summary(fit)), collapse = "\n")
  expect_match(out, "Rubin-Pooled Refit Coefficients")
  expect_match(out, "z1")
  expect_match(out, "z2")
})

test_that("validation-based MICE pooling still populates both inference and validation alias", {
  skip_if_not_installed("mice")

  imp <- make_mice_refit_pool_fixture(seed = 99, n = 100, m = 2)
  fit <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = paste0("x", 1:6),
    outcome_var = "y",
    covariate_vars = c("z1", "z2"),
    groups = list(g1 = paste0("x", 1:3), g2 = paste0("x", 4:6)),
    family = "gaussian",
    validation = TRUE,
    train_prop = 0.7,
    nfolds = 3,
    nlambda = 20,
    seed = 2,
    verbose = FALSE
  ))

  expect_false(is.null(fit$pooled$inference))
  expect_false(is.null(fit$pooled$validation))
  expect_identical(fit$pooled$inference$source_used, "validation_info")
  expect_equal(fit$pooled$inference, fit$pooled$validation)
})

test_that("sglwqs_mice accepts formal refit argument", {
  skip_if_not_installed("mice")

  imp <- make_mice_refit_pool_fixture(seed = 77, n = 90, m = 2)
  fit <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = paste0("x", 1:6),
    outcome_var = "y",
    refit = "full",
    nfolds = 3,
    nlambda = 15,
    seed = 3,
    verbose = FALSE
  ))

  expect_true(!is.null(fit$pooled$inference))
  expect_identical(fit$pooled$inference$source_used, "refit_info")
})
