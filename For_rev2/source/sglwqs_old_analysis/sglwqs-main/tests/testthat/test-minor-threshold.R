test_that("minor_threshold excludes negligible directions", {
  set.seed(42)
  n <- 2000
  X <- matrix(rnorm(n * 10), n, 10)
  colnames(X) <- paste0("x", 1:10)
  groups <- list(G1 = paste0("x", 1:5), G2 = paste0("x", 6:10))

  # G1: positive only effects; G2: mixed direction effects
  y <- X[, 1] * 0.3 + X[, 2] * 0.2 + X[, 3] * 0.1 -
    X[, 6] * 0.3 + X[, 7] * 0.2 + rnorm(n)

  fit <- suppressWarnings(sglwqs(
    X = X, y = y, groups = groups,
    bootstrap = TRUE, n_boot = 50,
    validation = TRUE, train_prop = 0.6,
    minor_threshold = 0.10, seed = 42,
    verbose = FALSE
  ))

  # G1 negative should be excluded (positive-only signal)
  expect_true("G1" %in% names(fit$excluded_directions))
  expect_equal(fit$excluded_directions[["G1"]], "negative")

  # Excluded direction should have NA in validation results
  g1_res <- fit$validation_info$group_results$G1
  expect_true(is.na(g1_res$neg_estimate))
  expect_true(is.na(g1_res$neg_se))
  expect_true(is.na(g1_res$neg_pvalue))
  expect_true(isTRUE(g1_res$neg_excluded))

  # Non-excluded direction should have values
  expect_false(isTRUE(g1_res$pos_excluded))

  # minor_threshold is stored in output
  expect_equal(fit$minor_threshold, 0.10)
})

test_that("minor_threshold = 0 disables exclusion", {
  set.seed(42)
  n <- 1000
  X <- matrix(rnorm(n * 10), n, 10)
  colnames(X) <- paste0("x", 1:10)
  groups <- list(G1 = paste0("x", 1:5), G2 = paste0("x", 6:10))
  y <- X[, 1] * 0.3 + X[, 2] * 0.2 + rnorm(n)

  fit <- suppressWarnings(sglwqs(
    X = X, y = y, groups = groups,
    bootstrap = TRUE, n_boot = 50,
    validation = TRUE, train_prop = 0.6,
    minor_threshold = 0, seed = 42,
    verbose = FALSE
  ))

  expect_length(fit$excluded_directions, 0)
})

test_that("minor_threshold has no effect without groups", {
  set.seed(42)
  n <- 500
  X <- matrix(rnorm(n * 5), n, 5)
  colnames(X) <- paste0("x", 1:5)
  y <- X[, 1] * 0.3 + rnorm(n)

  fit <- sglwqs(
    X = X, y = y,
    validation = TRUE, train_prop = 0.6,
    minor_threshold = 0.10, seed = 42,
    verbose = FALSE
  )

  expect_length(fit$excluded_directions, 0)
})

test_that("minor_threshold has no effect without refit", {
  set.seed(42)
  n <- 500
  X <- matrix(rnorm(n * 10), n, 10)
  colnames(X) <- paste0("x", 1:10)
  groups <- list(G1 = paste0("x", 1:5), G2 = paste0("x", 6:10))
  y <- X[, 1] * 0.3 + rnorm(n)

  fit <- sglwqs(
    X = X, y = y, groups = groups,
    minor_threshold = 0.10, seed = 42,
    verbose = FALSE
  )

  expect_length(fit$excluded_directions, 0)
})

test_that("minor_threshold also applies to full refit", {
  set.seed(42)
  n <- 1500
  X <- matrix(rnorm(n * 10), n, 10)
  colnames(X) <- paste0("x", 1:10)
  groups <- list(G1 = paste0("x", 1:5), G2 = paste0("x", 6:10))
  y <- X[, 1] * 0.35 + X[, 2] * 0.2 + X[, 3] * 0.1 -
    X[, 6] * 0.25 + X[, 7] * 0.2 + rnorm(n)

  fit <- suppressWarnings(sglwqs(
    X = X, y = y, groups = groups,
    refit = "full",
    minor_threshold = 0.10, seed = 42,
    verbose = FALSE
  ))

  expect_true("G1" %in% names(fit$excluded_directions))
  expect_equal(fit$excluded_directions[["G1"]], "negative")

  g1_res <- fit$refit_info$group_results$G1
  expect_true(isTRUE(g1_res$neg_excluded))
  expect_true(is.na(g1_res$neg_estimate))
})

test_that("minor_threshold propagates through sglwqs_mice", {
  skip_if_not_installed("mice")

  set.seed(42)
  n <- 500
  dat <- data.frame(matrix(rnorm(n * 10), nrow = n))
  names(dat) <- paste0("x", 1:10)
  dat$y <- dat$x1 * 0.35 + dat$x2 * 0.2 + dat$x3 * 0.1 -
    dat$x6 * 0.25 + dat$x7 * 0.2 + rnorm(n)
  dat$x1[sample.int(n, 25)] <- NA_real_

  groups <- list(G1 = paste0("x", 1:5), G2 = paste0("x", 6:10))
  imp <- suppressWarnings(mice::mice(dat, m = 2, maxit = 1, printFlag = FALSE))

  fit_mi <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = paste0("x", 1:10),
    outcome_var = "y",
    groups = groups,
    validation = TRUE,
    train_prop = 0.6,
    minor_threshold = 0.10,
    seed = 42,
    verbose = FALSE
  ))

  expect_equal(fit_mi$minor_threshold, 0.10)
  expect_true(all(vapply(
    Filter(Negate(is.null), fit_mi$fits),
    function(f) is.list(f$excluded_directions),
    logical(1)
  )))
})

test_that("minor_threshold rejects invalid inputs", {
  set.seed(42)
  n <- 100
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- paste0("x", 1:4)
  y <- rnorm(n)

  expect_error(sglwqs(X = X, y = y, minor_threshold = -0.1),
               "minor_threshold")
  expect_error(sglwqs(X = X, y = y, minor_threshold = 1.5),
               "minor_threshold")
  expect_error(sglwqs(X = X, y = y, minor_threshold = "a"),
               "minor_threshold")
  expect_error(sglwqs(X = X, y = y, minor_threshold = c(0.1, 0.2)),
               "minor_threshold")
  expect_error(sglwqs(X = X, y = y, minor_threshold = NA_real_),
               "minor_threshold")
  expect_error(sglwqs(X = X, y = y, minor_threshold = Inf),
               "minor_threshold")
})
