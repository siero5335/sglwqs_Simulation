library(testthat)

make_tiny_mids_fixture <- function(seed = 123, m = 2, n = 24) {
  set.seed(seed)
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    y = rnorm(n)
  )
  dat$x1[sample(n, 3)] <- NA
  dat$x2[sample(n, 3)] <- NA
  suppressWarnings(mice::mice(dat, m = m, maxit = 1, printFlag = FALSE))
}

test_that("mini MICE fit propagates canonical groups to all fits", {
  skip_if_not_installed("mice")

  groups <- list(G1 = "x1", G2 = "x2", Other = "x3")
  imp <- make_tiny_mids_fixture(seed = 123, m = 2)

  res <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3"),
    outcome_var = "y",
    groups = groups,
    group_by_compound = TRUE,
    nfolds = 3,
    seed = 42,
    verbose = FALSE
  ))

  check_canonical_groups(res, groups)
  invisible(lapply(Filter(Negate(is.null), res$fits), check_canonical_groups, canonical = groups))
})

test_that("mini bootstrap fit keeps group-level index sum names aligned", {
  dat <- make_simple_data(n = 100, seed = 123)
  groups <- list(G1 = "x1", G2 = "x2", Other = c("x3", "x4"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = groups,
    group_by_compound = TRUE,
    bootstrap = TRUE,
    n_boot = 3,
    nfolds = 3,
    seed = 42,
    verbose = FALSE
  ))

  expect_identical(names(fit$pos_index_sum_by_group), names(fit$groups))
  expect_identical(names(fit$neg_index_sum_by_group), names(fit$groups))
})
