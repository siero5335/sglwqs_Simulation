library(testthat)

make_minimal_mice_data <- function(seed = 123, n = 40) {
  set.seed(seed)
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    y = rnorm(n)
  )
  dat$x1[sample(n, 4)] <- NA
  dat$x2[sample(n, 4)] <- NA
  dat
}

make_minimal_mids <- function(seed = 123, m = 2) {
  dat <- make_minimal_mice_data(seed = seed)
  suppressWarnings(mice::mice(dat, m = m, maxit = 1, printFlag = FALSE))
}

test_that("explicit 'Other' triggers stop if unassigned variables remain", {
  dat <- make_minimal_mice_data(seed = 1, n = 30)
  broken_groups <- list(G1 = "x1", G2 = "x2", Other = character(0))

  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3")],
      y = dat$y,
      groups = broken_groups,
      group_by_compound = TRUE,
      nfolds = 3,
      verbose = FALSE
    ),
    "unassigned"
  )
})

test_that("sglwqs_mice propagates canonical groups to fits and top-level object", {
  skip_if_not_installed("mice")

  imp <- make_minimal_mids(seed = 2, m = 2)
  groups <- list(G1 = "x1", G2 = "x2", Other = "x3")

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

  top_groups <- res$groups
  check_canonical_groups(top_groups, groups)

  fits_groups <- lapply(Filter(Negate(is.null), res$fits), function(f) f$groups)
  expect_true(length(fits_groups) > 0)
  invisible(lapply(fits_groups, check_canonical_groups, canonical = groups))
})

test_that("plot/order helpers handle explicit 'Other' without duplication", {
  skip_if_not_installed("mice")

  imp <- make_minimal_mids(seed = 3, m = 2)
  groups <- list(G1 = "x1", G2 = "x2", Other = "x3")

  res <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3"),
    outcome_var = "y",
    groups = groups,
    group_by_compound = TRUE,
    nfolds = 3,
    seed = 99,
    verbose = FALSE
  ))

  g_eff <- sglwqs:::.get_effective_groups(res)
  group_order <- sglwqs:::.group_display_order(g_eff, c(names(g_eff), "Other"))
  expect_equal(sum(group_order == "Other"), 1)
})

test_that("single successful fit handles pooling", {
  skip_if_not_installed("mice")

  imp <- make_minimal_mids(seed = 4, m = 1)
  groups <- list(G1 = "x1", G2 = "x2", Other = "x3")

  res_single <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3"),
    outcome_var = "y",
    groups = groups,
    group_by_compound = TRUE,
    nfolds = 3,
    seed = 7,
    verbose = FALSE
  ))

  pooled <- sglwqs:::pool_sglwqs_results(Filter(Negate(is.null), res_single$fits))
  expect_false(is.null(pooled$pos_weights))
  expect_false(is.null(pooled$neg_weights))
})

test_that("bootstrap index sums preserve group names", {
  dat <- make_simple_data(n = 120, seed = 5)
  groups <- list(G1 = "x1", G2 = "x2", Other = c("x3", "x4"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = groups,
    group_by_compound = TRUE,
    bootstrap = TRUE,
    n_boot = 5,
    nfolds = 3,
    seed = 11,
    verbose = FALSE
  ))

  expect_identical(names(fit$pos_index_sum_by_group), names(fit$groups))
  expect_identical(names(fit$neg_index_sum_by_group), names(fit$groups))
})
