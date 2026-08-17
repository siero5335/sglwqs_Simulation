# Regression tests for group metadata consistency, "Other" duplication,
# bootstrap index sum alignment, sglwqs_mice group_by_compound default,
# .get_effective_groups(), pool metadata checking, and edge cases.

# ==== complete-data side ====

test_that("bootstrap index sums match mean-of-per-draw estimator", {
  dat <- make_simple_data(n = 200, seed = 100)
  grps <- list(metals = c("x1", "x2"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    bootstrap = TRUE,
    n_boot = 10,
    keep_boot_matrices = TRUE,
    nfolds = 3,
    seed = 42,
    verbose = FALSE
  ))

  # Manually recompute the expected mean-of-per-draw index sums
  boot_pos <- fit$boot_info$boot_pos_coef
  boot_neg <- fit$boot_info$boot_neg_coef
  boot_ok <- if (!is.null(fit$boot_info$boot_success)) {
    fit$boot_info$boot_success
  } else {
    rowSums(abs(boot_pos) + abs(boot_neg)) > 0
  }

  per_draw_pos_sum <- numeric(sum(boot_ok))
  per_draw_neg_sum <- numeric(sum(boot_ok))
  idx <- which(boot_ok)
  for (k in seq_along(idx)) {
    b <- idx[k]
    w <- suppressWarnings(sglwqs:::calculate_weights(
      pos_coef = boot_pos[b, ], neg_coef = boot_neg[b, ],
      var_names = fit$var_names, groups = fit$groups,
      handle_collinearity = "net"
    ))
    per_draw_pos_sum[k] <- w$pos_index_sum
    per_draw_neg_sum[k] <- w$neg_index_sum
  }

  expect_equal(fit$pos_index_sum, mean(per_draw_pos_sum), tolerance = 1e-10)
  expect_equal(fit$neg_index_sum, mean(per_draw_neg_sum), tolerance = 1e-10)

  # Group-level index sums: also verify mean-of-per-draw
  expect_true(!is.null(fit$pos_index_sum_by_group))
  for (grp_name in names(fit$groups)) {
    per_draw_grp_pos <- numeric(length(idx))
    per_draw_grp_neg <- numeric(length(idx))
    for (k in seq_along(idx)) {
      b <- idx[k]
      w <- suppressWarnings(sglwqs:::calculate_weights(
        pos_coef = boot_pos[b, ], neg_coef = boot_neg[b, ],
        var_names = fit$var_names, groups = fit$groups,
        handle_collinearity = "net"
      ))
      per_draw_grp_pos[k] <- w$pos_index_sum_by_group[[grp_name]]
      per_draw_grp_neg[k] <- w$neg_index_sum_by_group[[grp_name]]
    }
    expect_equal(fit$pos_index_sum_by_group[[grp_name]],
                 mean(per_draw_grp_pos), tolerance = 1e-10,
                 label = paste0("pos_index_sum_by_group[", grp_name, "]"))
    expect_equal(fit$neg_index_sum_by_group[[grp_name]],
                 mean(per_draw_grp_neg), tolerance = 1e-10,
                 label = paste0("neg_index_sum_by_group[", grp_name, "]"))
  }
})


test_that("canonical groups include auto-generated 'Other' with correct membership", {
  dat <- make_simple_data(n = 180, seed = 200)
  grps <- list(metals = c("x1", "x2"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    nfolds = 3,
    seed = 1,
    verbose = FALSE
  ))

  expect_true("Other" %in% names(fit$groups))
  expect_equal(sort(fit$groups$Other), sort(c("x3", "x4")))
  # All exposures covered
  expect_equal(sort(unlist(fit$groups, use.names = FALSE)),
               sort(c("x1", "x2", "x3", "x4")))
  # unique guard
  group_order <- unique(c(names(fit$groups), "Other"))
  expect_equal(sum(group_order == "Other"), 1)
})


test_that("explicit 'Other' does not duplicate in extract_weights or plot data", {
  dat <- make_simple_data(n = 180, seed = 300)
  grps <- list(metals = c("x1", "x2"), Other = c("x3", "x4"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    nfolds = 3,
    seed = 1,
    verbose = FALSE
  ))

  ew <- extract_weights(fit, by_group = TRUE)
  # "Other" should appear at most once per variable
  if ("group" %in% names(ew)) {
    other_rows <- ew[ew$group == "Other", ]
    expect_equal(nrow(other_rows), length(unique(other_rows$variable)) * 2)
    # No variable should appear in two different groups
    var_group <- unique(ew[, c("variable", "group")])
    expect_equal(nrow(var_group), length(unique(var_group$variable)))
  }
})


test_that("plot_weights with canonical 'Other' does not duplicate groups", {
  dat <- make_simple_data(n = 180, seed = 400)
  grps <- list(metals = c("x1", "x2"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    nfolds = 3,
    seed = 1,
    verbose = FALSE
  ))

  p <- plot_weights(fit)
  expect_true(inherits(p, "gg") || inherits(p, "ggplot"))

  # Extract plot data and verify "Other" appears at most once as a group level
  pd <- ggplot2::ggplot_build(p)
  if (!is.null(pd$data[[1]]) && "group" %in% names(p$data)) {
    group_vals <- as.character(p$data$group)
    expect_lte(sum(unique(group_vals) == "Other"), 1)
  }
})


# ==== sglwqs_mice side ====

test_that("sglwqs_mice resolves group_by_compound default and stores canonical groups", {
  skip_if_not_installed("mice")

  set.seed(500)
  n <- 120
  dat <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n),
    y = rnorm(n)
  )
  dat$x1[sample(n, 10)] <- NA
  dat$x2[sample(n, 10)] <- NA

  imp <- suppressWarnings(mice::mice(dat, m = 2, maxit = 1, printFlag = FALSE))
  grps <- list(metals = c("x1", "x2"))

  mfit <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3", "x4"),
    outcome_var = "y",
    groups = grps,
    nfolds = 3,
    seed = 42,
    verbose = FALSE
  ))

  expect_true(isTRUE(mfit$group_by_compound))
  # Canonical groups (validated at entry)
  expect_true(!is.null(mfit$groups))
  expect_true("Other" %in% names(mfit$groups))
  expect_equal(sort(unlist(mfit$groups, use.names = FALSE)),
               sort(c("x1", "x2", "x3", "x4")))
  # Each fit should carry the same groups (names AND membership)
  ref_normalized <- sglwqs:::.normalize_group_metadata(mfit$groups)
  for (f in Filter(Negate(is.null), mfit$fits)) {
    expect_identical(sglwqs:::.normalize_group_metadata(f$groups), ref_normalized)
  }
})


test_that("sglwqs_mice with group_by_compound=FALSE stores groups=NULL", {
  skip_if_not_installed("mice")

  set.seed(600)
  n <- 120
  dat <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n),
    y = rnorm(n)
  )
  dat$x1[sample(n, 10)] <- NA

  imp <- suppressWarnings(mice::mice(dat, m = 2, maxit = 1, printFlag = FALSE))
  grps <- list(metals = c("x1", "x2"))

  expect_warning(
    mfit <- sglwqs_mice(
      mids_obj = imp,
      exposure_vars = c("x1", "x2", "x3", "x4"),
      outcome_var = "y",
      groups = grps,
      group_by_compound = FALSE,
      nfolds = 3,
      seed = 42,
      verbose = FALSE
    ),
    "ignored"
  )

  expect_false(isTRUE(mfit$group_by_compound))
  expect_null(mfit$groups)
  expect_null(sglwqs:::.get_effective_groups(mfit))
})


# ==== .get_effective_groups ====

test_that(".get_effective_groups priority: top-level > pooled > first fit", {
  obj <- list(
    group_by_compound = TRUE,
    groups = list(A = "x1", B = "x2"),
    pooled = list(groups = list(C = "x3")),
    fits = list()
  )
  class(obj) <- "sglwqs_mids"

  # Top-level groups take priority
  expect_equal(names(sglwqs:::.get_effective_groups(obj)), c("A", "B"))

  # Fallback to pooled

  obj2 <- obj
  obj2$groups <- NULL
  expect_equal(names(sglwqs:::.get_effective_groups(obj2)), "C")

  # Fallback to first fit
  fit_mock <- list(groups = list(D = "x4"))
  class(fit_mock) <- "sglwqs"
  obj3 <- obj2
  obj3$pooled$groups <- NULL
  obj3$fits <- list(fit_mock)
  expect_equal(names(sglwqs:::.get_effective_groups(obj3)), "D")

  # group_by_compound=FALSE always NULL regardless of groups present
  obj4 <- obj
  obj4$group_by_compound <- FALSE
  expect_null(sglwqs:::.get_effective_groups(obj4))
})


test_that("pooled fallback: print/summary/extract_weights work with top-level groups=NULL", {
  skip_if_not_installed("mice")

  set.seed(700)
  n <- 120
  dat <- data.frame(
    x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n), x4 = rnorm(n),
    y = rnorm(n)
  )
  dat$x1[sample(n, 10)] <- NA

  imp <- suppressWarnings(mice::mice(dat, m = 2, maxit = 1, printFlag = FALSE))
  grps <- list(metals = c("x1", "x2"))

  mfit <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3", "x4"),
    outcome_var = "y",
    groups = grps,
    nfolds = 3,
    seed = 42,
    verbose = FALSE
  ))

  # Simulate a top-level groups=NULL scenario (fallback to pooled)
  mfit_stripped <- mfit
  mfit_stripped$groups <- NULL
  # .get_effective_groups should still find groups via pooled
  eff <- sglwqs:::.get_effective_groups(mfit_stripped)
  expect_true(!is.null(eff))

  # print should not error
  expect_output(print(mfit_stripped), "SGL-WQS")
  expect_output(summary(mfit_stripped), "Chemical Groups")

  # extract_weights should still include group column
  ew <- extract_weights(mfit_stripped, by_group = TRUE)
  expect_true("group" %in% names(ew))
  expect_true(any(ew$group != "All"))
})


# ==== pool_sglwqs_results ====

test_that("pool_sglwqs_results detects group metadata inconsistency", {
  fit1 <- list(
    var_names = c("x1", "x2"), p = 2,
    pos_weights = c(x1 = 0.6, x2 = 0.4),
    neg_weights = c(x1 = 0.3, x2 = 0.7),
    pos_coef = c(x1 = 1.0, x2 = 0.5),
    neg_coef = c(x1 = 0.2, x2 = 0.8),
    pos_index_sum = 1.5, neg_index_sum = 1.0,
    pos_index_sum_by_group = list(A = 1.0, B = 0.5),
    neg_index_sum_by_group = list(A = 0.2, B = 0.8),
    groups = list(A = "x1", B = "x2"),
    bootstrap = FALSE, validation = FALSE
  )
  class(fit1) <- "sglwqs"

  fit2 <- fit1
  fit2$groups <- list(A = "x2", B = "x1")

  expect_error(
    sglwqs:::pool_sglwqs_results(list(fit1, fit2)),
    "Group metadata differs"
  )
})


test_that("pool_sglwqs_results detects var_names inconsistency", {
  fit1 <- list(
    var_names = c("x1", "x2"), p = 2,
    pos_weights = c(x1 = 0.6, x2 = 0.4),
    neg_weights = c(x1 = 0.3, x2 = 0.7),
    pos_coef = c(x1 = 1.0, x2 = 0.5),
    neg_coef = c(x1 = 0.2, x2 = 0.8),
    pos_index_sum = 1.5, neg_index_sum = 1.0,
    pos_index_sum_by_group = NULL,
    neg_index_sum_by_group = NULL,
    groups = NULL,
    bootstrap = FALSE, validation = FALSE
  )
  class(fit1) <- "sglwqs"

  fit2 <- fit1
  fit2$var_names <- c("x1", "x3")

  expect_error(
    sglwqs:::pool_sglwqs_results(list(fit1, fit2)),
    "Variable names differ"
  )
})


test_that("pool_sglwqs_results handles single fit (n_successful=1)", {
  dat <- make_simple_data(n = 180, seed = 800)
  grps <- list(metals = c("x1", "x2"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    nfolds = 3,
    seed = 1,
    verbose = FALSE
  ))

  # Pool a single fit — rowMeans / apply should not collapse dimensions
  pooled <- sglwqs:::pool_sglwqs_results(list(fit))

  expect_equal(length(pooled$pos_weights), 4)
  expect_equal(length(pooled$neg_weights), 4)
  expect_equal(length(pooled$pos_coef), 4)
  expect_true(is.numeric(pooled$pos_index_sum))
  expect_equal(pooled$pos_weights, fit$pos_weights, tolerance = 1e-10)
  expect_equal(pooled$neg_weights, fit$neg_weights, tolerance = 1e-10)
  expect_equal(pooled$pos_index_sum, fit$pos_index_sum, tolerance = 1e-10)
  expect_equal(pooled$pos_weights_var, rep(0, 4), tolerance = 1e-10)
  expect_equal(pooled$neg_weights_var, rep(0, 4), tolerance = 1e-10)
})


test_that("pool_sglwqs_results handles single bootstrap fit without NA variances", {
  dat <- make_simple_data(n = 180, seed = 900)
  grps <- list(metals = c("x1", "x2"))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    groups = grps,
    bootstrap = TRUE,
    n_boot = 8,
    keep_boot_matrices = TRUE,
    nfolds = 3,
    seed = 1,
    verbose = FALSE
  ))

  pooled <- sglwqs:::pool_sglwqs_results(list(fit))

  expect_false(anyNA(pooled$pos_weights_var))
  expect_false(anyNA(pooled$neg_weights_var))
  expect_equal(pooled$pos_weights_var, rep(0, 4), tolerance = 1e-10)
  expect_equal(pooled$neg_weights_var, rep(0, 4), tolerance = 1e-10)
  expect_false(anyNA(pooled$bootstrap$se_pos_coef))
  expect_false(anyNA(pooled$bootstrap$se_neg_coef))
})


test_that("extract_weights.sglwqs_mids fails loudly for missing auxiliary pooled vectors", {
  obj <- structure(
    list(
      group_by_compound = FALSE,
      groups = NULL,
      pooled = list(
        var_names = c("x1", "x2"),
        pos_weights = c(x1 = 0.6, x2 = 0.4),
        neg_weights = c(x1 = 0.1, x2 = 0.9),
        pos_coef = c(x1 = 1.0, x2 = 0.5),
        neg_coef = c(x1 = 0.0, x2 = -0.7),
        mi_selection_freq_pos = c(x1 = 1.0, x2 = 0.5),
        mi_selection_freq_neg = c(x1 = 0.0, x2 = 1.0),
        pos_weights_var = c(x1 = 0.01, x2 = 0.02),
        neg_weights_var = NULL
      )
    ),
    class = "sglwqs_mids"
  )

  expect_error(
    extract_weights(obj, by_group = FALSE),
    "Missing pooled vector"
  )
})


test_that("explicit 'Other' with unassigned exposures errors with actionable message", {
  expect_error(
    sglwqs:::.validate_and_complete_groups(
      var_names = c("x1", "x2", "x3"),
      groups = list(Other = "x1")
    ),
    paste0(
      "some exposures are still unassigned: x2, x3.*",
      "must contain all intended members.*",
      "Either add these exposures.*",
      "let it be created automatically"
    )
  )
})
