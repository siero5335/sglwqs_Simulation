library(testthat)
library(sglwqs)

expect_mode_coverage <- function(fit, boot, mode, grouped = TRUE) {
  expect_true(!is.null(fit$diagnostics))
  expect_true(!is.null(fit$diagnostics$inference_mode))

  if (boot) {
    expect_true(!is.null(fit$boot_info))
    if (grouped) {
      expect_true(!is.null(fit$boot_info$mean_index_sum_by_group_pos))
      expect_true(!is.null(fit$boot_info$mean_index_sum_by_group_neg))
    }
  } else {
    expect_true(is.null(fit$boot_info))
  }

  if (identical(mode, "none") && !boot) {
    expect_error(summary_inference(fit), "No inference results available")
    return(invisible(NULL))
  }

  inf <- summary_inference(fit)
  expect_s3_class(inf, "sglwqs_validation_summary")

  if (identical(mode, "none")) {
    expect_identical(attr(inf, "source"), "boot_info")
  } else if (identical(mode, "full")) {
    expect_identical(attr(inf, "source"), "refit_info")
  } else {
    expect_identical(attr(inf, "source"), "validation_info")
  }
}

expect_mids_mode_coverage <- function(fit, boot, mode, grouped = TRUE) {
  expect_true(!is.null(fit$diagnostics))
  expect_true(!is.null(fit$diagnostics$imputation_count))

  if (boot) {
    expect_true(!is.null(fit$pooled$bootstrap))
    expect_true(!is.null(fit$pooled$bootstrap_inference))
    if (grouped) {
      expect_true(!is.null(fit$pooled$bootstrap_inference$group_results))
    }
  } else {
    expect_true(is.null(fit$pooled$bootstrap_inference))
  }

  if (identical(mode, "none") && !boot) {
    expect_true(is.null(fit$pooled$inference))
    expect_error(summary_inference(fit), "No inference results available")
    return(invisible(NULL))
  }

  inf <- summary_inference(fit)
  expect_s3_class(inf, "sglwqs_validation_summary")

  if (identical(mode, "none")) {
    expect_identical(attr(inf, "source"), "boot_info")
    expect_true(is.null(fit$pooled$inference))
  } else if (identical(mode, "full")) {
    expect_true(!is.null(fit$pooled$inference))
    expect_identical(fit$pooled$inference$source_used, "refit_info")
    expect_true(is.null(fit$pooled$validation))
  } else {
    expect_true(!is.null(fit$pooled$inference))
    expect_true(!is.null(fit$pooled$validation))
    expect_identical(fit$pooled$inference$source_used, "validation_info")
  }
}

test_that("all grouped non-MI mode combinations expose the expected outputs", {
  dat <- make_simple_data(n = 140, seed = 910)
  groups <- list(G1 = c("x1", "x2"), G2 = c("x3", "x4"))
  combos <- expand.grid(
    boot = c(FALSE, TRUE),
    mode = c("none", "full", "validation"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    boot <- combos$boot[[i]]
    mode <- combos$mode[[i]]

    fit <- suppressWarnings(sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      groups = groups,
      bootstrap = boot,
      n_boot = if (boot) 4 else 100,
      refit = if (identical(mode, "validation")) "validation" else mode,
      validation = identical(mode, "validation"),
      train_prop = 0.7,
      nfolds = 3,
      nlambda = 15,
      seed = 10 + i,
      verbose = FALSE
    ))

    expect_mode_coverage(fit, boot = boot, mode = mode, grouped = TRUE)
  }
})

test_that("all grouped MI mode combinations expose the expected outputs", {
  skip_if_not_installed("mice")

  dat <- make_simple_data(n = 120, seed = 911)
  dat$x1[sample(nrow(dat), 12)] <- NA
  dat$x2[sample(nrow(dat), 12)] <- NA
  dat$z1 <- rnorm(nrow(dat))
  dat$z1[sample(nrow(dat), 10)] <- NA

  imp <- mice::mice(dat, m = 2, printFlag = FALSE, seed = 21)
  groups <- list(G1 = c("x1", "x2"), G2 = c("x3", "x4"))
  combos <- expand.grid(
    boot = c(FALSE, TRUE),
    mode = c("none", "full", "validation"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(combos))) {
    boot <- combos$boot[[i]]
    mode <- combos$mode[[i]]

    fit <- suppressWarnings(sglwqs_mice(
      mids_obj = imp,
      exposure_vars = c("x1", "x2", "x3", "x4"),
      outcome_var = "y",
      covariate_vars = "z1",
      groups = groups,
      bootstrap = boot,
      n_boot = if (boot) 4 else 100,
      refit = if (identical(mode, "validation")) "validation" else mode,
      validation = identical(mode, "validation"),
      train_prop = 0.7,
      seed = 20 + i,
      nfolds = 3,
      nlambda = 15,
      verbose = FALSE
    ))

    expect_mids_mode_coverage(fit, boot = boot, mode = mode, grouped = TRUE)
  }
})

test_that("ungrouped bootstrap-only and downstream modes expose overall inference", {
  skip_if_not_installed("mice")

  dat <- make_simple_data(n = 130, seed = 912)

  fit_boot <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    bootstrap = TRUE,
    n_boot = 4,
    nfolds = 3,
    nlambda = 15,
    seed = 41,
    verbose = FALSE
  ))
  inf_boot <- summary_inference(fit_boot)
  expect_true(all(c("WQS Positive loading", "WQS Negative loading") %in% inf_boot$term))
  expect_true(all(is.na(inf_boot$p_value[inf_boot$type == "WQS Loading"])))

  fit_refit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    refit = "full",
    nfolds = 3,
    nlambda = 15,
    seed = 42,
    verbose = FALSE
  ))
  inf_refit <- summary_inference(fit_refit)
  expect_true(all(c("WQS Positive", "WQS Negative") %in% inf_refit$term))

  dat$x1[sample(nrow(dat), 10)] <- NA
  imp <- mice::mice(dat, m = 2, printFlag = FALSE, seed = 31)
  fit_mi_boot <- suppressWarnings(sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c("x1", "x2", "x3", "x4"),
    outcome_var = "y",
    bootstrap = TRUE,
    n_boot = 4,
    seed = 43,
    nfolds = 3,
    nlambda = 15,
    verbose = FALSE
  ))
  inf_mi_boot <- summary_inference(fit_mi_boot)
  expect_true(all(c("WQS Positive loading", "WQS Negative loading") %in% inf_mi_boot$term))
  expect_true(all(is.na(inf_mi_boot$p_value[inf_mi_boot$type == "WQS Loading"])))
})
