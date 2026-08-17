test_that("sglwqs returns expected object structure", {
  dat <- make_simple_data()
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 1,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "sglwqs")
  expect_true(is.numeric(fit$pos_weights))
  expect_true(is.numeric(fit$neg_weights))
  expect_equal(length(fit$var_names), 4)
})

test_that("sglwqs validates train_prop early", {
  dat <- make_simple_data(seed = 1001)

  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      validation = TRUE,
      train_prop = 1,
      nfolds = 3,
      nlambda = 10,
      verbose = FALSE
    ),
    "strictly between 0 and 1",
    fixed = TRUE
  )

  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      train_prop = 0,
      nfolds = 3,
      nlambda = 10,
      verbose = FALSE
    ),
    "strictly between 0 and 1",
    fixed = TRUE
  )
})


test_that("weights stay in [0,1] and sum to 1 when selected", {
  dat <- make_simple_data(seed = 456)
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 2,
    verbose = FALSE
  )
  
  expect_true(all(fit$pos_weights >= 0 & fit$pos_weights <= 1))
  expect_true(all(fit$neg_weights >= 0 & fit$neg_weights <= 1))
  
  if (sum(fit$pos_weights) > 0) {
    expect_equal(sum(fit$pos_weights), 1, tolerance = 1e-6)
  }
  if (sum(fit$neg_weights) > 0) {
    expect_equal(sum(fit$neg_weights), 1, tolerance = 1e-6)
  }
})


test_that("known signal direction is detected", {
  dat <- make_simple_data(seed = 789)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 3,
    verbose = FALSE
  ))
  
  expect_gt(fit$pos_weights["x1"], fit$neg_weights["x1"])
  expect_gt(fit$neg_weights["x2"], fit$pos_weights["x2"])
})


test_that("predict and coef methods work on new data", {
  dat <- make_simple_data(seed = 111)
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 4,
    verbose = FALSE
  )
  
  idx <- predict(fit, newdata = dat[1:10, c("x1", "x2", "x3", "x4")], type = "index")
  pred <- predict(fit, newdata = dat[1:10, c("x1", "x2", "x3", "x4")], type = "response")
  b <- coef(fit, what = "coefficients", direction = "both")
  w <- coef(fit, what = "weights", direction = "both")
  
  expect_equal(nrow(idx), 10)
  expect_true(all(c("wqs_pos", "wqs_neg") %in% names(idx)))
  expect_length(pred, 10)
  expect_true(is.numeric(b))
  expect_true(is.numeric(w))
  expect_true(any(grepl("^pos_b\\.", names(b))))
  expect_true(any(grepl("^neg_w\\.", names(w))))
})


test_that("bootstrap checkpoint retries failed iterations on resume", {
  checkpoint_dir <- tempfile("sglwqs-checkpoint-")
  dir.create(checkpoint_dir)
  
  X_quantile <- matrix(rnorm(24), ncol = 2)
  colnames(X_quantile) <- c("x1", "x2")
  y <- rnorm(12)
  fail_once <- local({
    attempted <- FALSE
    function(...) {
      if (!attempted) {
        attempted <<- TRUE
        stop("forced bootstrap failure")
      }
      list(pos_coef = c(0.5, 0), neg_coef = c(0, 0.25))
    }
  })
  
  suppressWarnings(
    expect_error(
      testthat::with_mocked_bindings(
        sglwqs:::bootstrap_sgl(
          X_quantile = X_quantile,
          y = y,
          cov_matrix = NULL,
          var_names = colnames(X_quantile),
          cov_names = NULL,
          groups = NULL,
          group_by_compound = FALSE,
          group_structure = "direction",
          penalize_covariates = FALSE,
          family = "gaussian",
          lambda = "lambda.min",
          nfolds = 2,
          n_boot = 1,
          seed = 999,
          verbose = FALSE,
          parallel = FALSE,
          checkpoint_dir = checkpoint_dir,
          checkpoint_interval = 1,
          cleanup_checkpoint = FALSE
        ),
        fit_sgl_core = fail_once,
        .package = "sglwqs"
      ),
      "All bootstrap iterations failed"
    )
  )
  
  checkpoint_file <- file.path(checkpoint_dir, "bootstrap_checkpoint.rds")
  checkpoint_data <- readRDS(checkpoint_file)
  expect_length(checkpoint_data$completed_boots, 0)
  
  resumed <- testthat::with_mocked_bindings(
    sglwqs:::bootstrap_sgl(
      X_quantile = X_quantile,
      y = y,
      cov_matrix = NULL,
      var_names = colnames(X_quantile),
      cov_names = NULL,
      groups = NULL,
      group_by_compound = FALSE,
      group_structure = "direction",
      penalize_covariates = FALSE,
      family = "gaussian",
      lambda = "lambda.min",
      nfolds = 2,
      n_boot = 1,
      seed = 999,
      verbose = FALSE,
      parallel = FALSE,
      checkpoint_dir = checkpoint_dir,
      checkpoint_interval = 1,
      cleanup_checkpoint = FALSE
    ),
    fit_sgl_core = function(...) {
      list(pos_coef = c(0.5, 0), neg_coef = c(0, 0.25))
    },
    .package = "sglwqs"
  )
  
  expect_equal(resumed$n_successful, 1)
  expect_equal(unname(resumed$mean_pos_coef), c(0.5, 0))
  expect_equal(unname(resumed$mean_neg_coef), c(0, 0.25))
})
