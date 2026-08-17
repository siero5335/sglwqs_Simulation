test_that("bootstrap_sgl returns aggregated covariate coefficients", {
  set.seed(123)
  n <- 200
  p <- 6
  q <- 2

  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  cov_mat <- matrix(rnorm(n * q), n, q)
  colnames(cov_mat) <- paste0("z", seq_len(q))

  beta_x <- c(0.8, 0.5, 0, 0, -0.6, 0)
  beta_z <- c(0.4, -0.3)
  y <- as.numeric(X %*% beta_x + cov_mat %*% beta_z + rnorm(n, sd = 0.5))

  X_q <- quantile_transform(X, n_quantiles = 4)$q

  res <- suppressWarnings(sglwqs:::bootstrap_sgl(
    X_quantile = X_q,
    y = y,
    cov_matrix = cov_mat,
    var_names = colnames(X),
    cov_names = colnames(cov_mat),
    groups = NULL,
    group_by_compound = FALSE,
    group_structure = "direction",
    penalize_covariates = FALSE,
    family = "gaussian",
    lambda = "lambda.min",
    nfolds = 5,
    n_boot = 50,
    seed = 1,
    verbose = FALSE
  ))

  expect_false(is.null(res$mean_cov_coef))
  expect_named(res$mean_cov_coef, colnames(cov_mat))
  expect_equal(length(res$se_cov_coef), q)
  expect_true(all(!is.na(res$mean_cov_coef)))
  expect_lt(abs(res$mean_cov_coef["z1"] - 0.4), 0.2)
  expect_lt(abs(res$mean_cov_coef["z2"] - (-0.3)), 0.2)
})

test_that("bootstrap_sgl handles zero-covariate case", {
  set.seed(124)
  X <- matrix(rnorm(120), ncol = 3)
  colnames(X) <- c("x1", "x2", "x3")
  y <- rnorm(nrow(X))

  res <- sglwqs:::bootstrap_sgl(
    X_quantile = quantile_transform(X, n_quantiles = 4)$q,
    y = y,
    cov_matrix = NULL,
    var_names = colnames(X),
    cov_names = NULL,
    groups = NULL,
    group_by_compound = FALSE,
    group_structure = "direction",
    penalize_covariates = FALSE,
    family = "gaussian",
    lambda = "lambda.min",
    nfolds = 3,
    n_boot = 5,
    seed = 1,
    verbose = FALSE
  )

  expect_null(res$mean_cov_coef)
  expect_null(res$boot_cov_coef)
})

test_that("bootstrap_sgl aggregates index sums for grouped fits", {
  set.seed(1241)
  X <- matrix(rnorm(180 * 4), ncol = 4)
  colnames(X) <- c("x1", "x2", "x3", "x4")
  y <- 0.8 * X[, 1] - 0.6 * X[, 4] + rnorm(nrow(X), sd = 0.7)

  res <- sglwqs:::bootstrap_sgl(
    X_quantile = quantile_transform(X, n_quantiles = 4)$q,
    y = y,
    cov_matrix = NULL,
    var_names = colnames(X),
    cov_names = NULL,
    groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
    group_by_compound = TRUE,
    group_structure = "direction",
    penalize_covariates = FALSE,
    family = "gaussian",
    lambda = "lambda.min",
    nfolds = 3,
    n_boot = 8,
    seed = 11,
    verbose = FALSE
  )

  expect_true(is.numeric(res$mean_index_sum_pos))
  expect_true(is.numeric(res$mean_index_sum_neg))
  expect_named(res$mean_index_sum_by_group_pos, c("G1", "G2"))
  expect_named(res$mean_index_sum_by_group_neg, c("G1", "G2"))
  expect_named(res$se_index_sum_by_group_pos, c("G1", "G2"))
  expect_named(res$se_index_sum_by_group_neg, c("G1", "G2"))
  expect_equal(length(res$boot_pos_index_sum), 8)
  expect_equal(length(res$boot_neg_index_sum), 8)
})

test_that("sglwqs bootstrap keeps aggregated covariate summaries but omits matrices when requested", {
  set.seed(125)
  n <- 160
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    age = rnorm(n),
    bmi = rnorm(n)
  )
  y <- 0.7 * dat$x1 - 0.5 * dat$x2 + 0.5 * dat$age - 0.3 * dat$bmi + rnorm(n, sd = 0.8)

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = y,
    covariates = dat[, c("age", "bmi")],
    family = "gaussian",
    bootstrap = TRUE,
    n_boot = 10,
    nfolds = 3,
    seed = 1,
    verbose = FALSE,
    keep_boot_matrices = FALSE
  ))

  expect_true(is.numeric(fit$cov_coef))
  expect_true(is.numeric(fit$boot_info$mean_cov_coef))
  expect_true(all(c("age", "bmi") %in% names(fit$boot_info$mean_cov_coef)))
  expect_null(fit$boot_info$boot_cov_coef)
})

test_that("bootstrap checkpoint restores boot_cov_coef and MI pooling prefers bootstrap means", {
  checkpoint_dir <- tempfile("sglwqs-cov-checkpoint-")
  dir.create(checkpoint_dir)

  X_quantile <- matrix(rnorm(24), ncol = 2)
  colnames(X_quantile) <- c("x1", "x2")
  y <- rnorm(12)
  cov_matrix <- matrix(rnorm(24), ncol = 2)
  colnames(cov_matrix) <- c("z1", "z2")

  fail_once <- local({
    attempted <- FALSE
    function(...) {
      if (!attempted) {
        attempted <<- TRUE
        stop("forced bootstrap failure")
      }
      list(
        pos_coef = c(0.5, 0),
        neg_coef = c(0, 0.25),
        cov_coef = c(z1 = 0.4, z2 = -0.3)
      )
    }
  })

  suppressWarnings(
    expect_error(
      testthat::with_mocked_bindings(
        sglwqs:::bootstrap_sgl(
          X_quantile = X_quantile,
          y = y,
          cov_matrix = cov_matrix,
          var_names = colnames(X_quantile),
          cov_names = colnames(cov_matrix),
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
  expect_true("boot_cov_coef" %in% names(checkpoint_data))

  resumed <- testthat::with_mocked_bindings(
    sglwqs:::bootstrap_sgl(
      X_quantile = X_quantile,
      y = y,
      cov_matrix = cov_matrix,
      var_names = colnames(X_quantile),
      cov_names = colnames(cov_matrix),
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
      list(pos_coef = c(0.5, 0), neg_coef = c(0, 0.25), cov_coef = c(z1 = 0.4, z2 = -0.3))
    },
    .package = "sglwqs"
  )

  expect_equal(unname(resumed$mean_cov_coef), c(0.4, -0.3))

  fit1 <- structure(
    list(
      cov_coef = c(z1 = 9, z2 = 9),
      boot_info = list(mean_cov_coef = c(z1 = 0.4, z2 = -0.3))
    ),
    class = "sglwqs"
  )
  fit2 <- structure(
    list(
      cov_coef = c(z1 = 8, z2 = 8),
      boot_info = list(mean_cov_coef = c(z1 = 0.2, z2 = -0.1))
    ),
    class = "sglwqs"
  )
  mids_obj <- structure(
    list(
      validation = FALSE,
      bootstrap = TRUE,
      m = 2,
      n_successful = 2,
      family = "gaussian",
      outcome_var = "y",
      exposure_vars = c("x1", "x2"),
      covariate_vars = c("z1", "z2"),
      penalize_covariates = FALSE,
      fits = list(fit1, fit2),
      pooled = list(
        pos_index_sum = 0.5,
        neg_index_sum = 0.2,
        pos_weights = c(x1 = 0.7, x2 = 0.1),
        neg_weights = c(x1 = 0.0, x2 = 0.2),
        mi_selection_freq_pos = c(x1 = 1, x2 = 0.5),
        mi_selection_freq_neg = c(x1 = 0, x2 = 1)
      )
    ),
    class = "sglwqs_mids"
  )

  coef_df <- .build_mi_training_coef_table(mids_obj)
  z1_row <- coef_df[coef_df$Term == "z1", , drop = FALSE]
  z2_row <- coef_df[coef_df$Term == "z2", , drop = FALSE]

  expect_equal(z1_row$Estimate, mean(c(0.4, 0.2)))
  expect_equal(z2_row$Estimate, mean(c(-0.3, -0.1)))
})

test_that("pool_bootstrap_inference pools grouped bootstrap summaries", {
  fit1 <- structure(
    list(
      groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
      cov_names = c("z1"),
      boot_info = list(
        mean_index_sum_by_group_pos = list(G1 = 0.7, G2 = 0.2),
        mean_index_sum_by_group_neg = list(G1 = 0.1, G2 = 0.4),
        se_index_sum_by_group_pos = list(G1 = 0.1, G2 = 0.05),
        se_index_sum_by_group_neg = list(G1 = 0.08, G2 = 0.06),
        mean_cov_coef = c(z1 = 0.4),
        se_cov_coef = c(z1 = 0.12)
      ),
      validation_info = NULL,
      refit_info = NULL
    ),
    class = "sglwqs"
  )
  fit2 <- structure(
    list(
      groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
      cov_names = c("z1"),
      boot_info = list(
        mean_index_sum_by_group_pos = list(G1 = 0.5, G2 = 0.3),
        mean_index_sum_by_group_neg = list(G1 = 0.2, G2 = 0.2),
        se_index_sum_by_group_pos = list(G1 = 0.11, G2 = 0.04),
        se_index_sum_by_group_neg = list(G1 = 0.07, G2 = 0.05),
        mean_cov_coef = c(z1 = 0.2),
        se_cov_coef = c(z1 = 0.10)
      ),
      validation_info = NULL,
      refit_info = NULL
    ),
    class = "sglwqs"
  )

  pooled <- sglwqs:::pool_bootstrap_inference(list(fit1, fit2), verbose = FALSE)

  expect_identical(pooled$source_used, "boot_info")
  expect_true(isTRUE(pooled$has_group_results))
  expect_true(!is.null(pooled$group_results$G1$positive))
  expect_true(!is.null(pooled$covariates$z1))
  expect_true(is.finite(pooled$covariates$z1$estimate))
})
