test_that("quantile_transform preserves empirical behavior without weights", {
  dat <- make_simple_data(n = 120, seed = 10)
  default_res <- quantile_transform(dat[, c("x1", "x2")], n_quantiles = 4)
  empirical_res <- quantile_transform(
    dat[, c("x1", "x2")],
    n_quantiles = 4,
    method = "empirical"
  )
  
  expect_equal(default_res$q, empirical_res$q)
})


test_that("equal weights are broadly consistent and extreme weights shift breaks", {
  dat <- data.frame(x = rnorm(200))
  empirical_res <- quantile_transform(dat, n_quantiles = 4, method = "empirical")
  equal_weight_res <- quantile_transform(
    dat,
    n_quantiles = 4,
    weights = rep(1, nrow(dat))
  )
  
  same_bin_prop <- mean(empirical_res$q[, 1] == equal_weight_res$q[, 1], na.rm = TRUE)
  expect_gte(same_bin_prop, 0.9)
  
  shifted <- data.frame(x = 1:100)
  weighted_res <- expect_warning(
    quantile_transform(
      shifted,
      n_quantiles = 4,
      weights = c(rep(1, 99), 100)
    ),
    "conservative fallback"
  )
  unweighted_res <- quantile_transform(
    shifted,
    n_quantiles = 4,
    method = "empirical"
  )
  
  expect_false(isTRUE(all.equal(weighted_res$breaks[[1]]$breaks, unweighted_res$breaks[[1]]$breaks)))
})


test_that("quantile_transform rejects non-numeric input columns", {
  bad_dat <- data.frame(
    x = rnorm(20),
    grp = rep(c("a", "b"), each = 10),
    stringsAsFactors = FALSE
  )
  
  expect_error(
    quantile_transform(bad_dat, n_quantiles = 4),
    "Non-numeric columns"
  )
})


test_that("obs_weights are accepted by selection and full glm refit works", {
  dat <- make_simple_data(n = 180, seed = 20)
  obs_weights <- seq(1, 2, length.out = nrow(dat))
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = data.frame(age = rnorm(nrow(dat))),
    obs_weights = obs_weights,
    refit = "full",
    refit_engine = "glm",
    nfolds = 3,
    nlambda = 20,
    seed = 99,
    verbose = FALSE
  )
  
  expect_s3_class(fit, "sglwqs")
  expect_equal(fit$obs_weights, obs_weights)
  expect_equal(fit$refit, "full")
  expect_equal(fit$refit_info$engine, "glm")
  expect_false(is.null(fit$refit_info$refit_fit))
  expect_no_error(capture.output(summary(fit)))
})


test_that("survey mode rejects unsupported combinations", {
  dat <- make_simple_data(n = 120, seed = 30)
  dummy_design <- list(variables = data.frame(id = seq_len(nrow(dat))))
  
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      refit_engine = "svyglm",
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "survey_design"
  )
  
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "validation",
      refit_engine = "svyglm",
      survey_design = dummy_design,
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "Survey refit"
  )
  
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      refit_engine = "svyglm",
      survey_design = dummy_design,
      bootstrap = TRUE,
      n_boot = 5,
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "Survey refit"
  )
  
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      refit_engine = "svyglm",
      survey_design = dummy_design,
      parallel = TRUE,
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "Survey refit"
  )
})


test_that("full svyglm refit works when survey is available", {
  testthat::skip_if_not_installed("survey")
  
  dat <- make_simple_data(n = 160, seed = 40)
  cov_df <- data.frame(age = rnorm(nrow(dat)))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 101,
    verbose = FALSE
  ))
  
  expect_equal(fit$refit_info$engine, "svyglm")
  expect_s3_class(fit$refit_info$refit_fit, "svyglm")
  expect_no_error(generics::glance(fit))
  expect_no_error(generics::augment(
    fit,
    data = dat[1:10, c("x1", "x2", "x3", "x4")],
    covariates = cov_df[1:10, , drop = FALSE]
  ))
  expect_no_error(capture.output(summary(fit)))
})


test_that("full glm refit predictions and augment use the refit model", {
  dat <- make_simple_data(n = 180, seed = 55)
  cov_df <- data.frame(age = rnorm(nrow(dat)))
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "glm",
    nfolds = 3,
    nlambda = 20,
    seed = 123,
    verbose = FALSE
  ))
  
  new_exp <- dat[1:12, c("x1", "x2", "x3", "x4")]
  new_cov <- cov_df[1:12, , drop = FALSE]
  qt <- quantile_transform(
    new_exp,
    n_quantiles = fit$n_quantiles,
    var_names = fit$var_names,
    breaks_list = fit$quantile_breaks
  )
  pred_data <- sglwqs:::.build_validation_prediction_data(
    fit,
    X_quantile = qt$q,
    cov_matrix = as.matrix(new_cov)
  )
  
  expected_link <- as.numeric(stats::predict(fit$refit_info$refit_fit, newdata = pred_data, type = "link"))
  expected_resp <- as.numeric(stats::predict(fit$refit_info$refit_fit, newdata = pred_data, type = "response"))
  
  expect_equal(
    predict(fit, newdata = new_exp, covariates = new_cov, type = "link"),
    expected_link
  )
  expect_equal(
    predict(fit, newdata = new_exp, covariates = new_cov, type = "response"),
    expected_resp
  )
  
  ag <- generics::augment(fit, data = new_exp, covariates = new_cov)
  expect_equal(ag$.fitted, expected_resp)
})


test_that("survey binomial refit uses quasibinomial and supports prediction uncertainty", {
  testthat::skip_if_not_installed("survey")
  
  dat <- make_binomial_data(n = 220, seed = 66)
  cov_df <- data.frame(age = rnorm(nrow(dat)))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    covariates = cov_df,
    family = "binomial",
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 321,
    verbose = FALSE
  )
  
  expect_equal(fit$refit_info$refit_fit$family$family, "quasibinomial")
  
  new_exp <- dat[1:15, c("x1", "x2", "x3")]
  new_cov <- cov_df[1:15, , drop = FALSE]
  qt <- quantile_transform(
    new_exp,
    n_quantiles = fit$n_quantiles,
    var_names = fit$var_names,
    breaks_list = fit$quantile_breaks
  )
  pred_data <- sglwqs:::.build_validation_prediction_data(
    fit,
    X_quantile = qt$q,
    cov_matrix = as.matrix(new_cov)
  )
  
  expected_link <- as.numeric(stats::predict(fit$refit_info$refit_fit, newdata = pred_data, type = "link"))
  expected_resp <- as.numeric(stats::predict(fit$refit_info$refit_fit, newdata = pred_data, type = "response"))
  
  got_link <- predict(fit, newdata = new_exp, covariates = new_cov, type = "link")
  got_resp <- predict(fit, newdata = new_exp, covariates = new_cov, type = "response")
  expect_equal(got_link, expected_link)
  expect_equal(got_resp, expected_resp)
  
  se_link <- predict(fit, newdata = new_exp, covariates = new_cov, type = "link", se.fit = TRUE)
  se_resp <- predict(fit, newdata = new_exp, covariates = new_cov, type = "response", se.fit = TRUE)
  expect_equal(se_link$source, "svyglm_refit")
  expect_equal(se_resp$source, "svyglm_refit")
  expect_length(se_link$fit, 15)
  expect_length(se_resp$fit, 15)
  expect_true(all(is.finite(se_link$se.fit)))
  expect_true(all(is.finite(se_resp$se.fit)))
  
  pred_int <- predict_interval(
    fit,
    newdata = new_exp,
    covariates = new_cov,
    type = "response",
    level = 0.95
  )
  expect_true(all(c("fit", "lwr", "upr", "source") %in% names(pred_int)))
  
  ag <- generics::augment(fit, data = new_exp, covariates = new_cov)
  expect_equal(ag$.fitted, expected_resp)
})


test_that("survey mode blocks unweighted metrics and calibration and checks row alignment", {
  testthat::skip_if_not_installed("survey")
  
  dat <- make_binomial_data(n = 160, seed = 77)
  rownames(dat) <- paste0("row", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    covariates = cov_df,
    family = "binomial",
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 88,
    verbose = FALSE
  )
  
  expect_error(validation_metrics(fit), "Survey-aware predictive metrics")
  expect_error(calibrate(fit), "Survey-aware calibration")
  
  bad_design_df <- design_df
  rownames(bad_design_df) <- rev(rownames(bad_design_df))
  bad_des <- survey::svydesign(ids = ~1, weights = ~w, data = bad_design_df)
  
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3")],
      y = dat$y,
      covariates = cov_df,
      family = "binomial",
      refit = "full",
      refit_engine = "svyglm",
      survey_design = bad_des,
      nfolds = 3,
      nlambda = 20,
      seed = 89,
      verbose = FALSE
    ),
    "Row names of `survey_design\\$variables` do not match"
  )
})


test_that("svyglm refit with stratified clustered design produces design-based SEs", {
  testthat::skip_if_not_installed("survey")

  set.seed(500)
  n <- 200
  dat <- make_simple_data(n = n, seed = 500)
  cov_df <- data.frame(age = rnorm(n))

  # Stratified clustered design (NHANES-like)
  strata <- rep(1:10, each = n / 10)
  cluster <- rep(1:(n / 5), each = 5)
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(n, 0.5, 3),
    strata = strata,
    cluster = cluster
  )
  des <- survey::svydesign(
    ids = ~cluster, strata = ~strata, weights = ~w,
    data = design_df, nest = TRUE
  )

  fit_svy <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 501,
    verbose = FALSE
  ))

  expect_s3_class(fit_svy$refit_info$refit_fit, "svyglm")
  expect_no_error(capture.output(summary(fit_svy)))

  # SE from svyglm should differ from ordinary glm due to design effect
  fit_glm <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "glm",
    nfolds = 3,
    nlambda = 20,
    seed = 501,
    verbose = FALSE
  ))

  svy_se <- summary(fit_svy$refit_info$refit_fit)$coefficients[, "Std. Error"]
  glm_se <- summary(fit_glm$refit_info$refit_fit)$coefficients[, "Std. Error"]

  # Design-based SEs should generally differ from model-based SEs
  expect_false(isTRUE(all.equal(svy_se, glm_se, tolerance = 0.01)))
})


test_that("svyglm prediction intervals use survey residual degrees of freedom", {
  testthat::skip_if_not_installed("survey")
  
  set.seed(700)
  n <- 180
  dat <- make_simple_data(n = n, seed = 700)
  cov_df <- data.frame(age = rnorm(n))
  strata <- rep(1:6, each = n / 6)
  cluster <- rep(1:(n / 5), each = 5)
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(n, 0.5, 2.5),
    strata = strata,
    cluster = cluster
  )
  des <- survey::svydesign(
    ids = ~cluster, strata = ~strata, weights = ~w,
    data = design_df, nest = TRUE
  )
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 701,
    verbose = FALSE
  ))
  
  new_exp <- dat[1:10, c("x1", "x2", "x3", "x4")]
  new_cov <- cov_df[1:10, , drop = FALSE]
  pred_se <- predict(
    fit,
    newdata = new_exp,
    covariates = new_cov,
    type = "link",
    se.fit = TRUE
  )
  pred_int <- predict_interval(
    fit,
    newdata = new_exp,
    covariates = new_cov,
    type = "link",
    level = 0.95
  )
  
  df_resid <- stats::df.residual(fit$refit_info$refit_fit)
  crit <- stats::qt(0.975, df = df_resid)
  
  expect_equal(pred_int$lwr, pred_se$fit - crit * pred_se$se.fit, tolerance = 1e-6)
  expect_equal(pred_int$upr, pred_se$fit + crit * pred_se$se.fit, tolerance = 1e-6)
})


test_that("rubin_pool honors supplied complete-data degrees of freedom", {
  pooled_small_df <- sglwqs:::rubin_pool(
    estimates = c(0.8, 1.1, 1.0),
    variances = c(0.04, 0.03, 0.05),
    m = 3,
    conf_level = 0.95,
    df_complete = 6
  )
  pooled_large_df <- sglwqs:::rubin_pool(
    estimates = c(0.8, 1.1, 1.0),
    variances = c(0.04, 0.03, 0.05),
    m = 3,
    conf_level = 0.95,
    df_complete = 1000
  )
  
  expect_true(is.finite(pooled_small_df$df))
  expect_lte(pooled_small_df$df, 6)
  expect_lt(pooled_small_df$df, pooled_large_df$df)
})


test_that("input validation catches bad y values", {
  dat <- make_simple_data(n = 100, seed = 600)

  # y length mismatch
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2")],
      y = dat$y[1:50],
      nfolds = 3, nlambda = 10, verbose = FALSE
    ),
    "same length"
  )

  # y contains Inf
  bad_y <- dat$y
  bad_y[1] <- Inf
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2")],
      y = bad_y,
      nfolds = 3, nlambda = 10, verbose = FALSE
    ),
    "Inf"
  )

  # y all NA
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2")],
      y = rep(NA_real_, 100),
      nfolds = 3, nlambda = 10, verbose = FALSE
    ),
    "entirely NA"
  )

  # binomial with non-binary y
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2")],
      y = dat$y,
      family = "binomial",
      nfolds = 3, nlambda = 10, verbose = FALSE
    ),
    "binary"
  )
})
