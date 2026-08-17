test_that("tidy/glance/augment return expected structures", {
  dat <- make_simple_data(seed = 987)
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    nfolds = 3,
    nlambda = 20,
    seed = 11,
    verbose = FALSE
  )
  
  td <- generics::tidy(fit, what = "weights")
  gl <- generics::glance(fit)
  ag <- generics::augment(fit, data = dat[1:20, c("x1", "x2", "x3", "x4")])
  
  expect_true(is.data.frame(td))
  expect_true(all(c("term", "weight", "direction") %in% names(td)))
  expect_true(is.data.frame(gl))
  expect_true(all(c("nobs", "n_exposures", "family") %in% names(gl)))
  expect_true(is.data.frame(ag))
  expect_true(all(c(".fitted", ".wqs_pos", ".wqs_neg") %in% names(ag)))
  expect_equal(nrow(ag), 20)
})


test_that("validation_metrics and calibrate work for binomial model", {
  dat <- make_binomial_data(seed = 1234)
  fit <- sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    family = "binomial",
    validation = TRUE,
    train_prop = 0.6,
    nfolds = 3,
    nlambda = 20,
    seed = 21,
    verbose = FALSE
  )
  
  m <- validation_metrics(fit)
  cal <- calibrate(fit, n_bins = 8)
  
  expect_true(is.data.frame(m))
  expect_true(all(c("auc", "brier", "accuracy") %in% names(m)))
  expect_true(m$auc >= 0 && m$auc <= 1)
  
  expect_true(is.data.frame(cal))
  expect_true(all(c("bin", "pred_mean", "obs_rate") %in% names(cal)))
  expect_gte(nrow(cal), 2)
})


test_that("autoplot and as.data.frame helpers run", {
  dat <- make_simple_data(seed = 222)
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    bootstrap = TRUE,
    n_boot = 5,
    nfolds = 3,
    nlambda = 20,
    seed = 31,
    verbose = FALSE
  ))
  
  p1 <- ggplot2::autoplot(fit, type = "cv")
  p2 <- ggplot2::autoplot(fit, type = "weights")
  df <- as.data.frame(fit, what = "weights")
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_true(is.data.frame(df))
  expect_true(all(c("term", "weight") %in% names(df)))
})
