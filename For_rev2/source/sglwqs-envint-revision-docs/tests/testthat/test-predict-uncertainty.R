test_that("predict returns se.fit with bootstrap uncertainty", {
  dat <- make_uncertainty_data(seed = 1001)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    bootstrap = TRUE,
    n_boot = 8,
    keep_boot_matrices = TRUE,
    nfolds = 3,
    nlambda = 20,
    seed = 123,
    verbose = FALSE
  ))
  
  pred <- predict(
    fit,
    newdata = dat[1:12, c("x1", "x2", "x3")],
    type = "response",
    se.fit = TRUE
  )
  
  expect_type(pred, "list")
  expect_true(all(c("fit", "se.fit", "source") %in% names(pred)))
  expect_length(pred$fit, 12)
  expect_length(pred$se.fit, 12)
  expect_true(all(is.finite(pred$se.fit)))
})


test_that("predict_interval returns fit and bounds", {
  dat <- make_uncertainty_data(seed = 1002)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    bootstrap = TRUE,
    n_boot = 8,
    keep_boot_matrices = TRUE,
    nfolds = 3,
    nlambda = 20,
    seed = 124,
    verbose = FALSE
  ))
  
  pi <- predict_interval(
    fit,
    newdata = dat[1:10, c("x1", "x2", "x3")],
    type = "response",
    level = 0.90
  )
  
  expect_true(is.data.frame(pi))
  expect_true(all(c("fit", "lwr", "upr", "source", "interval") %in% names(pi)))
  expect_equal(nrow(pi), 10)
  expect_true(all(pi$lwr <= pi$fit))
  expect_true(all(pi$fit <= pi$upr))
})


test_that("prediction interval for gaussian validation model is available", {
  dat <- make_uncertainty_data(seed = 1003)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    validation = TRUE,
    train_prop = 0.6,
    nfolds = 3,
    nlambda = 20,
    seed = 125,
    verbose = FALSE
  ))
  
  pi <- predict(
    fit,
    newdata = dat[1:9, c("x1", "x2", "x3")],
    type = "response",
    interval = "prediction",
    level = 0.95
  )
  
  expect_true(is.data.frame(pi))
  expect_equal(nrow(pi), 9)
  expect_true(all(c("fit", "lwr", "upr") %in% names(pi)))
})
