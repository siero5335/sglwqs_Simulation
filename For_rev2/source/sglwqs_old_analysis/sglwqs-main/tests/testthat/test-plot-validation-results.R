test_that("non-MI validation plots can include covariates", {
  set.seed(20260414)
  n <- 180
  dat <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n),
    x3 = rnorm(n),
    x4 = rnorm(n),
    age = rnorm(n, mean = 50, sd = 10),
    sex = factor(sample(c("F", "M"), n, replace = TRUE))
  )
  dat$y <- 0.8 * dat$x1 - 0.6 * dat$x2 + 0.15 * dat$age +
    0.5 * (dat$sex == "M") + rnorm(n, sd = 0.8)

  fit <- suppressWarnings(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      covariates = dat[, c("age", "sex")],
      groups = list(G1 = c("x1", "x2"), G2 = c("x3", "x4")),
      validation = TRUE,
      train_prop = 0.6,
      seed = 42,
      verbose = FALSE
    )
  )

  p_val <- plot_validation_results(
    fit,
    include_covariates = TRUE
  )

  expect_s3_class(p_val, "ggplot")
  expect_true("Covariate" %in% unique(as.character(p_val$data$group)))
  expect_true("age" %in% unique(as.character(p_val$data$term)))
  expect_true(any(grepl("^sex", unique(as.character(p_val$data$term)))))

  p_combined <- plot_combined_results(
    fit,
    include_covariates = TRUE,
    layout = "horizontal"
  )

  expect_true(
    inherits(p_combined, "patchwork") ||
      inherits(p_combined, "ggplot") ||
      is.list(p_combined)
  )
})
