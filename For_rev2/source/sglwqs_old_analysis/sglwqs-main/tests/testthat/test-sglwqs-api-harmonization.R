test_that("sglwqs accepts data plus named selectors", {
  data("sglwqs_example")

  exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))

  fit <- sglwqs(
    data = sglwqs_example[1:120, ],
    X = exp_vars,
    y = "outcome_cont",
    covariates = "sex",
    nfolds = 3,
    nlambda = 20,
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(fit, "sglwqs")
  expect_identical(fit$var_names, exp_vars)
  expect_identical(fit$exposure_vars, exp_vars)
  expect_identical(fit$outcome_var, "outcome_cont")
  expect_identical(fit$covariate_vars, "sex")
})

test_that("sglwqs formula interface resolves variables from data", {
  data("sglwqs_example")

  fit <- suppressWarnings(sglwqs(
    data = sglwqs_example[1:120, ],
    formula = outcome_cont ~ metal1 + metal2 + pesticide1 + pesticide2 | sex,
    nfolds = 3,
    nlambda = 20,
    seed = 321,
    verbose = FALSE
  ))

  expect_s3_class(fit, "sglwqs")
  expect_identical(fit$var_names, c("metal1", "metal2", "pesticide1", "pesticide2"))
  expect_identical(fit$outcome_var, "outcome_cont")
  expect_identical(fit$covariate_vars, "sex")
})
