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
  
  fit <- suppressWarnings(sglwqs(
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
  ))
  
  expect_s3_class(fit, "sglwqs")
  expect_equal(fit$obs_weights, obs_weights)
  expect_equal(fit$refit, "full")
  expect_equal(fit$refit_info$engine, "glm")
  expect_false(is.null(fit$refit_info$refit_fit))
  expect_no_error(capture.output(summary(fit)))
})


test_that("survey mode validates refit-engine inputs", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 30)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  design_df <- data.frame(
    y = dat$y,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
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
      refit = "full",
      refit_engine = "glm",
      survey_design = des,
      analysis_id = rownames(dat),
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "requires `refit_engine = \"svyglm\"`"
  )

  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      validation = TRUE,
      nfolds = 3,
      nlambda = 20,
      verbose = FALSE
    ),
    "validation = TRUE"
  )
})


test_that("full svyglm refit works when survey is available", {
  testthat::skip_if_not_installed("survey")
  
  dat <- make_simple_data(n = 160, seed = 40)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    analysis_id = rownames(dat),
    nfolds = 3,
    nlambda = 20,
    seed = 101,
    verbose = FALSE
  ))
  
  expect_equal(fit$refit_info$engine, "svyglm")
  expect_s3_class(fit$refit_info$refit_fit, "svyglm")
  expect_equal(fit$analysis_id, rownames(dat))
  expect_equal(as.numeric(fit$analysis_weights), design_df$w)
  expect_equal(as.numeric(fit$quantile_weights), design_df$w)
  expect_equal(fit$quantile_weights_source, "survey_design")
  expect_equal(as.numeric(fit$obs_weights), design_df$w / mean(design_df$w))
  expect_true(is.list(fit$survey_info))

  refit_data <- data.frame(
    fit$refit_info$wqs_indices,
    y = dat$y,
    age = cov_df$age,
    check.names = FALSE
  )
  manual <- survey::svyglm(
    stats::as.formula(fit$refit_info$formula),
    design = sglwqs:::.inject_design_variables(des, refit_data),
    family = stats::gaussian()
  )
  expect_equal(stats::coef(fit$refit_info$refit_fit), stats::coef(manual))
  expect_equal(summary(fit$refit_info$refit_fit)$coefficients, summary(manual)$coefficients)
  expect_equal(
    as.numeric(stats::predict(fit$refit_info$refit_fit, newdata = refit_data, type = "link")),
    as.numeric(stats::predict(manual, newdata = refit_data, type = "link"))
  )
  expect_no_error(generics::glance(fit))
  expect_no_error(generics::augment(
    fit,
    data = dat[1:10, c("x1", "x2", "x3", "x4")],
    covariates = cov_df[1:10, , drop = FALSE]
  ))
  expect_no_error(capture.output(summary(fit)))
})


test_that("weighted survey selection records lambda diagnostics", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 141)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    analysis_id = rownames(dat),
    nfolds = 3,
    nlambda = 10,
    seed = 141,
    verbose = FALSE
  ))

  diag <- fit$selection_diagnostics
  expect_true(is.list(diag))
  expect_true(isTRUE(diag$survey_mode))
  expect_true(isTRUE(diag$weighted_selection))
  expect_true(is.finite(diag$lambda_path_min))
  expect_true(is.finite(diag$lambda_path_max))
  expect_true(diag$lambda_path_length >= 1)
  expect_true(all(c(
    "selection_lambda_path",
    "selection_exposure_nonzero"
  ) %in% names(fit$diagnostics)))
})


test_that("explicit lambda_path is passed to sparsegl backend separately from lambda", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 143)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  lambda_path <- exp(seq(log(0.1), log(1e-3), length.out = 8))

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    covariates = cov_df,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    analysis_id = rownames(dat),
    lambda = "lambda.min",
    lambda_path = lambda_path,
    nfolds = 3,
    seed = 143,
    verbose = FALSE
  ))

  expect_true(length(fit$fit$lambda) <= length(lambda_path))
  expect_true(all(signif(fit$fit$lambda, 12) %in% signif(lambda_path, 12)))
  expect_equal(fit$selection_diagnostics$lambda_path_source, "user_explicit")
  expect_equal(
    fit$selection_diagnostics$lambda_path_length,
    length(fit$fit$lambda)
  )
})


test_that("lambda_path validates early and bootstrap errors are classifiable", {
  expect_error(
    sglwqs(
      X = matrix(rnorm(40), ncol = 2),
      y = rnorm(20),
      lambda_path = c(0.1, 0, 0.001),
      verbose = FALSE
    ),
    "lambda_path"
  )
  expect_equal(
    sglwqs:::.classify_boot_error("object 'boundary' not found"),
    "backend_convergence"
  )
  expect_equal(
    sglwqs:::.classify_boot_error("sparsegl_irls failed"),
    "backend_convergence"
  )
  expect_equal(
    sglwqs:::.classify_boot_error("lambda sequence failed"),
    "lambda_path"
  )
})


test_that("weighted survey boundary lambda selection warns", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 144)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)

  fit <- NULL
  expect_warning(
    fit <- sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      covariates = cov_df,
      refit = "full",
      refit_engine = "svyglm",
      survey_design = des,
      analysis_id = rownames(dat),
      lambda = 1e-5,
      lambda_path = c(1e-4, 1e-5),
      nfolds = 3,
      seed = 144,
      verbose = FALSE
    ),
    "path boundary"
  )
  expect_true(isTRUE(fit$selection_diagnostics$selected_lambda_at_path_boundary))
})


test_that("weighted survey all-zero exposure selection is explicit", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 142)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)

  fit <- NULL
  expect_warning(
    fit <- sglwqs(
        X = dat[, c("x1", "x2", "x3", "x4")],
        y = dat$y,
        covariates = cov_df,
        refit = "full",
        refit_engine = "svyglm",
        survey_design = des,
        analysis_id = rownames(dat),
        lambda = 1e6,
        nfolds = 3,
        nlambda = 10,
        seed = 142,
        verbose = FALSE
      ),
    "all-zero exposure coefficients"
  )

  expect_true(isTRUE(fit$selection_diagnostics$all_zero_exposure))
  expect_equal(fit$selection_diagnostics$nonzero_positive_coef, 0)
  expect_equal(fit$selection_diagnostics$nonzero_negative_coef, 0)
  expect_equal(fit$selection_diagnostics$positive_weight_sum, 0)
  expect_equal(fit$selection_diagnostics$negative_weight_sum, 0)
  expect_true("selection_all_zero_exposure" %in% names(fit$diagnostics))
})


test_that("survey bootstrap preparation selects and scales replicate weights", {
  testthat::skip_if_not_installed("survey")

  dat <- data.frame(
    y = rnorm(24),
    w = runif(24, 0.5, 2)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = dat)

  naive <- sglwqs:::.prepare_boot_weights(
    n = nrow(dat),
    n_boot = 4,
    boot_method = "auto",
    survey_design = NULL,
    family = "gaussian",
    y = dat$y,
    seed = 11
  )
  expect_equal(naive$method_used, "naive")
  expect_equal(dim(naive$index_matrix), c(nrow(dat), 4))

  svrep <- sglwqs:::.prepare_boot_weights(
    n = nrow(dat),
    n_boot = 3,
    boot_method = "auto",
    survey_design = des,
    svrep_type = "auto",
    family = "gaussian",
    y = dat$y,
    seed = 12
  )
  expect_equal(svrep$method_used, "svrep")
  expect_equal(svrep$svrep_type_used, "bootstrap")
  expect_equal(dim(svrep$weight_matrix), c(nrow(dat), 3))
  expect_equal(as.numeric(colMeans(svrep$weight_matrix)), rep(1, 3), tolerance = 1e-8)
  expect_true(is.finite(svrep$variance_scale))
  expect_length(svrep$variance_rscales, 3)
})


test_that("survey replicate bootstrap uses full-sample center and survey variance", {
  testthat::skip_if_not_installed("survey")

  n <- 20
  X_quantile <- matrix(runif(n * 2), ncol = 2)
  colnames(X_quantile) <- c("x1", "x2")
  dat <- data.frame(y = rnorm(n), w = runif(n, 0.5, 2))
  des <- survey::svydesign(ids = ~1, weights = ~w, data = dat)
  fit_calls <- 0L

  res <- testthat::with_mocked_bindings(
    sglwqs:::bootstrap_sgl(
      X_quantile = X_quantile,
      y = dat$y,
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
      n_boot = 2,
      seed = 13,
      verbose = FALSE,
      boot_method = "svrep",
      survey_design = des,
      svrep_type = "bootstrap",
      obs_weights = rep(1, n)
    ),
    fit_sgl_core = function(...) {
      fit_calls <<- fit_calls + 1L
      pos <- if (fit_calls <= 2L) {
        c(x1 = 10 + fit_calls, x2 = 0)
      } else {
        c(x1 = 1, x2 = 0)
      }
      list(
        pos_coef = pos,
        neg_coef = c(x1 = 0, x2 = 0),
        cov_coef = NULL
      )
    },
    .package = "sglwqs"
  )

  expect_equal(res$method, "svrep")
  expect_equal(res$n_successful, 2)
  expect_equal(res$n_boot_actual, 2)
  expect_true(isTRUE(res$weights_normalized))
  expect_equal(as.numeric(res$mean_pos_coef["x1"]), 1)
  expect_equal(as.numeric(res$svrep_center_pos_coef["x1"]), 1)
  expect_true(is.finite(res$se_pos_coef["x1"]))
  expect_equal(res$bootstrap_design, "survey_replicate_weights")
})


test_that("summary_bootstrap uses stored svrep center and SE when matrices are kept", {
  boot_info <- list(
    method = "svrep",
    mean_pos_coef = c(x1 = 1, x2 = 0),
    mean_neg_coef = c(x1 = 0, x2 = 0),
    se_pos_coef = c(x1 = 0.2, x2 = 0.1),
    se_neg_coef = c(x1 = 0.05, x2 = 0.05),
    selection_freq_pos = c(x1 = 1, x2 = 0),
    selection_freq_neg = c(x1 = 0, x2 = 0),
    boot_pos_coef = rbind(c(x1 = 10, x2 = 0), c(x1 = 12, x2 = 0)),
    boot_neg_coef = rbind(c(x1 = 0, x2 = 0), c(x1 = 0, x2 = 0)),
    boot_success = c(TRUE, TRUE),
    n_successful = 2,
    variance_df = 8
  )
  fit <- structure(
    list(
      bootstrap = TRUE,
      boot_info = boot_info,
      var_names = c("x1", "x2"),
      pos_weights = c(x1 = 1, x2 = 0),
      neg_weights = c(x1 = 0, x2 = 0),
      pos_coef = c(x1 = 99, x2 = 0),
      neg_coef = c(x1 = 0, x2 = 0),
      groups = NULL
    ),
    class = "sglwqs"
  )

  out <- summary_bootstrap(fit, conf_level = 0.95)
  pos_x1 <- out[out$variable == "x1" & out$direction == "positive", ]
  crit <- stats::qt(0.975, df = 8)

  expect_equal(pos_x1$mean_coef, 1)
  expect_equal(pos_x1$se_coef, 0.2)
  expect_equal(pos_x1$ci_lower, 1 - crit * 0.2)
  expect_equal(pos_x1$ci_upper, 1 + crit * 0.2)
})


test_that("validation_glm can use svyglm on a validation survey subset", {
  testthat::skip_if_not_installed("survey")

  n <- 30
  ids <- paste0("v", seq_len(n))
  X_quantile <- matrix(runif(n * 2), ncol = 2, dimnames = list(ids, c("x1", "x2")))
  y <- rnorm(n)
  des <- survey::svydesign(
    ids = ~1,
    weights = ~w,
    data = data.frame(
      y = y,
      w = runif(n, 0.5, 2),
      .analysis_id = ids,
      row.names = ids
    )
  )

  val <- sglwqs:::validation_glm(
    X_quantile_val = X_quantile,
    y_val = y,
    cov_matrix_val = NULL,
    pos_weights = c(x1 = 0.7, x2 = 0.3),
    neg_weights = c(x1 = 0, x2 = 0),
    family = "gaussian",
    groups = NULL,
    group_inference = FALSE,
    engine = "svyglm",
    survey_design = des,
    analysis_id = ids
  )

  expect_equal(val$engine, "svyglm")
  expect_s3_class(val$refit_fit, "svyglm")
  expect_equal(val$n_val, n)
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
  rownames(dat) <- paste0("bin", seq_len(nrow(dat)))
  cov_df <- data.frame(age = rnorm(nrow(dat)), row.names = rownames(dat))
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(nrow(dat), 0.5, 2),
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    covariates = cov_df,
    family = "binomial",
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    analysis_id = rownames(dat),
    nfolds = 3,
    nlambda = 20,
    seed = 321,
    verbose = FALSE
  ))
  
  expect_equal(fit$refit_info$refit_fit$family$family, "quasibinomial")
  expect_equal(as.numeric(fit$analysis_weights), design_df$w)

  refit_data <- data.frame(
    fit$refit_info$wqs_indices,
    y = dat$y,
    age = cov_df$age,
    check.names = FALSE
  )
  manual <- survey::svyglm(
    stats::as.formula(fit$refit_info$formula),
    design = sglwqs:::.inject_design_variables(des, refit_data),
    family = stats::quasibinomial()
  )
  expect_equal(stats::coef(fit$refit_info$refit_fit), stats::coef(manual))
  expect_equal(summary(fit$refit_info$refit_fit)$coefficients, summary(manual)$coefficients)
  
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
  expected_se_link <- sglwqs:::.extract_prediction_with_se(stats::predict(
    fit$refit_info$refit_fit,
    newdata = pred_data,
    type = "link",
    se.fit = TRUE
  ))
  expected_se_resp <- sglwqs:::.extract_prediction_with_se(stats::predict(
    fit$refit_info$refit_fit,
    newdata = pred_data,
    type = "response",
    se.fit = TRUE
  ))
  expect_equal(se_link$source, "svyglm_refit")
  expect_equal(se_resp$source, "svyglm_refit")
  expect_length(se_link$fit, 15)
  expect_length(se_resp$fit, 15)
  expect_equal(se_link$se.fit, expected_se_link$se.fit)
  expect_equal(se_resp$se.fit, expected_se_resp$se.fit)
  
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
  
  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3")],
    y = dat$y,
    covariates = cov_df,
    family = "binomial",
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    analysis_id = rownames(dat),
    nfolds = 3,
    nlambda = 20,
    seed = 88,
    verbose = FALSE
  ))
  
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
      analysis_id = rownames(dat),
      nfolds = 3,
      nlambda = 20,
      seed = 89,
      verbose = FALSE
    ),
    "`analysis_id` does not match"
  )
})


test_that("survey alignment guards analysis IDs and row-name fallback", {
  testthat::skip_if_not_installed("survey")

  dat <- make_simple_data(n = 120, seed = 91)
  rownames(dat) <- paste0("id", seq_len(nrow(dat)))
  design_df <- data.frame(
    y = dat$y,
    w = runif(nrow(dat), 0.5, 2),
    row.names = rownames(dat)
  )
  des <- survey::svydesign(ids = ~1, weights = ~w, data = design_df)

  fit <- suppressWarnings(sglwqs(
    X = dat[, c("x1", "x2", "x3", "x4")],
    y = dat$y,
    refit = "full",
    refit_engine = "svyglm",
    survey_design = des,
    nfolds = 3,
    nlambda = 20,
    seed = 91,
    verbose = FALSE
  ))
  expect_equal(fit$analysis_id, rownames(dat))
  expect_s3_class(fit$refit_info$refit_fit, "svyglm")

  shifted_id <- c(rownames(dat)[-1], rownames(dat)[1])
  shifted_design_df <- design_df
  shifted_design_df$.analysis_id <- shifted_id
  shifted_des <- survey::svydesign(ids = ~1, weights = ~w, data = shifted_design_df)
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      refit_engine = "svyglm",
      survey_design = shifted_des,
      analysis_id = shifted_id,
      nfolds = 3,
      nlambda = 20,
      seed = 92,
      verbose = FALSE
    ),
    "`analysis_id` does not match the row names"
  )

  empty_id <- rownames(dat)
  empty_id[1] <- ""
  expect_error(
    sglwqs(
      X = dat[, c("x1", "x2", "x3", "x4")],
      y = dat$y,
      refit = "full",
      refit_engine = "svyglm",
      survey_design = des,
      analysis_id = empty_id,
      nfolds = 3,
      nlambda = 20,
      seed = 93,
      verbose = FALSE
    ),
    "missing or empty"
  )

  data_api <- data.frame(
    id = paste0("api", seq_len(nrow(dat))),
    dat,
    check.names = FALSE
  )
  rownames(data_api) <- NULL
  data_api_design <- survey::svydesign(
    ids = ~1,
    weights = ~w,
    data = data.frame(
      .analysis_id = data_api$id,
      w = design_df$w
    )
  )
  fit_data_api <- suppressWarnings(sglwqs(
    data = data_api,
    exposure_vars = c("x1", "x2", "x3", "x4"),
    outcome_var = "y",
    refit = "full",
    refit_engine = "svyglm",
    survey_design = data_api_design,
    analysis_id = 1,
    nfolds = 3,
    nlambda = 20,
    seed = 94,
    verbose = FALSE
  ))
  expect_equal(fit_data_api$analysis_id, data_api$id)
  expect_s3_class(fit_data_api$refit_info$refit_fit, "svyglm")
})


test_that("svyglm refit with stratified clustered design produces design-based SEs", {
  testthat::skip_if_not_installed("survey")

  set.seed(500)
  n <- 200
  dat <- make_simple_data(n = n, seed = 500)
  rownames(dat) <- paste0("s", seq_len(n))
  cov_df <- data.frame(age = rnorm(n), row.names = rownames(dat))

  # Stratified clustered design (NHANES-like)
  strata <- rep(1:10, each = n / 10)
  cluster <- rep(1:(n / 5), each = 5)
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(n, 0.5, 3),
    strata = strata,
    cluster = cluster,
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
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
    analysis_id = rownames(dat),
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
  rownames(dat) <- paste0("pi", seq_len(n))
  cov_df <- data.frame(age = rnorm(n), row.names = rownames(dat))
  strata <- rep(1:6, each = n / 6)
  cluster <- rep(1:(n / 5), each = 5)
  design_df <- data.frame(
    y = dat$y,
    age = cov_df$age,
    w = runif(n, 0.5, 2.5),
    strata = strata,
    cluster = cluster,
    .analysis_id = rownames(dat),
    row.names = rownames(dat)
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
    analysis_id = rownames(dat),
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
