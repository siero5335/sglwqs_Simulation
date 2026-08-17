dense_group_signal_observed_rho <- function(within_rho = 0.45,
                                            cross_rho = 0.05,
                                            log_sd = 0.45) {
  latent_rho <- within_rho + cross_rho
  assert_or_stop(latent_rho >= 0 && latent_rho < 1,
    "The within-group latent correlation must be in [0, 1)."
  )
  (exp(log_sd^2 * latent_rho) - 1) / (exp(log_sd^2) - 1)
}

equicorrelation_signal_variance <- function(beta, rho) {
  (1 - rho) * sum(beta^2) + rho * sum(beta)^2
}

dense_mixed_signs <- function(n_components) {
  n_positive <- ceiling(2 * n_components / 3)
  c(rep(1, n_positive), rep(-1, n_components - n_positive))
}

make_dense_group_signal_truth <- function(within_rho = 0.45,
                                          cross_rho = 0.05) {
  sparse <- make_factorial_truth("unbalanced", "weak_heterogeneous")
  groups <- sparse$groups
  beta <- setNames(rep(0, length(sparse$vars)), sparse$vars)
  observed_rho <- dense_group_signal_observed_rho(within_rho, cross_rho)
  active_plan <- data.frame(
    Group = c("G01", "G05", "G10"),
    Direction = c("positive", "negative", "mixed"),
    stringsAsFactors = FALSE
  )
  signal_rows <- vector("list", nrow(active_plan))

  for (i in seq_len(nrow(active_plan))) {
    group <- active_plan$Group[[i]]
    group_vars <- groups[[group]]
    sparse_beta <- sparse$beta[group_vars]
    target_variance <- equicorrelation_signal_variance(sparse_beta, observed_rho)
    signs <- switch(active_plan$Direction[[i]],
      positive = rep(1, length(group_vars)),
      negative = rep(-1, length(group_vars)),
      mixed = dense_mixed_signs(length(group_vars))
    )
    unit_variance <- equicorrelation_signal_variance(signs, observed_rho)
    magnitude <- sqrt(target_variance / unit_variance)
    beta[group_vars] <- signs * magnitude
    signal_rows[[i]] <- data.frame(
      Group = group,
      Group_Direction = active_plan$Direction[[i]],
      Group_Size = length(group_vars),
      Sparse_Active_Components = sum(abs(sparse_beta) > r2r2_truth_tolerance()),
      Dense_Active_Components = length(group_vars),
      Target_Signal_Variance = target_variance,
      Dense_Signal_Variance = equicorrelation_signal_variance(beta[group_vars], observed_rho),
      Component_Abs_Beta = magnitude,
      stringsAsFactors = FALSE
    )
  }

  truth <- truth_from_beta(beta, groups)
  active_groups <- active_plan$Group
  truth$Group_Size_Tier <- ifelse(truth$Group == "G01", "small",
    ifelse(truth$Group == "G05", "medium",
      ifelse(truth$Group == "G10", "large", "null")
    )
  )
  truth$Effect_Tier <- ifelse(truth$Group %in% active_groups, "diffuse", "null")
  truth$Effect_Profile <- "dense_group_matched"
  truth$Group_Structure <- "unbalanced"
  signal_spec <- do.call(rbind, signal_rows)

  assert_or_stop(sum(lengths(groups)) == 100L, "Dense group sizes do not sum to p=100.")
  assert_or_stop(sum(truth$IsActive) == 30L, "Dense truth must contain exactly 30 active components.")
  assert_or_stop(all(beta[unlist(groups[setdiff(names(groups), active_groups)], use.names = FALSE)] == 0),
    "Null groups must have exactly zero coefficients."
  )
  assert_or_stop(all(abs(signal_spec$Target_Signal_Variance - signal_spec$Dense_Signal_Variance) < 1e-12),
    "Dense and sparse group signal variances do not match."
  )

  list(
    groups = groups,
    vars = sparse$vars,
    beta = beta,
    truth = truth,
    signal_spec = signal_spec,
    observed_within_rho = observed_rho,
    sparse_reference = sparse
  )
}

generate_dense_group_signal_data <- function(n,
                                             data_seed,
                                             family = c("binomial", "gaussian"),
                                             within_rho = 0.45,
                                             cross_rho = 0.05,
                                             gaussian_sd = 1) {
  family <- match.arg(family)
  spec <- make_dense_group_signal_truth(within_rho, cross_rho)
  x <- generate_latent_exposure_matrix(n, spec$groups, data_seed, within_rho, cross_rho)
  covars <- make_covariates(n, data_seed)
  eta_exposure <- as.numeric(scale(x) %*% spec$beta)
  eta_cov <- covars$sex * 0.3 + (covars$age - 50) / 10 * 0.02 + (covars$bmi - 25) / 5 * 0.05
  eta_linear <- eta_exposure + eta_cov
  set.seed(as.integer(data_seed) + 39001L)
  y <- if (identical(family, "binomial")) {
    stats::rbinom(n, 1, stats::plogis(-3.0 + eta_linear))
  } else {
    as.numeric(eta_linear + stats::rnorm(n, 0, gaussian_sd))
  }
  dat <- data.frame(x, covars, Y = y, check.names = FALSE)
  attr(dat, "groups") <- spec$groups
  attr(dat, "truth") <- spec$truth
  attr(dat, "beta") <- spec$beta
  attr(dat, "family") <- family
  attr(dat, "group_structure") <- "unbalanced"
  attr(dat, "effect_profile") <- "dense_group_matched"
  attr(dat, "signal_profile") <- "dense_group_matched"
  attr(dat, "within_rho") <- within_rho
  attr(dat, "cross_rho") <- cross_rho
  attr(dat, "gaussian_residual_sd") <- if (identical(family, "gaussian")) gaussian_sd else NA_real_
  attr(dat, "eta_linear") <- eta_linear
  attr(dat, "signal_spec") <- spec$signal_spec
  dat
}
