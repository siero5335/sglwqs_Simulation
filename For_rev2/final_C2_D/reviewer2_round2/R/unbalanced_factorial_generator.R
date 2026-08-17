make_group_definitions_equal <- function(p, n_groups = 10L) {
  stopifnot(p >= n_groups)
  vars <- sprintf("X%03d", seq_len(p))
  group_ids <- cut(seq_len(p),
    breaks = n_groups,
    labels = sprintf("G%02d", seq_len(n_groups)),
    include.lowest = TRUE
  )
  split(vars, group_ids)
}

make_group_definitions_from_sizes <- function(group_sizes) {
  group_sizes <- as.integer(group_sizes)
  stopifnot(all(group_sizes > 0L))
  vars <- sprintf("X%03d", seq_len(sum(group_sizes)))
  split(vars, rep(sprintf("G%02d", seq_along(group_sizes)), group_sizes))
}

factorial_group_sizes <- function(group_structure = c("balanced", "unbalanced")) {
  group_structure <- match.arg(group_structure)
  if (identical(group_structure, "balanced")) {
    rep(10L, 10L)
  } else {
    c(4L, 6L, 8L, 8L, 10L, 10L, 12L, 12L, 14L, 16L)
  }
}

active_group_plan <- function(group_structure = c("balanced", "unbalanced")) {
  group_structure <- match.arg(group_structure)
  if (identical(group_structure, "unbalanced")) {
    data.frame(
      Group = c("G01", "G05", "G10"),
      Group_Size_Tier = c("small", "medium", "large"),
      Group_Direction = c("positive", "negative", "mixed"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Group = c("G01", "G05", "G10"),
      Group_Size_Tier = c("small", "medium", "large"),
      Group_Direction = c("positive", "negative", "mixed"),
      stringsAsFactors = FALSE
    )
  }
}

effect_profile_abs_beta <- function(effect_profile = c("uniform_baseline", "heterogeneous", "weak_heterogeneous")) {
  effect_profile <- match.arg(effect_profile)
  if (identical(effect_profile, "uniform_baseline")) {
    c(strong = 0.20, medium = 0.20, weak = 0.20)
  } else if (identical(effect_profile, "heterogeneous")) {
    c(strong = 0.30, medium = 0.17, weak = 0.05)
  } else {
    0.5 * c(strong = 0.30, medium = 0.17, weak = 0.05)
  }
}

make_factorial_truth <- function(group_structure = c("balanced", "unbalanced"),
                                 effect_profile = c("uniform_baseline", "heterogeneous", "weak_heterogeneous")) {
  group_structure <- match.arg(group_structure)
  effect_profile <- match.arg(effect_profile)
  groups <- make_group_definitions_from_sizes(factorial_group_sizes(group_structure))
  stopifnot(sum(lengths(groups)) == 100L)
  plan <- active_group_plan(group_structure)
  abs_beta <- effect_profile_abs_beta(effect_profile)
  vars <- unlist(groups, use.names = FALSE)
  beta <- setNames(rep(0, length(vars)), vars)
  effect_tier <- setNames(rep("null", length(vars)), vars)
  group_size_tier <- setNames(rep("null", length(vars)), vars)
  group_type <- setNames(rep("null", length(vars)), vars)

  for (i in seq_len(nrow(plan))) {
    group <- plan$Group[[i]]
    group_vars <- groups[[group]]
    active_vars <- group_vars[seq_len(3L)]
    group_size_tier[group_vars] <- plan$Group_Size_Tier[[i]]
    group_type[group_vars] <- if (identical(plan$Group_Direction[[i]], "mixed")) {
      "mixed-direction"
    } else {
      "unidirectional"
    }
    effect_tier[active_vars] <- names(abs_beta)

    signed <- switch(plan$Group_Direction[[i]],
      positive = abs_beta,
      negative = -abs_beta,
      mixed = c(strong = abs_beta[["strong"]], medium = -abs_beta[["medium"]], weak = abs_beta[["weak"]])
    )
    beta[active_vars] <- as.numeric(signed)
  }

  truth <- truth_from_beta(beta, groups, active_threshold = 0)
  truth$Group_Size <- as.integer(lengths(groups)[truth$Group])
  truth$Group_Size_Tier <- unname(group_size_tier[truth$Variable])
  truth$Effect_Tier <- unname(effect_tier[truth$Variable])
  truth$Group_Type <- unname(group_type[truth$Variable])
  truth$Effect_Profile <- effect_profile
  truth$Group_Structure <- group_structure
  list(groups = groups, vars = vars, beta = beta, truth = truth)
}

make_highdim_p_truth <- function(p,
                                 n_groups = 10L,
                                 signal_profile = c("baseline", "weak")) {
  signal_profile <- match.arg(signal_profile)
  groups <- make_group_definitions_equal(p, n_groups)
  vars <- unlist(groups, use.names = FALSE)
  beta <- setNames(rep(0, length(vars)), vars)

  g1 <- groups[[1]]
  n_g1 <- min(length(g1), max(3L, ceiling(length(g1) * 0.15)))
  beta[g1[seq_len(n_g1)]] <- seq(0.08, 0.12, length.out = n_g1)

  g2 <- groups[[2]]
  n_g2_neg <- min(length(g2), max(2L, floor(length(g2) * 0.08)))
  n_g2_pos <- min(max(length(g2) - n_g2_neg, 0L), max(2L, floor(length(g2) * 0.06)))
  if (n_g2_neg > 0L) {
    beta[g2[seq_len(n_g2_neg)]] <- seq(-0.16, -0.11, length.out = n_g2_neg)
  }
  if (n_g2_pos > 0L) {
    pos_idx <- seq.int(n_g2_neg + 1L, n_g2_neg + n_g2_pos)
    beta[g2[pos_idx]] <- seq(0.10, 0.14, length.out = n_g2_pos)
  }

  g3 <- groups[[3]]
  n_g3 <- min(length(g3), max(2L, ceiling(length(g3) * 0.08)))
  beta[g3[seq_len(n_g3)]] <- seq(-0.14, -0.10, length.out = n_g3)

  if (identical(signal_profile, "weak")) {
    beta <- beta * 0.5
  }
  truth <- truth_from_beta(beta, groups, active_threshold = 0.01)
  truth$Group_Size_Tier <- "equal"
  truth$Effect_Tier <- ifelse(abs(truth$True_Beta) >= 0.14, "strong",
    ifelse(abs(truth$True_Beta) >= 0.08, "medium",
      ifelse(truth$IsActive, "weak", "null")
    )
  )
  truth$Group_Type <- ave(truth$True_Direction, truth$Group, FUN = function(x) {
    active <- unique(x[x != "None"])
    if (!length(active)) "null" else if (length(active) == 1L) "unidirectional" else "mixed-direction"
  })
  list(groups = groups, vars = vars, beta = beta, truth = truth)
}

generate_latent_exposure_matrix <- function(n,
                                            groups,
                                            seed,
                                            within_rho = 0.45,
                                            cross_rho = 0.05) {
  set.seed(as.integer(seed))
  vars <- unlist(groups, use.names = FALSE)
  if (cross_rho < 0 || within_rho < cross_rho || within_rho >= 1) {
    stop("Require 0 <= cross_rho <= within_rho < 1", call. = FALSE)
  }
  global_factor <- stats::rnorm(n)
  group_factors <- matrix(stats::rnorm(n * length(groups)), nrow = n, ncol = length(groups))
  colnames(group_factors) <- names(groups)

  xz <- matrix(NA_real_, nrow = n, ncol = length(vars), dimnames = list(NULL, vars))
  for (g in seq_along(groups)) {
    group_vars <- groups[[g]]
    idx <- match(group_vars, vars)
    eps <- matrix(stats::rnorm(n * length(group_vars)), nrow = n, ncol = length(group_vars))
    # `within_rho` is the target total within-group correlation. The global
    # factor contributes `cross_rho` to every pair, so the group-specific
    # contribution is their difference.
    latent <- sqrt(cross_rho) * global_factor +
      sqrt(within_rho - cross_rho) * group_factors[, g] +
      sqrt(1 - within_rho) * eps
    xz[, idx] <- latent
  }
  group_shift <- rep(seq(-0.3, 0.3, length.out = length(groups)), lengths(groups))
  x <- exp(0.45 * sweep(xz, 2, group_shift, "+"))
  colnames(x) <- vars
  x
}

make_covariates <- function(n, seed) {
  set.seed(as.integer(seed) + 19001L)
  sex <- stats::rbinom(n, 1, 0.5)
  age <- pmax(pmin(round(stats::rnorm(n, 50, 12)), 80), 20)
  bmi <- pmax(pmin(round(22 + 0.05 * age + 1.5 * sex + stats::rnorm(n, 0, 4), 1), 45), 15)
  data.frame(sex = sex, age = age, bmi = bmi)
}

factorial_binomial_target_prevalence <- function() {
  0.20
}

calibrate_logistic_intercept <- function(eta_linear,
                                         target_prevalence = factorial_binomial_target_prevalence()) {
  target_prevalence <- as.numeric(target_prevalence)
  if (length(target_prevalence) != 1L || !is.finite(target_prevalence) ||
      target_prevalence <= 0 || target_prevalence >= 1) {
    stop("target_prevalence must be one finite value strictly between 0 and 1.", call. = FALSE)
  }
  if (!length(eta_linear) || any(!is.finite(eta_linear))) {
    stop("eta_linear must contain only finite values.", call. = FALSE)
  }
  objective <- function(intercept) {
    mean(stats::plogis(intercept + eta_linear)) - target_prevalence
  }
  stats::uniroot(objective, interval = c(-30, 30), tol = 1e-12)$root
}

generate_highdim_p_data_matched <- function(p,
                                            n,
                                            data_seed,
                                            family = c("binomial", "gaussian"),
                                            signal_profile = c("baseline", "weak"),
                                            n_groups = 10L,
                                            within_rho = 0.45,
                                            cross_rho = 0.05,
                                            gaussian_sd = 1) {
  family <- match.arg(family)
  signal_profile <- match.arg(signal_profile)
  spec <- make_highdim_p_truth(p = p, n_groups = n_groups, signal_profile = signal_profile)
  x <- generate_latent_exposure_matrix(n, spec$groups, data_seed, within_rho, cross_rho)
  covars <- make_covariates(n, data_seed)
  x_scaled <- scale(x)
  eta_exposure <- as.numeric(x_scaled %*% spec$beta)
  eta_cov <- covars$sex * 0.3 + (covars$age - 50) / 10 * 0.02 + (covars$bmi - 25) / 5 * 0.05
  eta_linear <- eta_exposure + eta_cov
  set.seed(as.integer(data_seed) + 29001L)
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
  attr(dat, "signal_profile") <- signal_profile
  attr(dat, "within_rho") <- within_rho
  attr(dat, "cross_rho") <- cross_rho
  attr(dat, "gaussian_residual_sd") <- if (identical(family, "gaussian")) gaussian_sd else NA_real_
  attr(dat, "eta_linear") <- eta_linear
  dat
}

generate_factorial_data <- function(n,
                                    data_seed,
                                    family = c("binomial", "gaussian"),
                                    group_structure = c("balanced", "unbalanced"),
                                    effect_profile = c("uniform_baseline", "heterogeneous", "weak_heterogeneous"),
                                    within_rho = 0.45,
                                    cross_rho = 0.05,
                                    gaussian_sd = 1,
                                    binomial_target_prevalence = factorial_binomial_target_prevalence()) {
  family <- match.arg(family)
  group_structure <- match.arg(group_structure)
  effect_profile <- match.arg(effect_profile)
  spec <- make_factorial_truth(group_structure, effect_profile)
  assert_or_stop(sum(lengths(spec$groups)) == 100L, "Factorial group sizes do not sum to p=100.")
  x <- generate_latent_exposure_matrix(n, spec$groups, data_seed, within_rho, cross_rho)
  covars <- make_covariates(n, data_seed)
  eta_exposure <- as.numeric(scale(x) %*% spec$beta)
  eta_cov <- covars$sex * 0.3 + (covars$age - 50) / 10 * 0.02 + (covars$bmi - 25) / 5 * 0.05
  eta_linear <- eta_exposure + eta_cov
  binomial_intercept <- NA_real_
  binomial_expected_prevalence <- NA_real_
  set.seed(as.integer(data_seed) + 39001L)
  y <- if (identical(family, "binomial")) {
    binomial_intercept <- calibrate_logistic_intercept(
      eta_linear,
      target_prevalence = binomial_target_prevalence
    )
    event_probability <- stats::plogis(binomial_intercept + eta_linear)
    binomial_expected_prevalence <- mean(event_probability)
    stats::rbinom(n, 1, event_probability)
  } else {
    as.numeric(eta_linear + stats::rnorm(n, 0, gaussian_sd))
  }
  dat <- data.frame(x, covars, Y = y, check.names = FALSE)
  attr(dat, "groups") <- spec$groups
  attr(dat, "truth") <- spec$truth
  attr(dat, "beta") <- spec$beta
  attr(dat, "family") <- family
  attr(dat, "group_structure") <- group_structure
  attr(dat, "effect_profile") <- effect_profile
  attr(dat, "within_rho") <- within_rho
  attr(dat, "cross_rho") <- cross_rho
  attr(dat, "gaussian_residual_sd") <- if (identical(family, "gaussian")) gaussian_sd else NA_real_
  attr(dat, "eta_linear") <- eta_linear
  attr(dat, "binomial_intercept") <- binomial_intercept
  attr(dat, "binomial_target_prevalence") <- if (identical(family, "binomial")) {
    as.numeric(binomial_target_prevalence)
  } else {
    NA_real_
  }
  attr(dat, "binomial_expected_prevalence") <- binomial_expected_prevalence
  dat
}
