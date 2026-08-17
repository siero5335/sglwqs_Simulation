legacy13_groups <- function() {
  list(
    PCBs = c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_187"),
    Metals = c("cd", "pb", "hg", "se", "mg"),
    PFASs = c("PFOA", "PFNA", "PFOS")
  )
}

legacy13_exposures <- function() {
  unlist(legacy13_groups(), use.names = FALSE)
}

legacy13_covariates <- function() {
  c("sex", "age", "bmi")
}

legacy13_base_beta <- function() {
  c(
    pcb_118 = 0.04, pcb_138 = 0.15, pcb_153 = 0.20, pcb_180 = 0.08, pcb_187 = 0.05,
    cd = -0.12, pb = -0.18, hg = 0.10, se = 0.14, mg = -0.02,
    PFOA = 0, PFNA = 0, PFOS = 0
  )
}

legacy13_beta <- function(signal_profile = c("baseline", "weak", "global_null")) {
  signal_profile <- match.arg(signal_profile)
  beta <- legacy13_base_beta()
  if (identical(signal_profile, "weak")) {
    beta <- beta * 0.5
  }
  if (identical(signal_profile, "global_null")) {
    beta[] <- 0
  }
  beta
}

r2r2_truth_tolerance <- function() 1e-12

truth_from_beta <- function(beta, groups, active_threshold = r2r2_truth_tolerance()) {
  if (length(active_threshold) != 1L || !is.finite(active_threshold) || active_threshold < 0) {
    stop("active_threshold must be one finite non-negative value.", call. = FALSE)
  }
  vars <- unlist(groups, use.names = FALSE)
  out <- do.call(rbind, lapply(names(groups), function(group) {
    data.frame(
      Variable = groups[[group]],
      Group = group,
      True_Beta = as.numeric(beta[groups[[group]]]),
      stringsAsFactors = FALSE
    )
  }))
  out$True_Direction <- ifelse(out$True_Beta > active_threshold, "Positive",
    ifelse(out$True_Beta < -active_threshold, "Negative", "None")
  )
  out$IsActive <- out$True_Direction != "None"
  out$IsActiveGroup <- ave(out$IsActive, out$Group, FUN = function(x) any(x))
  out$Group_Size <- as.integer(lengths(groups)[out$Group])
  out$Group_Size_Tier <- "legacy13"
  out$Effect_Tier <- ifelse(out$IsActive, "active", "null")
  out$Group_Type <- ave(out$True_Direction, out$Group, FUN = function(x) {
    active <- unique(x[x != "None"])
    if (!length(active)) "null" else if (length(active) == 1L) "unidirectional" else "mixed-direction"
  })
  out[match(vars, out$Variable), , drop = FALSE]
}

legacy13_truth <- function(signal_profile = "baseline") {
  truth_from_beta(legacy13_beta(signal_profile), legacy13_groups())
}

generate_legacy13_data <- function(n,
                                   data_seed,
                                   family = c("binomial", "gaussian"),
                                   signal_profile = c("baseline", "weak", "global_null"),
                                   gaussian_sd = 1) {
  family <- match.arg(family)
  signal_profile <- match.arg(signal_profile)
  set.seed(as.integer(data_seed))

  sex <- stats::rbinom(n, 1, 0.5)
  age <- pmax(pmin(round(stats::rnorm(n, 50, 12)), 80), 20)
  bmi <- pmax(pmin(round(22 + 0.05 * age + 1.5 * sex + stats::rnorm(n, 0, 4), 1), 45), 15)

  cor_pcb <- matrix(c(
    1.00, 0.60, 0.70, 0.65, 0.60,
    0.60, 1.00, 0.75, 0.70, 0.65,
    0.70, 0.75, 1.00, 0.80, 0.70,
    0.65, 0.70, 0.80, 1.00, 0.75,
    0.60, 0.65, 0.70, 0.75, 1.00
  ), 5, 5)
  pcb <- MASS::mvrnorm(
    n,
    c(3.5, 3.5, 4.0, 4.5, 3.0),
    diag(c(1.5, 1.5, 1.8, 2.0, 1.5)) %*% cor_pcb %*% diag(c(1.5, 1.5, 1.8, 2.0, 1.5))
  )
  colnames(pcb) <- legacy13_groups()$PCBs

  cor_metal <- matrix(c(
    1.00, 0.60, 0.15, 0.10, 0.20,
    0.60, 1.00, 0.10, 0.15, 0.25,
    0.15, 0.10, 1.00, 0.55, 0.20,
    0.10, 0.15, 0.55, 1.00, 0.30,
    0.20, 0.25, 0.20, 0.30, 1.00
  ), 5, 5, byrow = TRUE)
  metal_raw <- MASS::mvrnorm(
    n,
    log(c(0.5, 2, 1, 100, 20)),
    diag(c(0.6, 0.8, 0.5, 0.3, 0.25)) %*% cor_metal %*%
      diag(c(0.6, 0.8, 0.5, 0.3, 0.25))
  )
  metal <- exp(metal_raw)
  colnames(metal) <- legacy13_groups()$Metals

  pfas <- pmax(MASS::mvrnorm(n, c(3, 1.5, 5), diag(c(2, 1, 3.5))), 0.1)
  colnames(pfas) <- legacy13_groups()$PFASs

  groups <- legacy13_groups()
  vars <- legacy13_exposures()
  x <- data.frame(pcb, metal, pfas, check.names = FALSE)
  beta <- legacy13_beta(signal_profile)
  x_scaled <- scale(as.matrix(x[, vars, drop = FALSE]))
  eta_exposure <- as.numeric(x_scaled %*% beta)
  eta_cov <- sex * 0.3 + (age - 50) / 10 * 0.02 + (bmi - 25) / 5 * 0.05
  eta_linear <- eta_exposure + eta_cov

  if (identical(family, "binomial")) {
    eta <- -3.0 + eta_linear
    y <- stats::rbinom(n, 1, stats::plogis(eta))
  } else {
    eta <- eta_linear
    y <- as.numeric(eta + stats::rnorm(n, 0, gaussian_sd))
  }

  dat <- data.frame(x, sex = sex, age = age, bmi = bmi, Y = y, check.names = FALSE)
  attr(dat, "groups") <- groups
  attr(dat, "truth") <- truth_from_beta(beta, groups)
  attr(dat, "beta") <- beta
  attr(dat, "family") <- family
  attr(dat, "signal_profile") <- signal_profile
  attr(dat, "gaussian_residual_sd") <- if (identical(family, "gaussian")) gaussian_sd else NA_real_
  attr(dat, "eta_linear") <- eta_linear
  dat
}

block_cor_from_groups <- function(groups, within_r, between_r = 0) {
  p <- sum(lengths(groups))
  cor <- matrix(between_r, nrow = p, ncol = p)
  diag(cor) <- 1
  offset <- 0L
  for (g in groups) {
    idx <- seq.int(offset + 1L, offset + length(g))
    cor[idx, idx] <- within_r
    diag(cor)[idx] <- 1
    offset <- offset + length(g)
  }
  cor <- as.matrix(Matrix::nearPD(cor, corr = TRUE)$mat)
  dimnames(cor) <- list(unlist(groups, use.names = FALSE), unlist(groups, use.names = FALSE))
  cor
}

generate_validation_split_data <- function(n,
                                           data_seed,
                                           family = c("binomial", "gaussian"),
                                           effect_profile = c("global_null", "partial_null"),
                                           within_r = 0.6,
                                           between_r = 0,
                                           gaussian_sd = 1) {
  family <- match.arg(family)
  effect_profile <- match.arg(effect_profile)
  signal_profile <- if (identical(effect_profile, "global_null")) "global_null" else "baseline"
  set.seed(as.integer(data_seed))

  groups <- legacy13_groups()
  vars <- legacy13_exposures()
  beta <- legacy13_beta(signal_profile)

  sex <- stats::rbinom(n, 1, 0.5)
  age <- pmax(pmin(round(stats::rnorm(n, 50, 12)), 80), 20)
  bmi <- pmax(pmin(round(22 + 0.05 * age + 1.5 * sex + stats::rnorm(n, 0, 4), 1), 45), 15)

  z <- MASS::mvrnorm(n, mu = rep(0, length(vars)), Sigma = block_cor_from_groups(groups, within_r, between_r))
  colnames(z) <- vars

  pcb <- sweep(z[, groups$PCBs, drop = FALSE], 2, c(1.5, 1.5, 1.8, 2.0, 1.5), `*`)
  pcb <- sweep(pcb, 2, c(3.5, 3.5, 4.0, 4.5, 3.0), `+`)
  metal_log <- sweep(z[, groups$Metals, drop = FALSE], 2, c(0.6, 0.8, 0.5, 0.3, 0.25), `*`)
  metal_log <- sweep(metal_log, 2, log(c(0.5, 2, 1, 100, 20)), `+`)
  metal <- exp(metal_log)
  pfas <- sweep(z[, groups$PFASs, drop = FALSE], 2, sqrt(c(2, 1, 3.5)), `*`)
  pfas <- pmax(sweep(pfas, 2, c(3, 1.5, 5), `+`), 0.1)

  x <- data.frame(pcb, metal, pfas, check.names = FALSE)
  x_scaled <- scale(as.matrix(x[, vars, drop = FALSE]))
  eta_exposure <- as.numeric(x_scaled %*% beta)
  eta_cov <- sex * 0.30 + (age - 50) / 10 * 0.02 + (bmi - 25) / 5 * 0.05
  eta_linear <- eta_exposure + eta_cov
  if (identical(family, "binomial")) {
    y <- stats::rbinom(n, 1, stats::plogis(-3.0 + eta_linear))
  } else {
    y <- as.numeric(eta_linear + stats::rnorm(n, 0, gaussian_sd))
  }

  dat <- data.frame(x, sex = sex, age = age, bmi = bmi, Y = y, check.names = FALSE)
  attr(dat, "groups") <- groups
  attr(dat, "truth") <- truth_from_beta(beta, groups)
  attr(dat, "beta") <- beta
  attr(dat, "family") <- family
  attr(dat, "effect_profile") <- effect_profile
  attr(dat, "within_r") <- within_r
  attr(dat, "between_r") <- between_r
  attr(dat, "gaussian_residual_sd") <- if (identical(family, "gaussian")) gaussian_sd else NA_real_
  attr(dat, "eta_linear") <- eta_linear
  dat
}
