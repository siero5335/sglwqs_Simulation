validation_groups <- function() {
  list(
    PCBs = c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_187"),
    Metals = c("cd", "pb", "hg", "se", "mg"),
    PFASs = c("PFOA", "PFNA", "PFOS")
  )
}

validation_exposures <- function() {
  unlist(validation_groups(), use.names = FALSE)
}

validation_covariates <- function() {
  c("sex", "age", "bmi")
}

make_block_cor <- function(groups, within_r, between_r = 0) {
  p <- sum(lengths(groups))
  cor <- matrix(between_r, nrow = p, ncol = p)
  diag(cor) <- 1
  idx <- 0L
  for (g in groups) {
    ind <- seq.int(idx + 1L, idx + length(g))
    cor[ind, ind] <- within_r
    diag(cor)[ind] <- 1
    idx <- idx + length(g)
  }
  cor <- as.matrix(Matrix::nearPD(cor, corr = TRUE)$mat)
  dimnames(cor) <- list(validation_exposures(), validation_exposures())
  cor
}

truth_beta <- function(effect = c("global_null", "partial_null", "gaussian_null")) {
  effect <- match.arg(effect)
  vars <- validation_exposures()
  beta <- setNames(rep(0, length(vars)), vars)
  if (identical(effect, "partial_null")) {
    beta[c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_187")] <-
      c(0.04, 0.15, 0.20, 0.08, 0.05)
    beta[c("cd", "pb", "hg", "se", "mg")] <-
      c(-0.12, -0.18, 0.10, 0.14, -0.02)
  }
  beta
}

truth_table <- function(effect = "global_null") {
  groups <- validation_groups()
  beta <- truth_beta(effect)
  out <- do.call(rbind, lapply(names(groups), function(group) {
    vars <- groups[[group]]
    data.frame(
      group = group,
      variable = vars,
      true_beta = unname(beta[vars]),
      stringsAsFactors = FALSE
    )
  }))
  out$true_direction <- ifelse(out$true_beta > 0, "positive",
    ifelse(out$true_beta < 0, "negative", "none")
  )
  out
}

truth_group_direction <- function(effect = "global_null") {
  groups <- validation_groups()
  beta <- truth_beta(effect)
  do.call(rbind, lapply(names(groups), function(group) {
    b <- beta[groups[[group]]]
    data.frame(
      group = group,
      direction = c("positive", "negative"),
      is_true_null = c(!any(b > 0), !any(b < 0)),
      active_variable_count = c(sum(b > 0), sum(b < 0)),
      true_beta_sum = c(sum(b[b > 0]), sum(b[b < 0])),
      stringsAsFactors = FALSE
    )
  }))
}

generate_validation_data <- function(n, seed, effect = "global_null",
                                     outcome_family = c("binomial", "gaussian"),
                                     within_r = 0.6, between_r = 0) {
  outcome_family <- match.arg(outcome_family)
  set.seed(seed)
  groups <- validation_groups()
  vars <- validation_exposures()

  sex <- stats::rbinom(n, 1, 0.5)
  age <- pmax(pmin(round(stats::rnorm(n, 50, 12)), 80), 20)
  bmi <- pmax(pmin(round(22 + 0.05 * age + 1.5 * sex + stats::rnorm(n, 0, 4), 1), 45), 15)

  z <- MASS::mvrnorm(n, mu = rep(0, length(vars)), Sigma = make_block_cor(groups, within_r, between_r))
  colnames(z) <- vars

  pcb_sds <- c(1.5, 1.5, 1.8, 2.0, 1.5)
  pcb_mu <- c(3.5, 3.5, 4.0, 4.5, 3.0)
  pcb <- sweep(z[, groups$PCBs, drop = FALSE], 2, pcb_sds, `*`)
  pcb <- sweep(pcb, 2, pcb_mu, `+`)

  metal_sd <- c(0.6, 0.8, 0.5, 0.3, 0.25)
  metal_mu <- log(c(0.5, 2, 1, 100, 20))
  metal_log <- sweep(z[, groups$Metals, drop = FALSE], 2, metal_sd, `*`)
  metal_log <- sweep(metal_log, 2, metal_mu, `+`)
  metal <- exp(metal_log)

  pfas_sd <- sqrt(c(2, 1, 3.5))
  pfas_mu <- c(3, 1.5, 5)
  pfas <- sweep(z[, groups$PFASs, drop = FALSE], 2, pfas_sd, `*`)
  pfas <- pmax(sweep(pfas, 2, pfas_mu, `+`), 0.1)

  x <- data.frame(pcb, metal, pfas, check.names = FALSE)
  x_scaled <- scale(as.matrix(x[, vars, drop = FALSE]))
  beta <- truth_beta(effect)
  eta_exposure <- as.numeric(x_scaled %*% beta)
  eta_cov <- sex * 0.30 + (age - 50) / 10 * 0.02 + (bmi - 25) / 5 * 0.05

  if (identical(outcome_family, "binomial")) {
    eta <- -3.0 + eta_exposure + eta_cov
    y <- stats::rbinom(n, 1, stats::plogis(eta))
  } else {
    eta <- eta_exposure + eta_cov
    y <- as.numeric(eta + stats::rnorm(n, 0, 1))
  }

  dat <- data.frame(x, sex = sex, age = age, bmi = bmi, Y = y, check.names = FALSE)
  attr(dat, "true_beta") <- beta
  attr(dat, "truth_table") <- truth_table(effect)
  attr(dat, "truth_group_direction") <- truth_group_direction(effect)
  attr(dat, "within_r") <- within_r
  attr(dat, "between_r") <- between_r
  attr(dat, "effect") <- effect
  attr(dat, "outcome_family") <- outcome_family
  dat
}

make_validation_scenarios <- function(profile) {
  global <- expand.grid(
    n = c(500L, 5000L),
    within_r = c(0.20, 0.95),
    between_r = 0,
    KEEP.OUT.ATTRS = FALSE
  )
  global$scenario_family <- "global_binary"
  global$effect <- "global_null"
  global$outcome_family <- "binomial"

  partial <- expand.grid(
    n = c(500L, 5000L),
    between_r = c(0, 0.30),
    KEEP.OUT.ATTRS = FALSE
  )
  partial$within_r <- 0.60
  partial$scenario_family <- "partial_binary"
  partial$effect <- "partial_null"
  partial$outcome_family <- "binomial"

  gaussian <- data.frame(
    n = 1000L,
    within_r = c(0.20, 0.95),
    between_r = 0,
    scenario_family = "gaussian_null",
    effect = "gaussian_null",
    outcome_family = "gaussian",
    stringsAsFactors = FALSE
  )

  scenarios <- rbind(
    global[, names(gaussian)],
    partial[, names(gaussian)],
    gaussian
  )
  scenarios$scenario_id <- sprintf(
    "%s_n%s_w%s_b%s",
    scenarios$scenario_family,
    scenarios$n,
    gsub("\\.", "p", scenarios$within_r),
    gsub("\\.", "p", scenarios$between_r)
  )
  scenarios$scenario_label <- sprintf(
    "%s; n=%s; within=%s; between=%s",
    scenarios$scenario_family,
    scenarios$n,
    scenarios$within_r,
    scenarios$between_r
  )
  rownames(scenarios) <- NULL
  scenarios
}

make_replicate_jobs <- function(profile) {
  scenarios <- make_validation_scenarios(profile)
  jobs <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    sc <- scenarios[i, , drop = FALSE]
    n_rep <- profile_reps(profile, sc$scenario_family)
    data.frame(
      sc[rep(1, n_rep), , drop = FALSE],
      replicate = seq_len(n_rep),
      data_seed = as.integer(profile$seed_start) + i * 100000L + seq_len(n_rep),
      fit_seed = as.integer(profile$seed_start) + i * 100000L + 50000L + seq_len(n_rep),
      job_type = "replicate",
      stringsAsFactors = FALSE
    )
  }))
  jobs$job_id <- sprintf("%s_rep%04d", jobs$scenario_id, jobs$replicate)
  rownames(jobs) <- NULL
  jobs
}

make_split_stability_jobs <- function(profile) {
  base <- rbind(
    data.frame(
      scenario_family = "split_global_binary",
      effect = "global_null",
      outcome_family = "binomial",
      n = c(500L, 5000L),
      within_r = 0.20,
      between_r = 0,
      data_seed = as.integer(profile$seed_start) + c(9001L, 9002L),
      stringsAsFactors = FALSE
    ),
    data.frame(
      scenario_family = "split_partial_binary",
      effect = "partial_null",
      outcome_family = "binomial",
      n = c(500L, 5000L),
      within_r = 0.60,
      between_r = 0.30,
      data_seed = as.integer(profile$seed_start) + c(9101L, 9102L),
      stringsAsFactors = FALSE
    )
  )
  split_reps <- as.integer(profile$split_reps)
  jobs <- do.call(rbind, lapply(seq_len(nrow(base)), function(i) {
    sc <- base[i, , drop = FALSE]
    out <- data.frame(
      sc[rep(1, split_reps), , drop = FALSE],
      replicate = seq_len(split_reps),
      fit_seed = as.integer(profile$seed_start) + i * 200000L + 70000L + seq_len(split_reps),
      job_type = "split_stability",
      stringsAsFactors = FALSE
    )
    out$scenario_id <- sprintf(
      "%s_n%s_w%s_b%s",
      out$scenario_family,
      out$n,
      gsub("\\.", "p", out$within_r),
      gsub("\\.", "p", out$between_r)
    )
    out$scenario_label <- sprintf(
      "%s fixed data; n=%s; within=%s; between=%s",
      out$scenario_family,
      out$n,
      out$within_r,
      out$between_r
    )
    out
  }))
  jobs$job_id <- sprintf("%s_split%04d", jobs$scenario_id, jobs$replicate)
  rownames(jobs) <- NULL
  jobs
}
