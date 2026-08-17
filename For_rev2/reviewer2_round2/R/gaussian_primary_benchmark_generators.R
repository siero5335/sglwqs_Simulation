gaussian_correlation_benchmark_beta <- function() {
  c(
    pcb_118 = 0.00, pcb_138 = 0.10, pcb_153 = 0.20, pcb_180 = 0.00, pcb_187 = 0.00,
    cd = -0.08, pb = -0.18, hg = 0.00, se = 0.00, mg = 0.00,
    PFOA = 0.00, PFNA = 0.00, PFOS = 0.00
  )
}

gaussian_active_component_beta <- function() {
  c(
    pcb_118 = 0.00, pcb_138 = 0.00, pcb_153 = 0.25, pcb_180 = 0.00, pcb_187 = 0.00,
    cd = 0.00, pb = -0.20, hg = 0.00, se = 0.00, mg = 0.00,
    PFOA = 0.00, PFNA = 0.00, PFOS = 0.00
  )
}

generate_primary_benchmark_covariates <- function(n) {
  sex <- stats::rbinom(n, 1, 0.5)
  age <- pmax(pmin(round(stats::rnorm(n, 50, 12)), 80), 20)
  bmi <- pmax(pmin(round(22 + 0.05 * age + 1.5 * sex + stats::rnorm(n, 0, 4), 1), 45), 15)
  list(sex = sex, age = age, bmi = bmi)
}

finalize_gaussian_primary_benchmark <- function(x,
                                                covariates,
                                                beta,
                                                gaussian_sd,
                                                design_source,
                                                signal_profile,
                                                within_rho) {
  groups <- legacy13_groups()
  vars <- legacy13_exposures()
  x_scaled <- scale(as.matrix(x[, vars, drop = FALSE]))
  eta_exposure <- as.numeric(x_scaled %*% beta)
  eta_cov <- covariates$sex * 0.30 +
    (covariates$age - 50) / 10 * 0.02 +
    (covariates$bmi - 25) / 5 * 0.05

  # Preserve the original logistic benchmark's complete linear predictor.
  eta_linear <- -3.0 + eta_exposure + eta_cov
  y <- as.numeric(eta_linear + stats::rnorm(nrow(x), 0, gaussian_sd))

  dat <- data.frame(
    x,
    sex = covariates$sex,
    age = covariates$age,
    bmi = covariates$bmi,
    Y = y,
    check.names = FALSE
  )
  attr(dat, "groups") <- groups
  attr(dat, "truth") <- truth_from_beta(beta, groups)
  attr(dat, "beta") <- beta
  attr(dat, "family") <- "gaussian"
  attr(dat, "signal_profile") <- signal_profile
  attr(dat, "effect_profile") <- signal_profile
  attr(dat, "gaussian_residual_sd") <- gaussian_sd
  attr(dat, "eta_linear") <- eta_linear
  attr(dat, "within_rho") <- within_rho
  attr(dat, "cross_rho") <- 0
  attr(dat, "design_source") <- design_source
  dat
}

generate_gaussian_correlation_benchmark_data <- function(n = 10000L,
                                                         data_seed = 71L,
                                                         correlation = 0.5,
                                                         gaussian_sd = 1) {
  if (length(correlation) != 1L || !is.finite(correlation) || correlation < 0 || correlation >= 1) {
    stop("correlation must be one finite value in [0, 1).", call. = FALSE)
  }
  set.seed(as.integer(data_seed))
  covariates <- generate_primary_benchmark_covariates(n)

  cor_pcb <- matrix(correlation, 5, 5)
  diag(cor_pcb) <- 1
  cor_pcb <- as.matrix(Matrix::nearPD(cor_pcb)$mat)
  sd_pcb <- diag(c(1.5, 1.8, 1.6, 1.7, 1.4))
  pcb <- MASS::mvrnorm(
    n,
    mu = c(3.5, 4.0, 3.8, 4.2, 3.6),
    Sigma = sd_pcb %*% cor_pcb %*% sd_pcb
  )
  colnames(pcb) <- legacy13_groups()$PCBs

  cor_metal <- matrix(correlation, 5, 5)
  diag(cor_metal) <- 1
  cor_metal <- as.matrix(Matrix::nearPD(cor_metal)$mat)
  sd_metal <- diag(c(0.6, 0.8, 0.5, 0.3, 0.25))
  metal_raw <- MASS::mvrnorm(
    n,
    mu = log(c(0.5, 2, 1, 100, 20)),
    Sigma = sd_metal %*% cor_metal %*% sd_metal
  )
  metal <- exp(metal_raw)
  colnames(metal) <- legacy13_groups()$Metals

  pfas_rho <- min(correlation, 0.4)
  cor_pfas <- matrix(pfas_rho, 3, 3)
  diag(cor_pfas) <- 1
  pfas <- pmax(
    MASS::mvrnorm(n, mu = c(3, 1.5, 5), Sigma = cor_pfas * 2),
    0.1
  )
  colnames(pfas) <- legacy13_groups()$PFASs

  x <- data.frame(pcb, metal, pfas, check.names = FALSE)
  finalize_gaussian_primary_benchmark(
    x = x,
    covariates = covariates,
    beta = gaussian_correlation_benchmark_beta(),
    gaussian_sd = gaussian_sd,
    design_source = "sglwqs_Simulation-main/compare_methods_correlation_robustness_fixed_30seeds.Rmd",
    signal_profile = "sparse_correlation",
    within_rho = correlation
  )
}

generate_gaussian_active_component_data <- function(n = 10000L,
                                                    data_seed = 71L,
                                                    pcb_correlation = 0.70,
                                                    gaussian_sd = 1) {
  set.seed(as.integer(data_seed))
  covariates <- generate_primary_benchmark_covariates(n)

  cor_pcb <- matrix(pcb_correlation, 5, 5)
  diag(cor_pcb) <- 1
  cor_pcb[1, 5] <- cor_pcb[5, 1] <- pcb_correlation - 0.05
  cor_pcb[2, 5] <- cor_pcb[5, 2] <- pcb_correlation - 0.03
  sd_pcb <- diag(c(1.5, 1.8, 1.6, 1.7, 1.4))
  pcb <- MASS::mvrnorm(
    n,
    mu = c(3.5, 4.0, 3.8, 4.2, 3.6),
    Sigma = sd_pcb %*% cor_pcb %*% sd_pcb
  )
  colnames(pcb) <- legacy13_groups()$PCBs

  cor_metal <- matrix(c(
    1.00, 0.45, 0.40, 0.35, 0.30,
    0.45, 1.00, 0.45, 0.40, 0.35,
    0.40, 0.45, 1.00, 0.50, 0.25,
    0.35, 0.40, 0.50, 1.00, 0.30,
    0.30, 0.35, 0.25, 0.30, 1.00
  ), 5, 5, byrow = TRUE)
  sd_metal <- diag(c(0.6, 0.8, 0.5, 0.3, 0.25))
  metal_raw <- MASS::mvrnorm(
    n,
    mu = log(c(0.5, 2, 1, 100, 20)),
    Sigma = sd_metal %*% cor_metal %*% sd_metal
  )
  metal <- exp(metal_raw)
  colnames(metal) <- legacy13_groups()$Metals

  cor_pfas <- matrix(c(
    1.00, 0.30, 0.25,
    0.30, 1.00, 0.35,
    0.25, 0.35, 1.00
  ), 3, 3, byrow = TRUE)
  pfas <- pmax(
    MASS::mvrnorm(n, mu = c(3, 1.5, 5), Sigma = cor_pfas * 2),
    0.1
  )
  colnames(pfas) <- legacy13_groups()$PFASs

  x <- data.frame(pcb, metal, pfas, check.names = FALSE)
  finalize_gaussian_primary_benchmark(
    x = x,
    covariates = covariates,
    beta = gaussian_active_component_beta(),
    gaussian_sd = gaussian_sd,
    design_source = "compare_methods_Active_component_identification_30seeds.Rmd",
    signal_profile = "single_active_component",
    within_rho = pcb_correlation
  )
}
