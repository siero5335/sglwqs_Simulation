find_review_root <- function(start = getwd()) {
  candidates <- unique(c(
    normalizePath(start, mustWork = FALSE),
    normalizePath(file.path(start, ".."), mustWork = FALSE),
    normalizePath(file.path(start, "../.."), mustWork = FALSE)
  ))
  for (path in candidates) {
    if (dir.exists(file.path(path, "source", "sglwqs-main")) ||
        dir.exists(file.path(path, "final_res"))) {
      return(path)
    }
  }
  stop("Could not locate the SGLWQS_rev project root.", call. = FALSE)
}

review_run_mode <- function(default = "smoke") {
  mode <- Sys.getenv("SGLWQS_RUN_MODE", Sys.getenv("RUN_MODE", default))
  mode <- tolower(trimws(mode))
  if (!mode %in% c("smoke", "full")) {
    warning("Unknown run mode '", mode, "'; using 'smoke'.", call. = FALSE)
    mode <- "smoke"
  }
  mode
}

review_workers <- function(run_mode = review_run_mode()) {
  default <- if (identical(run_mode, "full")) 8L else 1L
  workers <- suppressWarnings(as.integer(Sys.getenv("SGLWQS_N_WORKERS", default)))
  if (!is.finite(workers) || workers < 1L) {
    workers <- default
  }
  min(workers, parallel::detectCores(logical = TRUE))
}

review_lapply <- function(x, fun, workers = 1L, seed = TRUE) {
  if (workers > 1L && length(x) > 1L &&
      requireNamespace("future", quietly = TRUE) &&
      requireNamespace("future.apply", quietly = TRUE)) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = workers)
    return(future.apply::future_lapply(x, fun, future.seed = seed))
  }
  lapply(x, fun)
}

review_seeds <- function(n = 30L) {
  seeds <- c(
    71, 42, 123, 256, 314, 500, 617, 789, 888, 999,
    1001, 1123, 1234, 1357, 1500, 1618, 1729, 1847, 1963, 2048,
    2222, 2345, 2500, 2718, 2801, 3001, 3141, 3333, 3500, 3777
  )
  seeds[seq_len(min(n, length(seeds)))]
}

load_review_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required package(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(lapply(pkgs, function(pkg) {
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE)
    )
  }))
}

load_review_sglwqs <- function(source = Sys.getenv("SGLWQS_SOURCE", "installed"),
                               root = find_review_root()) {
  source <- tolower(source)
  if (source %in% c("", "main", "sglwqs-main")) {
    path <- Sys.getenv(
      "SGLWQS_SOURCE_PATH",
      file.path(root, "source", "sglwqs-main")
    )
  } else if (source == "installed") {
    suppressPackageStartupMessages(library(sglwqs))
    return(invisible(list(source = "installed", path = find.package("sglwqs"))))
  } else {
    stop("Unknown SGLWQS_SOURCE: ", source, call. = FALSE)
  }

  if (!dir.exists(path)) {
    stop("Requested sglwqs source directory is missing: ", path, call. = FALSE)
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Package 'pkgload' is required to load local sglwqs source.",
         call. = FALSE)
  }
  pkgload::load_all(path, quiet = TRUE, export_all = FALSE,
                    helpers = FALSE, attach_testthat = FALSE)
  invisible(list(source = basename(path), path = normalizePath(path)))
}

review_groups <- function() {
  list(
    PCBs = c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_187"),
    Metals = c("cd", "pb", "hg", "se", "mg"),
    PFASs = c("PFOA", "PFNA", "PFOS")
  )
}

review_exposures <- function(groups = review_groups()) {
  unname(unlist(groups, use.names = FALSE))
}

near_equicor <- function(p, r) {
  mat <- matrix(r, p, p)
  diag(mat) <- 1
  as.matrix(Matrix::nearPD(mat, corr = TRUE)$mat)
}

simulate_review_mixture <- function(n = 1000, seed = 1, correlation = 0.5,
                                    outcome = c("binomial", "gaussian"),
                                    scenario = c("active", "all_null"),
                                    intercept = -1.4,
                                    noise_sd = 1.0) {
  outcome <- match.arg(outcome)
  scenario <- match.arg(scenario)
  set.seed(seed)

  sex <- rbinom(n, 1, 0.5)
  age <- pmin(pmax(round(rnorm(n, 50, 12)), 20), 80)
  bmi <- pmin(pmax(round(22 + 0.05 * age + 1.5 * sex + rnorm(n, 0, 4), 1), 15), 45)

  cor_pcb <- near_equicor(5, correlation)
  pcb <- MASS::mvrnorm(
    n,
    mu = c(3.5, 4.0, 3.8, 4.2, 3.6),
    Sigma = diag(c(1.5, 1.8, 1.6, 1.7, 1.4)) %*% cor_pcb %*%
      diag(c(1.5, 1.8, 1.6, 1.7, 1.4))
  )
  colnames(pcb) <- review_groups()$PCBs

  cor_metal <- near_equicor(5, correlation)
  metal_raw <- MASS::mvrnorm(
    n,
    mu = log(c(0.5, 2, 1, 100, 20)),
    Sigma = diag(c(0.6, 0.8, 0.5, 0.3, 0.25)) %*% cor_metal %*%
      diag(c(0.6, 0.8, 0.5, 0.3, 0.25))
  )
  metal <- exp(metal_raw)
  colnames(metal) <- review_groups()$Metals

  cor_pfas <- near_equicor(3, min(correlation, 0.4))
  pfas <- pmax(MASS::mvrnorm(n, mu = c(3, 1.5, 5), Sigma = 2 * cor_pfas), 0.1)
  colnames(pfas) <- review_groups()$PFASs

  beta_pcb <- c(pcb_118 = 0, pcb_138 = 0.10, pcb_153 = 0.20,
                pcb_180 = 0, pcb_187 = 0)
  beta_metal <- c(cd = -0.08, pb = -0.18, hg = 0, se = 0, mg = 0)
  beta_pfas <- c(PFOA = 0, PFNA = 0, PFOS = 0)
  if (identical(scenario, "all_null")) {
    beta_pcb[] <- 0
    beta_metal[] <- 0
  }

  eta <- as.numeric(
    scale(pcb) %*% beta_pcb +
      scale(metal) %*% beta_metal +
      scale(pfas) %*% beta_pfas +
      0.25 * sex + 0.02 * ((age - 50) / 10) + 0.04 * ((bmi - 25) / 5)
  )

  if (identical(outcome, "binomial")) {
    prob <- stats::plogis(intercept + eta)
    y <- rbinom(n, 1, prob)
  } else {
    y <- eta + 0.25 * sex + 0.02 * ((age - 50) / 10) +
      0.04 * ((bmi - 25) / 5) + rnorm(n, 0, noise_sd)
  }

  out <- data.frame(pcb, metal, pfas, sex = sex, age = age, bmi = bmi, Y = y)
  attr(out, "true_beta") <- c(beta_pcb, beta_metal, beta_pfas)
  attr(out, "true_group") <- c(
    PCBs = sum(beta_pcb),
    Metals = sum(beta_metal),
    PFASs = sum(beta_pfas)
  )
  out
}

simulate_minor_ratio_data <- function(n = 800, seed = 1, ratio = 0.10,
                                      correlation = 0.5, noise_sd = 1.0) {
  set.seed(seed)
  p <- 8L
  sigma <- near_equicor(p, correlation)
  x <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma)
  colnames(x) <- paste0("mix_", seq_len(p))

  z1 <- rnorm(n)
  z2 <- rbinom(n, 1, 0.5)
  pos_beta <- c(0.18, 0.14, 0.10, 0.08)
  neg_beta <- -ratio * sum(pos_beta) * c(0.40, 0.30, 0.20, 0.10)
  beta <- c(pos_beta, neg_beta)
  eta <- as.numeric(scale(x) %*% beta) + 0.20 * z1 - 0.15 * z2
  y <- eta + rnorm(n, 0, noise_sd)
  out <- data.frame(x, z1 = z1, z2 = z2, Y = y)
  attr(out, "true_beta") <- stats::setNames(beta, colnames(x))
  attr(out, "true_minor_major_ratio") <- ratio
  out
}

simulate_quantile_estimand_data <- function(n = 800, seed = 1, correlation = 0.5,
                                            noise_sd = 2.5) {
  base <- simulate_review_mixture(
    n = n,
    seed = seed,
    correlation = correlation,
    outcome = "gaussian",
    scenario = "all_null",
    noise_sd = 1
  )
  groups <- review_groups()
  exposures <- review_exposures(groups)
  qx <- sapply(base[, exposures, drop = FALSE], function(x) {
    as.numeric(cut(
      x,
      breaks = unique(stats::quantile(x, probs = seq(0, 1, length.out = 5),
                                      na.rm = TRUE, type = 7)),
      include.lowest = TRUE,
      labels = FALSE
    )) - 1
  })
  qx <- as.matrix(qx)
  colnames(qx) <- exposures
  beta <- c(
    pcb_118 = 0.04, pcb_138 = 0.10, pcb_153 = 0.16, pcb_180 = 0.04, pcb_187 = 0.02,
    cd = -0.08, pb = -0.14, hg = 0.06, se = 0.08, mg = 0.00,
    PFOA = 0, PFNA = 0, PFOS = 0
  )
  eta <- as.numeric(qx %*% beta) +
    0.25 * base$sex + 0.03 * ((base$age - 50) / 10) +
    0.08 * ((base$bmi - 25) / 5)
  base$Y <- eta + rnorm(n, 0, noise_sd)
  attr(base, "true_beta_quantile") <- beta
  attr(base, "true_all_exposure_contrast") <- sum(beta)
  attr(base, "true_positive_contrast") <- sum(beta[beta > 0])
  attr(base, "true_negative_contrast") <- sum(beta[beta < 0])
  base
}

safe_sglwqs <- function(data, exposures, groups = NULL, family = "binomial",
                        covariates = c("sex", "age", "bmi"),
                        validation = TRUE, train_prop = 0.6,
                        bootstrap = FALSE, n_boot = 50, parallel = FALSE,
                        seed = 1, nfolds = 3, nlambda = 12,
                        minor_threshold = 0.10, asparse = NULL,
                        verbose = FALSE) {
  warnings <- character(0)
  fit <- withCallingHandlers(
    tryCatch({
      args <- list(
        X = data[, exposures, drop = FALSE],
        y = data$Y,
        covariates = if (length(covariates) > 0) data[, covariates, drop = FALSE] else NULL,
        groups = groups,
        family = family,
        n_quantiles = 4,
        validation = validation,
        train_prop = train_prop,
        bootstrap = bootstrap,
        n_boot = n_boot,
        parallel = parallel,
        seed = seed,
        nfolds = nfolds,
        nlambda = nlambda,
        minor_threshold = minor_threshold,
        verbose = verbose
      )
      if (!is.null(asparse)) {
        args$asparse <- asparse
      }
      do.call(sglwqs::sglwqs, args)
    }, error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {
    list(ok = FALSE, fit = NULL, warnings = warnings, error = conditionMessage(fit))
  } else {
    list(ok = TRUE, fit = fit, warnings = warnings, error = NA_character_)
  }
}

extract_validation_df <- function(fit, conf_level = 0.95) {
  out <- tryCatch(as.data.frame(sglwqs::summary_validation(fit, conf_level = conf_level)),
                  error = function(e) data.frame())
  if (nrow(out) == 0L) {
    return(out)
  }
  if (!"p_value" %in% names(out) && "p.value" %in% names(out)) {
    out$p_value <- out$p.value
  }
  if (!"std_error" %in% names(out) && "std.error" %in% names(out)) {
    out$std_error <- out$std.error
  }
  out$n_train <- attr(out, "n_train") %||% NA_integer_
  out$n_val <- attr(out, "n_val") %||% NA_integer_
  out
}

summarize_split_pvalues <- function(split_df) {
  if (nrow(split_df) == 0L) {
    return(data.frame())
  }
  split_df %>%
    dplyr::group_by(group, direction) %>%
    dplyr::summarise(
      n_estimable = sum(is.finite(p_value)),
      estimate_median = stats::median(estimate, na.rm = TRUE),
      estimate_sign_positive = mean(estimate > 0, na.rm = TRUE),
      estimate_sign_negative = mean(estimate < 0, na.rm = TRUE),
      sign_stability = max(mean(estimate > 0, na.rm = TRUE),
                           mean(estimate < 0, na.rm = TRUE), na.rm = TRUE),
      p_median = stats::median(p_value, na.rm = TRUE),
      p_cv = stats::sd(p_value, na.rm = TRUE) / mean(p_value, na.rm = TRUE),
      p_iqr = stats::IQR(p_value, na.rm = TRUE),
      p_min = suppressWarnings(min(p_value, na.rm = TRUE)),
      p_max = suppressWarnings(max(p_value, na.rm = TRUE)),
      p_range = p_max - p_min,
      sig_rate_0_05 = mean(p_value < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
}

extract_sglwqs_diagnostics_row <- function(fit) {
  diag <- tryCatch(sglwqs::compute_diagnostics(fit), error = function(e) list())
  diag_value <- function(name, default = NA_real_) {
    item <- diag[[name]]
    if (is.null(item) || is.null(item$value)) default else item$value
  }
  lambda_path <- fit$fit$lambda %||% numeric(0)
  selected_lambda <- fit$fit$lambda.min %||% NA_real_
  if (length(lambda_path) > 0L && is.finite(selected_lambda)) {
    lambda_idx <- which.min(abs(lambda_path - selected_lambda))
    lambda_at_edge <- lambda_idx %in% c(1L, length(lambda_path))
  } else {
    lambda_idx <- NA_integer_
    lambda_at_edge <- NA
  }
  data.frame(
    backend_jerr = NA_character_,
    bootstrap_success_rate = diag_value("bootstrap_success_rate"),
    all_zero_weights = !is.null(diag$weights_all_zero),
    validation_estimable = !is.null(fit$validation_info),
    rank_deficient = !is.null(diag$rank_deficient_coefficients),
    selected_lambda = selected_lambda,
    lambda_index = lambda_idx,
    lambda_at_path_edge = lambda_at_edge,
    n_nonzero = diag_value("weights_nonzero_count"),
    stringsAsFactors = FALSE
  )
}

classification_from_error <- function(error) {
  if (is.na(error) || !nzchar(error)) {
    return(NA_character_)
  }
  dplyr::case_when(
    grepl("rank|singular|collinear|not defined", error, ignore.case = TRUE) ~ "rank_deficiency",
    grepl("conver|iterate|jerr|failed to converge", error, ignore.case = TRUE) ~ "backend_convergence",
    grepl("bootstrap", error, ignore.case = TRUE) ~ "bootstrap_failure",
    grepl("validation|glm|contrasts|fitted probabilities|separation", error, ignore.case = TRUE) ~ "validation_refit",
    TRUE ~ "other_error"
  )
}

lazy_load_cache <- function(base) {
  if (!file.exists(paste0(base, ".rdx")) || !file.exists(paste0(base, ".rdb"))) {
    return(NULL)
  }
  env <- new.env(parent = emptyenv())
  lazyLoad(base, envir = env)
  as.list(env)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
