`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

r2r2_repo_root <- function() {
  normalizePath(file.path(Sys.getenv("R2R2_REPO_ROOT", unset = getwd())),
    mustWork = TRUE
  )
}

r2r2_root <- function() {
  normalizePath(file.path(Sys.getenv("R2R2_ROOT", unset = file.path(r2r2_repo_root(), "reviewer2_round2"))),
    mustWork = TRUE
  )
}

r2r2_file <- function(...) {
  file.path(r2r2_root(), ...)
}

r2r2_result_file <- function(...) {
  r2r2_file("results", ...)
}

r2r2_make_dirs <- function() {
  dirs <- c(
    "config", "R", "scripts",
    file.path("results", "raw"),
    file.path("results", "summaries"),
    file.path("results", "tables"),
    file.path("results", "figures"),
    file.path("results", "logs"),
    file.path("results", "quarantine")
  )
  invisible(lapply(file.path(r2r2_root(), dirs), dir.create,
    showWarnings = FALSE, recursive = TRUE
  ))
}

r2r2_required_packages <- function() {
  c(
    "MASS", "Matrix", "qgcomp", "gWQS", "groupWQS", "pkgload",
    "dplyr", "tidyr", "purrr", "readr", "ggplot2", "knitr", "future.apply",
    "digest"
  )
}

r2r2_load_packages <- function() {
  missing <- r2r2_required_packages()[
    !vapply(r2r2_required_packages(), requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(MASS)
    library(Matrix)
    library(qgcomp)
    library(gWQS)
    library(groupWQS)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(readr)
    library(ggplot2)
  })
}

r2r2_sglwqs_source_path <- function() {
  path <- Sys.getenv(
    "R2R2_SGLWQS_SOURCE",
    unset = file.path(r2r2_repo_root(), "source", "sglwqs_old_analysis", "sglwqs-main")
  )
  if (!dir.exists(path)) {
    fallback <- file.path(r2r2_repo_root(), "submitted_review_zip", "sglwqs-main")
    if (dir.exists(fallback)) {
      path <- fallback
    }
  }
  normalizePath(path, mustWork = TRUE)
}

r2r2_load_local_sglwqs <- function() {
  pkgload::load_all(r2r2_sglwqs_source_path(), quiet = TRUE, export_all = FALSE)
  invisible(TRUE)
}

r2r2_settings <- function(n_boot = NULL, smoke = FALSE) {
  list(
    quantiles = as.integer(Sys.getenv("R2R2_QUANTILES", "4")),
    train_prop = as.numeric(Sys.getenv("R2R2_TRAIN_PROP", "0.6")),
    n_boot = as.integer(n_boot %||% Sys.getenv("R2R2_N_BOOT", if (smoke) "10" else "200")),
    nfolds = as.integer(Sys.getenv("R2R2_NFOLDS", "10")),
    lambda = Sys.getenv("R2R2_LAMBDA", "lambda.min"),
    asparse = as.numeric(Sys.getenv("R2R2_ASPARSE", "0.05")),
    minor_threshold = as.numeric(Sys.getenv("R2R2_MINOR_THRESHOLD", "0.10")),
    alpha = as.numeric(Sys.getenv("R2R2_ALPHA", "0.05")),
    save_raw_fit = Sys.getenv("R2R2_SAVE_RAW_FIT", "false") %in% c("1", "true", "TRUE"),
    sglwqs_source = r2r2_sglwqs_source_path()
  )
}

r2r2_seeds <- function(n = 30L) {
  seeds <- c(
    71, 42, 123, 256, 314, 500, 617, 789, 888, 999,
    1001, 1123, 1234, 1357, 1500, 1618, 1729, 1847, 1963, 2048,
    2222, 2345, 2500, 2718, 2801, 3001, 3141, 3333, 3500, 3777
  )
  seeds[seq_len(min(as.integer(n), length(seeds)))]
}

r2r2_methods <- function() {
  c("SGL-WQS", "gWQS", "groupWQS", "qgcomp")
}

r2r2_set_thread_env <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
}

safe_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("_+", "_", x)
}

job_output_path <- function(job) {
  scenario <- safe_filename(job$scenario_id)
  method <- safe_filename(job$method)
  data_seed <- sprintf("seed%s", as.integer(job$data_seed))
  split <- if ("split_replicate" %in% names(job) && !is.na(job$split_replicate)) {
    sprintf("_split%04d", as.integer(job$split_replicate))
  } else {
    ""
  }
  dir <- r2r2_result_file("raw", safe_filename(job$battery), scenario, method)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  file.path(dir, paste0(data_seed, split, ".rds"))
}

atomic_save_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) {
    file.copy(tmp, path, overwrite = TRUE)
    unlink(tmp)
  }
  invisible(path)
}

is_complete_job_file <- function(path) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  ok <- tryCatch({
    x <- readRDS(path)
    is.list(x) && isTRUE(x$schema_version == "reviewer2_round2_v1") &&
      !is.null(x$method_metrics) && nrow(x$method_metrics) == 1L
  }, error = function(e) FALSE)
  if (!ok) {
    quarantine_dir <- r2r2_result_file("quarantine")
    dir.create(quarantine_dir, showWarnings = FALSE, recursive = TRUE)
    file.rename(path, file.path(quarantine_dir, paste0(basename(path), ".bad.", format(Sys.time(), "%Y%m%d%H%M%S"))))
  }
  ok
}

atomic_read_results <- function(paths) {
  paths <- paths[file.exists(paths)]
  out <- vector("list", length(paths))
  keep <- logical(length(paths))
  for (i in seq_along(paths)) {
    out[[i]] <- tryCatch(readRDS(paths[[i]]), error = function(e) NULL)
    keep[[i]] <- is.list(out[[i]]) && isTRUE(out[[i]]$schema_version == "reviewer2_round2_v1")
  }
  out[keep]
}

make_train_validation_split <- function(n, split_seed, train_prop = 0.6) {
  set.seed(as.integer(split_seed))
  idx <- seq_len(as.integer(n))
  train_n <- round(length(idx) * train_prop)
  train_idx <- sort(sample(idx, size = train_n, replace = FALSE))
  valid_idx <- setdiff(idx, train_idx)
  list(train = train_idx, validation = valid_idx)
}

assert_or_stop <- function(ok, message) {
  if (!isTRUE(ok)) {
    stop(message, call. = FALSE)
  }
}

classify_error <- function(message) {
  message <- paste(message, collapse = " | ")
  if (!nzchar(message)) return(NA_character_)
  if (grepl("lambda|path|boundary", message, ignore.case = TRUE)) return("lambda_path")
  if (grepl("conver|jerr|maxit|irls|separation|fitted probabilities", message, ignore.case = TRUE)) return("backend_convergence")
  if (grepl("rank|singular|alias|contrasts", message, ignore.case = TRUE)) return("rank_deficiency")
  if (grepl("class|levels|0.*1|both", message, ignore.case = TRUE)) return("outcome_class")
  if (grepl("parallel|future|fork|cluster", message, ignore.case = TRUE)) return("parallel_backend")
  "other"
}

finite_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.finite(x), x, NA_real_)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else NA_real_
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

package_version_table <- function() {
  pkgs <- unique(c(r2r2_required_packages(), "sglwqs", "sparsegl", "glmnet"))
  data.frame(
    package = pkgs,
    version = vapply(pkgs, function(pkg) {
      if (requireNamespace(pkg, quietly = TRUE)) {
        as.character(utils::packageVersion(pkg))
      } else {
        NA_character_
      }
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

git_commit_for_path <- function(path) {
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  if (!dir.exists(path)) return(NA_character_)
  setwd(path)
  out <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE), error = function(e) NA_character_)
  if (length(out) == 1L && grepl("^[0-9a-f]{7,40}$", out)) out else NA_character_
}

write_session_artifacts <- function(prefix = "production") {
  r2r2_make_dirs()
  writeLines(capture.output(sessionInfo()),
    con = r2r2_result_file("logs", paste0(prefix, "_sessionInfo.txt"))
  )
  atomic_write_csv(package_version_table(),
    r2r2_result_file("logs", paste0(prefix, "_package_versions.csv"))
  )
  source_meta <- data.frame(
    item = c("repo_root", "reviewer2_round2_root", "sglwqs_source", "sglwqs_source_git_commit", "workspace_git_commit"),
    value = c(
      r2r2_repo_root(),
      r2r2_root(),
      r2r2_sglwqs_source_path(),
      git_commit_for_path(r2r2_sglwqs_source_path()),
      git_commit_for_path(r2r2_repo_root())
    ),
    stringsAsFactors = FALSE
  )
  atomic_write_csv(source_meta,
    r2r2_result_file("logs", paste0(prefix, "_source_metadata.csv"))
  )
  invisible(source_meta)
}

source_r2r2 <- function(file) {
  source(r2r2_file("R", file), chdir = TRUE)
}
