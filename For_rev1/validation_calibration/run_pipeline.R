#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list(profile = "smoke", workers = 1L, targets = character(0))
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--profile")) {
      i <- i + 1L
      out$profile <- args[[i]]
    } else if (identical(arg, "--workers")) {
      i <- i + 1L
      out$workers <- as.integer(args[[i]])
    } else if (identical(arg, "--targets")) {
      i <- i + 1L
      out$targets <- strsplit(args[[i]], ",", fixed = TRUE)[[1]]
    } else if (arg %in% c("-h", "--help")) {
      cat("Usage: Rscript analysis/validation_calibration/run_pipeline.R --profile smoke --workers 8 [--targets summary_tables]\n")
      quit(status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
if (is.na(file_arg) || !nzchar(file_arg)) {
  file_arg <- "analysis/validation_calibration/run_pipeline.R"
}
root <- normalizePath(dirname(file_arg), mustWork = TRUE)

source(file.path(root, "R", "00_paths.R"))
Sys.setenv(
  VALIDATION_CALIBRATION_ROOT = root,
  VALIDATION_CALIBRATION_PROFILE = args$profile,
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
activate_local_lib()
ensure_vc_dirs()

if (!requireNamespace("targets", quietly = TRUE)) {
  stop(
    "The targets package is required. Install it into this project with: ",
    "install.packages('targets', lib = '", file.path(root, "_r_libs"), "')",
    call. = FALSE
  )
}
if (!requireNamespace("future", quietly = TRUE)) {
  stop("The future package is required for --workers > 1.", call. = FALSE)
}

workers <- max(as.integer(args$workers %||% 1L), 1L)
store <- file.path(root, paste0("_targets_", args$profile))
message("Validation calibration profile: ", args$profile)
message("targets workers: ", workers)
message("targets store: ", store)
message("SGL-WQS source: ", sglwqs_source_path())

target_names <- if (length(args$targets)) args$targets else NULL

if (workers > 1L) {
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  if (is.null(target_names)) {
    targets::tar_make_future(
      workers = workers,
      callr_function = NULL,
      script = file.path(root, "_targets.R"),
      store = store
    )
  } else {
    targets::tar_make_future(
      names = tidyselect::all_of(target_names),
      workers = workers,
      callr_function = NULL,
      script = file.path(root, "_targets.R"),
      store = store
    )
  }
} else {
  if (is.null(target_names)) {
    targets::tar_make(
      script = file.path(root, "_targets.R"),
      store = store
    )
  } else {
    targets::tar_make(
      names = tidyselect::all_of(target_names),
      script = file.path(root, "_targets.R"),
      store = store
    )
  }
}
