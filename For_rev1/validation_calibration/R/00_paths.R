`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

vc_root <- function() {
  root <- Sys.getenv("VALIDATION_CALIBRATION_ROOT", unset = NA_character_)
  if (!is.na(root) && nzchar(root)) {
    return(normalizePath(root, mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

vc_file <- function(...) {
  file.path(vc_root(), ...)
}

repo_root <- function() {
  normalizePath(file.path(vc_root(), "..", ".."), mustWork = TRUE)
}

activate_local_lib <- function() {
  lib <- vc_file("_r_libs")
  if (dir.exists(lib)) {
    .libPaths(unique(c(normalizePath(lib, mustWork = TRUE), .libPaths())))
  }
  invisible(.libPaths())
}

ensure_vc_dirs <- function() {
  dirs <- c(
    "output",
    "output/html",
    "output/tables",
    "output/figures",
    "output/manifests",
    "output/checkpoints",
    "cache"
  )
  for (d in dirs) {
    dir.create(vc_file(d), showWarnings = FALSE, recursive = TRUE)
  }
  invisible(TRUE)
}

sglwqs_source_path <- function() {
  path <- Sys.getenv("SGLWQS_SOURCE", "")
  if (!nzchar(path)) {
    return(NA_character_)
  }
  normalizePath(path, mustWork = TRUE)
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp-", Sys.getpid(), "-", sample.int(1e8, 1))
  saveRDS(object, tmp)
  ok <- file.rename(tmp, path)
  if (!ok) {
    stop("Failed to move temporary RDS into place: ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  tmp <- paste0(path, ".tmp-", Sys.getpid(), "-", sample.int(1e8, 1))
  readr::write_csv(x, tmp, na = "")
  ok <- file.rename(tmp, path)
  if (!ok) {
    stop("Failed to move temporary CSV into place: ", path, call. = FALSE)
  }
  invisible(path)
}

source_vc <- function(file) {
  source(vc_file("R", file), chdir = TRUE)
}
