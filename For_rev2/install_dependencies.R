required <- c(
  "MASS", "Matrix", "qgcomp", "gWQS", "groupWQS", "pkgload",
  "dplyr", "tidyr", "purrr", "readr", "ggplot2", "knitr",
  "future.apply", "digest", "remotes"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing)
}

message("Required packages are available.")
message("Exact SGL-WQS source snapshots are bundled under source/ and loaded with pkgload.")

