#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

export SGLWQS_SPLIT_WORKERS="${SGLWQS_SPLIT_WORKERS:-8}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

Rscript -e 'rmarkdown::render(
  "validation_split_stability_standalone.Rmd",
  params = list(
    profile = "production",
    workers = as.integer(Sys.getenv("SGLWQS_SPLIT_WORKERS", "8")),
    run_analysis = TRUE,
    force = FALSE
  ),
  output_file = "validation_split_stability_standalone_production.html",
  output_dir = ".",
  quiet = FALSE
)'
