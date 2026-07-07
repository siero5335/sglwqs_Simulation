#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p logs

export GT_N_BOOT=100
export GT_N_QUANTILES=4
export GT_N_SEEDS=10
export GT_N_WORKERS=15
export GT_FORCE_RERUN=false

echo "[launcher] started at $(date)"
echo "[launcher] checkpointed group/truth sensitivity run"
echo "[launcher] workers=${GT_N_WORKERS}, boot=${GT_N_BOOT}, seeds=${GT_N_SEEDS}"

Rscript -e 'rmarkdown::render("compare_mixture_methods_highdim_group_truth_checkpointed.Rmd", output_file = "compare_mixture_methods_highdim_group_truth_checkpointed.html")'

echo "[launcher] finished at $(date) status=$?"
