#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"
mkdir -p logs

export HD_N_SAMPLE=1000
export HD_N_BOOT=200
export HD_N_QUANTILES=4
export HD_N_GROUPS=10
export HD_P_VALUES=50,100,200
export HD_N_SEEDS=30
export HD_N_WORKERS=15
export HD_FORCE_RERUN=false

echo "[launcher] started at $(date)"
echo "[launcher] checkpointed high-dimensional p-sensitivity run, p<=200"
echo "[launcher] workers=${HD_N_WORKERS}, n=${HD_N_SAMPLE}, boot=${HD_N_BOOT}, seeds=${HD_N_SEEDS}, p=${HD_P_VALUES}"

Rscript -e 'rmarkdown::render("compare_mixture_methods_highdim_p_checkpointed.Rmd", output_file = "compare_mixture_methods_highdim_p_checkpointed_p200.html")'

echo "[launcher] finished at $(date) status=$?"
