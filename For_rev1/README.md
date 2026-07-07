# SGL-WQS Revision 1 Analysis Files

This folder contains the main executable analysis files and compact summary
tables used for the Rev1 reviewer response.

## Contents

- `original_simulations/`
  - Rmd/R files for the original simulation analyses and regenerated
    supplementary figure summaries.
- `validation_calibration/`
  - Targets-based validation-stage calibration pipeline.
  - Includes full and non-split target definitions, runner scripts, report Rmds,
    DGP notes, and compact production summaries.
- `validation_split_stability/`
  - Standalone split-stability Rmd and runner script.
  - Includes production split-stability summary tables, not raw split manifests.
- `n100_rare_event_sensitivity/`
  - Targeted n = 100 rare-event sensitivity Rmd.
- `additional_sensitivity/`
  - asparse high-correlation sensitivity.
  - minor-direction threshold sensitivity.
  - effect-estimate bias/coverage.
  - convergence/failure summary.
- `nhanes_voc_demo/`
  - NHANES VOC demonstration Rmd, DAG/variable/group summary, and compact output
    CSVs.
- `highdimensional_simulations/`
  - High-dimensional simulation Rmd/R/sh scripts and compact summary tables.

## Intended Use

Install `sglwqs` and the required R packages, then run the relevant Rmd or
runner script from the corresponding subfolder. Most scripts default to using an
installed `sglwqs` package. To test a local source tree, set the appropriate
environment variable documented in the script, such as `SGLWQS_SOURCE` or
`SGLWQS_SOURCE_PATH`.

## Notes

The files here are meant to be GitHub-friendly. Larger archival result bundles
and manuscript submission files are kept separately.
