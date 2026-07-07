# Validation Split-Stability Production Bundle

Generated: 2026-07-01

## Contents

- `html/validation_split_stability_standalone_production.html`: rendered R Markdown report.
- `tables/split_stability_summary_production.csv`: direction-level split-stability summary.
- `tables/split_stability_direction_results_production.csv`: all direction-level results for 400 split jobs.
- `tables/split_stability_diagnostics_production.csv`: fit diagnostics and runtime per split job.
- `tables/split_stability_completion_production.csv`: completion table from the Rmd.
- `raw_manifests/split_stability/`: raw per-split RDS manifests.
- `source/`: Rmd, production run script, progress script, and conditions note.
- `logs/`: production run log.
- `progress/`: final progress table.

## Production Conditions

- Four fixed-data scenarios.
- 100 train/validation splits per scenario.
- 200 bootstrap resamples per split.
- `train_prop = 0.60`, `n_quantiles = 4`, `nfolds = 10`, `nlambda = 100`.
- SGL-WQS source: `source/sglwqs_old_analysis/sglwqs-main`.
- Workers: 15.

## Key Results

- All 400 split jobs completed successfully.
- Bootstrap success rate was 1.0 in all scenarios.
- No lambda-boundary or all-zero-exposure failures were observed.
- Global-null false-positive rates were 25/600 (4.2%) for n=500 and 17/600 (2.8%) for n=5000.
- Partial-null active-direction rejection rates were 8/300 (2.7%) for n=500 and 118/300 (39.3%) for n=5000.
- Partial-null true-null rejection rates were 9/300 (3.0%) for n=500 and 6/300 (2.0%) for n=5000.

Interpretation: the production split-stability run supports computational stability of the implementation. The main limitation is not fitting failure, but reduced validation-stage information under low-event small-sample binary settings.
