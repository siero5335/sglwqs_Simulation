# Validation Split Stability: Production Conditions

This note describes the split-stability jobs that can be run independently on
another machine.

## Purpose

For fixed generated data sets, repeatedly change only the train/validation
split and summarize validation-stage p-value stability:

- p-value mean, SD, CV, IQR, min, and max
- estimable split count
- sign stability
- validation-stage rejection rate
- retained direction rate

## SGL-WQS version

Use the same SGL-WQS analysis source used for the validation-calibration run.
Set `SGLWQS_SOURCE` to the local path of that source tree before rendering.

Example:

```sh
export SGLWQS_SOURCE=/path/to/sglwqs-main-or-analysis-sglwqs
```

Do not point this to the private development tree unless intentionally testing
that implementation.

## Split-Stability Scenarios

Production uses 4 fixed-data scenarios, each with 100 random train/validation
splits. Each split fit uses 200 bootstrap resamples.

| scenario_id | effect | outcome | n | within_r | between_r | split_reps | n_boot |
|---|---|---|---:|---:|---:|---:|---:|
| `split_global_binary_n500_w0p2_b0` | global null | binomial | 500 | 0.20 | 0.00 | 100 | 200 |
| `split_global_binary_n5000_w0p2_b0` | global null | binomial | 5000 | 0.20 | 0.00 | 100 | 200 |
| `split_partial_binary_n500_w0p6_b0p3` | partial null | binomial | 500 | 0.60 | 0.30 | 100 | 200 |
| `split_partial_binary_n5000_w0p6_b0p3` | partial null | binomial | 5000 | 0.60 | 0.30 | 100 | 200 |

Other production settings:

- `train_prop = 0.60`
- `n_quantiles = 4`
- `nfolds = 10`
- `nlambda = 100`
- `lambda = "lambda.min"`
- `bootstrap = TRUE`
- `parallel_inside_fit = FALSE`
- `minor_threshold = 0.10`
- `asparse = 0.05`
- `alpha = 0.05`
- `checkpoint_interval = 25`
- `cleanup_checkpoint = FALSE`

## Files To Copy

Copy these from the repo to the other machine:

- `analysis/validation_split_stability_standalone.Rmd`
- `analysis/validation_calibration/`

The standalone Rmd sources the existing DGM and fitting functions from
`analysis/validation_calibration/R/`.

## Run Command

From the repo root:

```sh
export SGLWQS_SOURCE=/path/to/sglwqs-main-or-analysis-sglwqs
export SGLWQS_SPLIT_WORKERS=8

Rscript -e 'rmarkdown::render("analysis/validation_split_stability_standalone.Rmd", params = list(profile = "production", workers = 8, run_analysis = TRUE), output_file = "validation_split_stability_standalone.html", output_dir = "analysis", quiet = FALSE)'
```

The Rmd is resumable. It writes one manifest per split job to:

```text
analysis/validation_calibration/output/manifests/production/split_stability/
```

If the render stops, run the same command again. Existing manifest files are
skipped unless `force = TRUE` is passed.

## Output Files

The Rmd writes:

- `analysis/validation_split_stability_standalone.html`
- `analysis/validation_calibration/output/tables/split_stability_direction_results_production.csv`
- `analysis/validation_calibration/output/tables/split_stability_diagnostics_production.csv`
- `analysis/validation_calibration/output/tables/split_stability_summary_production.csv`
- `analysis/validation_calibration/output/tables/split_stability_completion_production.csv`

