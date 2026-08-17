# sglwqs 0.8.13

## Improvements

* Bootstrap aggregation now retains covariate coefficients (`cov_coef`) across
  iterations and stores aggregated summaries in `fit$boot_info$mean_cov_coef`,
  `se_cov_coef`, and percentile intervals.
* Bootstrap aggregation now also stores overall and group-level WQS index-sum
  summaries, making bootstrap-only two-stage summaries available through
  `summary_inference()` and `plot_inference_results()`.
* Checkpointed bootstrap runs now persist and restore covariate coefficient
  matrices alongside exposure coefficient matrices.
* MI pooling helpers prefer bootstrap-aggregated covariate coefficients when
  available, while preserving the existing top-level `fit$cov_coef` contract
  for backward compatibility.
* `sglwqs_mice()` now Rubin-pools downstream GLM coefficients from
  `refit = "full"` as well as `validation = TRUE`, exposing pooled WQS-index
  inference via `fit$pooled$inference` and pooled covariate coefficients via
  `fit$pooled$covariates`.
* MI fits with `bootstrap = TRUE` and no downstream GLM now expose pooled
  bootstrap-derived inference via `fit$pooled$bootstrap_inference`.
* Added `compute_diagnostics()` plus unified `summary_inference()` and
  `plot_inference_results()` helpers so validation, refit, and bootstrap-only
  paths can be accessed through the same entry points.
* Diagnostic messages and bootstrap failure warnings now use fact-only wording
  consistent with the diagnostics API.
