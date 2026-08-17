# sglwqs 0.8.13.9001

## Bug fixes

* Bootstrap-only inference summaries now label WQS positive/negative rows as
  loading magnitudes and suppress p-values for these loading rows, avoiding
  confusion with signed downstream refit coefficients.
* Covariate bootstrap summaries and MI bootstrap inference summaries now respect
  the requested confidence level when intervals are displayed.
* Stratified binomial bootstrap resampling now handles singleton strata without
  falling into R's scalar `sample()` behavior.
* MI bootstrap pooling now preserves matrix dimensions when there is only one
  exposure variable.
* `tidy(what = "coefficients", conf.int = TRUE)` now merges bootstrap
  confidence intervals for WQS coefficient rows when intercept or covariate rows
  are present.
* `sglwqs()` and `sglwqs_mice()` now validate that `train_prop` is a single
  finite value strictly between 0 and 1.
* Parallel bootstrap batches now fall back to sequential retry when the future
  backend fails at the batch level, and representative bootstrap errors are
  retained for diagnostics.
* Bootstrap diagnostics now retain per-iteration error messages and
  batch-level future backend failures in returned objects; failed parallel
  backends degrade remaining batches to sequential execution. Multiple
  imputation fits also retain imputation-level error messages.
* Added `vcov.sglwqs()` for downstream refit and bootstrap-only covariance
  extraction.
* Complex-survey support for Gaussian and binomial refits now
  uses `survey_design` sampling weights for quantile cutpoints, mean-normalized
  design weights for sparse-group selection, and the original survey design for
  downstream `survey::svyglm()` inference. Survey mode validates row alignment
  via `analysis_id` or stable row names and supports full-data and
  validation-stage refits. The downstream survey refit is supported, while
  survey-weighted sparse-group selection remains limited by the current
  `sparsegl` weighted backend; fits now retain selection diagnostics and warn
  when weighted selection returns all-zero exposure coefficients or selects a
  lambda value at the edge of the backend path. These all-zero and boundary
  indicators are backend diagnostics, not fit failures. An experimental
  `lambda_path` argument allows sensitivity checks with an explicit backend
  lambda sequence while leaving `lambda` as the coefficient extraction point.
* Survey-aware bootstrap now supports `boot_method = "svrep"` for replicate
  weight resampling. With `survey_design`, `boot_method = "auto"` chooses
  survey replicate weights, normalizes each replicate column for sparse-group
  selection, and computes bootstrap standard errors with the survey design's
  `scale` and `rscales`. `svrep_type = "auto"` uses bootstrap replicate
  weights; jackknife replicates remain available by explicitly setting
  `svrep_type = "JK1"` or `"JKn"`. `boot_method = "naive"` remains available
  as an explicit row-resampling fallback and warns when it ignores survey
  cluster/strata structure.
* Revision follow-up: survey bootstrap summaries now use the stored
  full-sample replicate center and survey-scale standard errors even when
  bootstrap matrices are retained, avoiding ordinary percentile summaries for
  survey replicate matrices.
* Revision follow-up: validation-demo wording now avoids the phrase "valid
  inference", and survey-mode alignment checks now reject empty analysis IDs,
  compare stable analysis row names against `analysis_id`, and carry stable
  row-name IDs into downstream `svyglm` refits.

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
