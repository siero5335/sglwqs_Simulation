# sglwqs: Weighted Quantile Sum Regression with Sparse Group Lasso

## Overview

`sglwqs` is an R package for Weighted Quantile Sum (WQS) regression with Sparse Group Lasso. It evaluates mixture effects of multiple exposure variables and simultaneously estimates positive and negative direction weights.

## Features

- **Quantile Transformation**: Automatically converts continuous variables into quantile-based categories
- **Sparse Group Lasso**: Performs variable selection and weight estimation simultaneously
- **Bidirectional Weights**: Separately estimates positive and negative effects
- **Grouping**: Applies penalties by variable groups (e.g., PCBs, Phthalates)
- **Bootstrap**: Bootstrap aggregation for weight stabilization
- **Train/Validation Split**: Data splitting for valid p-values
- **Multiple Imputation (mice)**: Handles missing data with pooling via Rubin's rules
- **Categorical Covariates**: Automatic dummy coding for factor variables
- **Rich Visualization**: Multiple plot types for visualizing weights

## Installation

### Dependencies (install first)

```r
install.packages(c("sparsegl", "ggplot2", "tidyr", "dplyr"))

# Optional (recommended)
install.packages(c("future", "future.apply"))  # Parallel processing
install.packages("patchwork")                   # Combined plots
install.packages("mice")                        # Multiple imputation
```

### Package Installation

```r
# Method 1: Install from tar.gz (recommended)
install.packages("sglwqs_0.8.13.tar.gz", repos = NULL, type = "source")

# Method 2: Install from directory
install.packages("/path/to/sglwqs", repos = NULL, type = "source")

# Method 3: Using devtools
devtools::install("/path/to/sglwqs")

# Method 4: From GitHub (after publication)
# devtools::install_github("username/sglwqs")
```

### Using Without Installing

```r
# Manually source each file (order matters)
source("R/quantile_transform.R")
source("R/utils.R")
source("R/core_functions.R")
source("R/sglwqs_main.R")
source("R/model_methods.R")
source("R/inference_functions.R")
source("R/mice_functions.R")
source("R/plot_weights.R")  # Visualization functions
```

## Usage

### Choosing an Inference Path

`sglwqs` separates weight estimation from downstream effect-size estimation.
The main configuration choice is therefore not just "which model", but also
"which uncertainty path".

| Goal | Suggested settings | What you get |
|---|---|---|
| Quick prototype / point estimates only | `bootstrap = FALSE`, `refit = "none"` | Weights and sparse coefficients only |
| Weight stability without downstream GLM | `bootstrap = TRUE`, `refit = "none"` | Bootstrap summaries for exposures, covariates, and WQS index sums |
| Full-data downstream GLM (small `n`) | `refit = "full"` | In-sample downstream GLM inference |
| Held-out downstream GLM | `validation = TRUE` | Validation-split downstream GLM inference |
| Missing data + downstream GLM | `sglwqs_mice(...)` with `refit = "full"` or `validation = TRUE` | Rubin-pooled downstream inference |
| Missing data + bootstrap-only summaries | `sglwqs_mice(...)` with `bootstrap = TRUE`, `refit = "none"` | Rubin-pooled bootstrap summaries |

Unified accessors:

- `summary_inference(fit)`
- `plot_inference_results(fit)`
- `compute_diagnostics(fit)`

### Basic Usage

```r
library(sglwqs)
data("sglwqs_example")

# Define variable groups
metals <- paste0("metal", 1:4)
pesticides <- paste0("pesticide", 1:4)

# Fit the model
fit <- sglwqs(
  X = sglwqs_example[, c(metals, pesticides)],
  y = sglwqs_example$outcome_cont,
  covariates = sglwqs_example[, c("age", "sex")],
  groups = list(Metals = metals, Pesticides = pesticides),
  n_quantiles = 4
)

# View results
summary(fit)
```

The same can be written using `data` with column name specification. This syntax is convenient when aligning with `sglwqs_mice()`.

```r
fit_named <- sglwqs(
  data = sglwqs_example,
  X = c(metals, pesticides),
  y = "outcome_cont",
  covariates = c("age", "sex"),
  groups = list(Metals = metals, Pesticides = pesticides),
  n_quantiles = 4
)
```

### Standard Methods (Publication-Ready)

```r
# Coefficients and weights (standard interface)
coef(fit, what = "coefficients", direction = "both")
coef(fit, what = "weights", direction = "positive")

# Prediction on new data (reuses quantile breaks from training)
pred_idx <- predict(fit, newdata = sglwqs_example[1:5, ], type = "index")
pred_y <- predict(fit, newdata = sglwqs_example[1:5, ], type = "response")

# CV curve (lambda selection visualization)
plot_cv(fit)

# Confidence intervals (bootstrap / validation)
fit_boot <- sglwqs(
  X = sglwqs_example[, c(metals, pesticides)],
  y = sglwqs_example$outcome_cont,
  bootstrap = TRUE,
  n_boot = 50,
  seed = 123
)
confint(fit_boot, type = "bootstrap", direction = "positive")

# Prediction uncertainty (requires bootstrap or validation)
pred_se <- predict(fit_boot, newdata = sglwqs_example[1:5, ], se.fit = TRUE)
pred_ci <- predict_interval(fit_boot, newdata = sglwqs_example[1:5, ], level = 0.95)
```

### Integration and Evaluation (broom / diagnostics)

```r
# broom-compatible (generics::tidy / glance / augment)
generics::tidy(fit, what = "weights")
generics::glance(fit)
generics::augment(fit, data = sglwqs_example[1:20, c(metals, pesticides)])

# Predictive performance (validation model, or explicit outcome)
fit_val <- sglwqs(
  X = sglwqs_example[, c(metals, pesticides)],
  y = sglwqs_example$outcome_bin,
  family = "binomial",
  validation = TRUE,
  train_prop = 0.6,
  seed = 123
)
validation_metrics(fit_val)

# Calibration
cal <- calibrate(fit_val, n_bins = 10)
plot(cal)
```

### Bootstrap Aggregation (Weight Stabilization)

When sample sizes are small or correlations among variables are high, bootstrap can stabilize weights:

```r
fit_boot <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  bootstrap = TRUE,
  n_boot = 100,
  seed = 123
)

# Get weights including selection frequencies
weights_df <- extract_weights(fit_boot)

# Unified bootstrap-only inference summary
summary_inference(fit_boot)
plot_inference_results(fit_boot)
```

### Train/Validation Split (Inference)

Split data to avoid post-selection inference issues and obtain valid p-values:

```r
fit_val <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  validation = TRUE,
  train_prop = 0.6,  # 60% training, 40% validation
  seed = 123
)

# Check p-values (refer to group_results when groups are specified)
summary_validation(fit_val)
summary_inference(fit_val)

# Example: p-values for PCBs
fit_val$validation_info$group_results$PCBs$pos_pvalue
fit_val$validation_info$group_results$PCBs$neg_pvalue
```

### Bootstrap + Validation (Recommended)

Combine both for the most robust results:

```r
fit_full <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  bootstrap = TRUE,
  n_boot = 100,
  validation = TRUE,
  train_prop = 0.6,
  seed = 123
)

summary(fit_full)
```

When grouped validation or full refit is used, `minor_threshold` can be used
to exclude a group's minor direction from the downstream GLM when its
coefficient-mass ratio is very small. The default `minor_threshold = 0.10`
excludes a direction when `minor / major < 10%`; set `minor_threshold = 0`
to disable this behavior.

### Full Refit and Complex Survey

To explicitly request a full-sample refit, use `refit = "full"`. `obs_weights` are used in the sparse-group selection loss, and `quantile_weights` are used for quantile cutpoints. `refit_engine = "glm"` performs a standard model-based full refit and does not automatically incorporate survey weights. The intended inference route for NHANES in v1 is `refit_engine = "svyglm"`, where the user passes a pre-created `survey_design`.

```r
# Weighted selection + full-sample glm refit
fit_refit <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  obs_weights = wqs_data$sample_weight,
  refit = "full",
  refit_engine = "glm"
)

# survey::svydesign must be created by the user beforehand
fit_svy <- sglwqs(
  X = nhanes_panel[, exposure_vars],
  y = nhanes_panel$hba1c,
  covariates = nhanes_panel[, c("age", "sex", "bmi")],
  obs_weights = nhanes_panel$wt_subsample,
  refit = "full",
  refit_engine = "svyglm",
  survey_design = nhanes_design
)
```

In survey mode with `refit_engine = "svyglm"`, only `refit = "full"` is supported in v1. Using `validation = TRUE`, `bootstrap = TRUE`, or `parallel = TRUE` (which would expect survey-aware resampling) will raise explicit errors. Weighted point summaries are available through `validation_metrics()` and `calibrate()` using analysis weights, but weighted AUC is intentionally returned as `NA` rather than using an unvalidated approximation.

### Categorical Covariates

Factor and character covariates are automatically dummy-coded:

```r
# Covariates including factor variables
wqs_data$race <- factor(sample(c("White", "Black", "Asian"), nrow(wqs_data), replace = TRUE))

fit <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data[, c("sex", "race")],  # race is auto dummy-coded
  groups = list(PCBs = PCBs, Phthalates = phthalates)
)
```

### Group Structure

The current implementation only supports `group_structure = "direction"`.
Each compound group's **positive (beta+) and negative (beta-) directions are treated as separate groups**, allowing mixed-direction effects within the same group.

> Note: Specifying anything other than `"direction"` will produce a warning and fall back to `"direction"`.

```r
fit_direction <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  groups = list(PCBs = PCBs, Phthalates = phthalates)
  # group_structure = "direction"  # default
)
```

### Weight Visualization

```r
# Butterfly chart of weights (positive/negative)
plot_weights(fit)

# Cross-validation results
plot_cv(fit)

# (When bootstrap is used)
# plot_selection_frequency(fit_boot, top_n = 10)
# plot_bootstrap_ci(fit_boot, top_n = 10)
```

### Bootstrap/Validation Output and Visualization

```r
# Bootstrap summary
fit <- sglwqs(..., bootstrap = TRUE, n_boot = 100)
boot_summary <- summary_bootstrap(fit)
print(boot_summary)

# Selection frequency plot
plot_selection_frequency(fit, top_n = 15)

# Bootstrap confidence interval plot
plot_bootstrap_ci(fit)

# Validation summary
fit <- sglwqs(..., validation = TRUE)
val_summary <- summary_validation(fit)
print(val_summary)

# Validation results visualization
plot_validation_results(fit)
plot_inference_results(fit)

# Comprehensive inference table (for publications)
inf_table <- get_inference_table(fit, top_n = 10)
print(inf_table)

# Combined plot (patchwork package recommended)
plot_combined_results(fit)

# Fact-only diagnostics
diag <- compute_diagnostics(fit)
print(diag)
```

### Parallel Processing

```r
library(future)
plan(multisession, workers = 4)  # Use 4 cores

fit <- sglwqs(
  X = exposures,
  y = outcome,
  bootstrap = TRUE,
  n_boot = 100,
  parallel = TRUE,  # Enable parallel processing
  seed = 123
)

plan(sequential)  # Reset after completion
```

### Multiple Imputation (mice)

Use multiply imputed datasets generated by the mice package for data with missing values.
Results are pooled according to Rubin's rules.

```r
library(mice)
library(sglwqs)

# Multiple imputation
imp <- mice(data_with_missing, m = 10, method = "pmm", seed = 123)

# Apply sglwqs to MI data
fit_mi <- sglwqs_mice(
  mids_obj = imp,
  exposure_vars = c("x1", "x2", "x3", "x4"),
  outcome_var = "y",
  covariate_vars = c("age", "sex"),
  groups = list(GroupA = c("x1", "x2"), GroupB = c("x3", "x4")),
  validation = TRUE,
  seed = 123
)

# View results
summary(fit_mi)

# Pooled p-values via Rubin's rules
fit_mi$pooled$validation$wqs_pos$p_value
fit_mi$pooled$validation$wqs_neg$p_value

# Visualize pooled weights
plot_pooled_weights(fit_mi)

# Pooled validation forest plot
plot_validation_results(fit_mi)

# Combined plot also available
plot_combined_results(fit_mi)
```

**Pooling via Rubin's Rules:**
- Estimates: Mean of estimates across imputed datasets
- Variance: Combination of within-imputation and between-imputation variance
- Degrees of freedom: Barnard-Rubin adjustment
- FMI (Fraction of Missing Information): Indicator of information loss due to missingness

## Function Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `sglwqs()` | Fit an SGL-WQS model |
| `quantile_transform()` | Quantile transformation of data |
| `extract_weights()` | Extract weights (including bootstrap info) |
| `extract_weights_grouped()` | Extract weights with group information |
| `summary_by_group()` | Summary by group |

### Inference Functions

| Function | Description |
|----------|-------------|
| `summary_bootstrap()` | Bootstrap summary (selection frequency, CI) |
| `summary_validation()` | Validation summary (p-values, estimates) |
| `get_inference_table()` | Comprehensive inference table for publications |

### Multiple Imputation Functions

| Function | Description |
|----------|-------------|
| `sglwqs_mice()` | Run SGL-WQS with mice objects |
| `extract_pooled_weights()` | Extract pooled weights |
| `plot_pooled_weights()` | Visualize pooled weights |

### Visualization Functions

| Function | Description |
|----------|-------------|
| `plot_weights()` | Butterfly chart of weights (positive/negative) |
| `plot_cv()` | Cross-validation results |
| `plot_selection_frequency()` | Bootstrap selection frequency plot |
| `plot_bootstrap_ci()` | Bootstrap coefficient confidence interval plot |
| `plot_validation_results()` | Validation inference results plot |
| `plot_combined_results()` | Combined plot (requires patchwork) |

## Key Arguments

| Argument | Description | Default |
|----------|-------------|---------|
| `groups` | List of compound groups | NULL |
| `group_structure` | "direction" (only option currently) | "direction" |
| `bootstrap` | Whether to use bootstrap | FALSE |
| `n_boot` | Number of bootstrap iterations | 100 |
| `parallel` | Whether to use parallel processing | FALSE |
| `validation` | Whether to split data | FALSE |
| `train_prop` | Proportion of training data | 0.6 |

## Output Object

The return value of `sglwqs()` includes:

- `fit`: cv.sparsegl fit object
- `pos_weights`, `neg_weights`: Normalized weights
- `pos_coef`, `neg_coef`: Raw coefficients
- `pos_index_sum`, `neg_index_sum`: Index sums
- `quantile_breaks`: Quantile boundaries (applicable to new data)
- `boot_info`: Bootstrap information (when bootstrap=TRUE)
  - `selection_freq_pos`, `selection_freq_neg`: Selection frequencies
  - `se_pos_coef`, `se_neg_coef`: Standard errors of coefficients
- `validation_info`: Validation information (when validation=TRUE)
  - `group_results`: Estimates, SE, and p-values by group (when groups are specified)
  - `wqs_pos_pvalue`, `wqs_neg_pvalue`: Overall WQS p-values (when groups are not specified)
  - `wqs_pos_estimate`, `wqs_neg_estimate`: Overall WQS estimates (when groups are not specified)

## Interpretation of Weights and Indices

### Weight Definition

In SGL-WQS, weights are computed as follows:

1. **Within-group normalization**: Weights are normalized so that the sum within each compound group equals 1
2. **Equal group contribution**: Each group contributes equally to the overall WQS index

```
WQS_pos = (1/G) * Sigma_g [ Sigma_{j in g} w_j^{(g)} * X_j ]

where:
- G = number of groups
- w_j^{(g)} = normalized weight of variable j in group g (within-group sum = 1)
- X_j = quantile-transformed exposure variable
```

### Index Sum Meaning

- `pos_index_sum` / `neg_index_sum`: **Sum of raw SGL coefficients**
- This serves as an indicator of the "effect magnitude" for each direction
- Different from the sum of weights (which always equals 1)

### Group-Level Interpretation

```r
# Index sum by group
fit$pos_index_sum_by_group  # Coefficient sum per group (contribution)

# Within-group weights (normalized to sum to 1)
summary_by_group(fit)
```

**Example:**
```
Group          Index Sum    Interpretation
PCBs           0.35         Contributes 35% of mixture effect
Phthalates     0.15         Contributes 15% of mixture effect
```

## License

MIT License

## References

- Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2013). A sparse-group lasso. *Journal of Computational and Graphical Statistics*.
- Carrico, C., Gennings, C., Wheeler, D. C., & Factor-Litvak, P. (2015). Characterization of weighted quantile sum regression for highly correlated data in a risk analysis setting. *Journal of Agricultural, Biological, and Environmental Statistics*.
- Keil, A. P., et al. (2020). A quantile-based g-computation approach to addressing the effects of exposure mixtures. *Environmental Health Perspectives*.
- Renzetti, S., Gennings, C., & Calza, S. (2023). A weighted quantile sum regression with penalized weights and two indices. *Frontiers in Public Health*.
