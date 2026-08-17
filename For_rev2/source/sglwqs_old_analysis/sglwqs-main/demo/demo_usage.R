# =============================================================================
# sglwqs Package Demo
# =============================================================================
# This script demonstrates the features of the sglwqs package:
# - Basic model fitting
# - Chemical group-based Group Lasso
# - Bootstrap aggregation for stable weights
# - Train/Validation split for inference
# - Parallel processing for bootstrap
# - Inference summary functions
# - Visualization functions for bootstrap/validation results
# - Multiple imputation with mice

# =============================================================================
# INSTALLATION
# =============================================================================
# Option 1: Install from tar.gz (recommended)
# install.packages("sglwqs_0.8.10.tar.gz", repos = NULL, type = "source")

# Option 2: Install from directory
# install.packages("/path/to/sglwqs", repos = NULL, type = "source")

# Option 3: For development/testing without installation
# source all R files in order:
# source("R/quantile_transform.R")
# source("R/utils.R")
# source("R/core_functions.R")
# source("R/sglwqs_main.R")
# source("R/inference_functions.R")
# source("R/mice_functions.R")
# source("R/plot_weights.R")

# Load the package (after installation)
library(sglwqs)

# =============================================================================
# Load Example Data
# =============================================================================
if (!require(gWQS, quietly = TRUE)) {
  install.packages("gWQS")
}
library(gWQS)
data("wqs_data")

# Define chemical groups
PCBs <- names(wqs_data)[grepl("^LBX", names(wqs_data))]
phthalates <- names(wqs_data)[grepl("^URX", names(wqs_data))]

cat("PCBs:", length(PCBs), "variables\n")
cat("Phthalates:", length(phthalates), "variables\n")

# =============================================================================
# 1. Basic Model Fitting
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("1. BASIC MODEL FITTING\n")
cat(strrep("=", 60), "\n")

fit_basic <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  group_structure = "direction"
)

print(fit_basic)

# =============================================================================
# 2. Bootstrap Aggregation
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("2. BOOTSTRAP AGGREGATION\n")
cat(strrep("=", 60), "\n")

fit_boot <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  bootstrap = TRUE,
  n_boot = 50,  # Use 100+ for production
  seed = 123,
  verbose = TRUE
)

print(fit_boot)

# --- Bootstrap Summary ---
cat("\n--- Bootstrap Summary ---\n")
boot_summary <- summary_bootstrap(fit_boot, conf_level = 0.95, min_freq = 0.1)
print(boot_summary)

# --- Selection Frequency Plot ---
cat("\n--- Generating Selection Frequency Plot ---\n")
p_sel_freq <- plot_selection_frequency(fit_boot, top_n = 10, min_freq = 0.05)
print(p_sel_freq)

# --- Bootstrap CI Plot ---
cat("\n--- Generating Bootstrap CI Plot ---\n")
p_boot_ci <- plot_bootstrap_ci(fit_boot, top_n = 10)
print(p_boot_ci)

# =============================================================================
# 3. Train/Validation Split for Inference
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("3. TRAIN/VALIDATION SPLIT FOR INFERENCE\n")
cat(strrep("=", 60), "\n")

fit_val <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  validation = TRUE,
  train_prop = 0.6,
  seed = 123,
  verbose = TRUE
)

print(fit_val)

# --- Validation Summary ---
cat("\n--- Validation Summary ---\n")
val_summary <- summary_validation(fit_val, conf_level = 0.95)
print(val_summary)

# --- Validation Results Plot ---
cat("\n--- Generating Validation Results Plot ---\n")
p_val <- plot_validation_results(fit_val, include_covariates = TRUE)
print(p_val)

# =============================================================================
# 4. Bootstrap + Validation Combined (RECOMMENDED)
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("4. BOOTSTRAP + VALIDATION COMBINED (RECOMMENDED)\n")
cat(strrep("=", 60), "\n")

fit_full <- sglwqs(
  X = wqs_data[, c(PCBs, phthalates)],
  y = wqs_data$y,
  covariates = wqs_data["sex"],
  groups = list(PCBs = PCBs, Phthalates = phthalates),
  bootstrap = TRUE,
  n_boot = 50,
  validation = TRUE,
  train_prop = 0.6,
  seed = 123,
  verbose = TRUE
)

# --- GLM-style Summary ---
cat("\n--- GLM-style Summary ---\n")
summary(fit_full)

# The summary now shows:
# - Coefficient table like GLM output
# - Group-specific effects (PCBs positive/negative, Phthalates positive/negative)
# - Standard errors, z/t values, and p-values
# - Significance codes

# --- Group-wise Validation Results ---
cat("\n--- Group-wise Validation Results ---\n")
val_summary <- summary_validation(fit_full)
print(val_summary)

# This shows inference results separately for each chemical group:
# - PCBs (positive) / PCBs (negative)
# - Phthalates (positive) / Phthalates (negative)
# - Each with estimate, SE, CI, and p-value

# --- Bootstrap Summary ---
cat("\n--- Bootstrap Summary ---\n")
boot_summary <- summary_bootstrap(fit_full, min_freq = 0.5)
print(boot_summary)

# --- Combined Results Plot ---
cat("\n--- Generating Combined Results Plot ---\n")
# Note: Requires 'patchwork' package for combined layout
if (requireNamespace("patchwork", quietly = TRUE)) {
  p_combined <- plot_combined_results(fit_full, top_n = 8)
  print(p_combined)
} else {
  cat("Install 'patchwork' package for combined plot: install.packages('patchwork')\n")
  # Individual plots instead
  print(plot_selection_frequency(fit_full, top_n = 8))
  print(plot_validation_results(fit_full))
}

# =============================================================================
# 5. Parallel Processing for Bootstrap
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("5. PARALLEL PROCESSING FOR BOOTSTRAP\n")
cat(strrep("=", 60), "\n")

# Check if future.apply is available
if (requireNamespace("future.apply", quietly = TRUE) && 
    requireNamespace("future", quietly = TRUE)) {
  
  library(future)
  
  cat("Setting up parallel backend with 2 workers...\n")
  plan(multisession, workers = 2)
  
  # Time comparison
  cat("\n--- Sequential Bootstrap ---\n")
  t1 <- system.time({
    fit_seq <- sglwqs(
      X = wqs_data[, c(PCBs, phthalates)],
      y = wqs_data$y,
      groups = list(PCBs = PCBs, Phthalates = phthalates),
      bootstrap = TRUE,
      n_boot = 20,
      parallel = FALSE,  # Sequential
      seed = 123,
      verbose = FALSE
    )
  })
  cat("Sequential time:", round(t1["elapsed"], 2), "seconds\n")
  
  cat("\n--- Parallel Bootstrap ---\n")
  t2 <- system.time({
    fit_par <- sglwqs(
      X = wqs_data[, c(PCBs, phthalates)],
      y = wqs_data$y,
      groups = list(PCBs = PCBs, Phthalates = phthalates),
      bootstrap = TRUE,
      n_boot = 20,
      parallel = TRUE,  # Parallel
      seed = 123,
      verbose = FALSE
    )
  })
  cat("Parallel time:", round(t2["elapsed"], 2), "seconds\n")
  cat("Speedup:", round(t1["elapsed"] / t2["elapsed"], 2), "x\n")
  
  # Reset to sequential
  plan(sequential)
  cat("\nParallel backend reset to sequential.\n")
  
} else {
  cat("Parallel processing requires 'future' and 'future.apply' packages.\n")
  cat("Install with: install.packages(c('future', 'future.apply'))\n")
  cat("\nExample usage:\n")
  cat("  library(future)\n")
  cat("  plan(multisession, workers = 4)  # Use 4 cores\n")
  cat("  fit <- sglwqs(..., bootstrap = TRUE, parallel = TRUE)\n")
  cat("  plan(sequential)  # Reset when done\n")
}

# =============================================================================
# 6. Extracting Results for Publication
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("6. EXTRACTING RESULTS FOR PUBLICATION\n")
cat(strrep("=", 60), "\n")

# Using the full model from section 4
cat("Extracting publishable results from fit_full...\n\n")

# 6.1 Weight table with bootstrap info
cat("--- Weight Table with Bootstrap Statistics ---\n")
weights_df <- extract_weights(fit_full, include_bootstrap = TRUE)
weights_df <- weights_df[weights_df$weight > 0, ]
weights_df <- weights_df[order(-weights_df$selection_freq), ]
print(head(weights_df, 10))

# 6.2 Group-wise summary
cat("\n--- Group-wise Summary ---\n")
summary_by_group(fit_full)

# 6.3 Validation p-values for abstract/results section
cat("\n--- Key Statistics for Reporting ---\n")
cat("Positive mixture effect:\n")
cat(sprintf("  beta = %.3f (95%% CI: %.3f to %.3f), p = %s\n",
            fit_full$validation_info$wqs_pos_estimate,
            fit_full$validation_info$wqs_pos_estimate - 1.96 * fit_full$validation_info$wqs_pos_se,
            fit_full$validation_info$wqs_pos_estimate + 1.96 * fit_full$validation_info$wqs_pos_se,
            format.pval(fit_full$validation_info$wqs_pos_pvalue, digits = 3)))

cat("\nNegative mixture effect:\n")
cat(sprintf("  beta = %.3f (95%% CI: %.3f to %.3f), p = %s\n",
            fit_full$validation_info$wqs_neg_estimate,
            fit_full$validation_info$wqs_neg_estimate - 1.96 * fit_full$validation_info$wqs_neg_se,
            fit_full$validation_info$wqs_neg_estimate + 1.96 * fit_full$validation_info$wqs_neg_se,
            format.pval(fit_full$validation_info$wqs_neg_pvalue, digits = 3)))

# =============================================================================
# Summary: Recommended Analysis Workflow
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("RECOMMENDED ANALYSIS WORKFLOW\n")
cat(strrep("=", 60), "\n")
cat("
1. EXPLORATORY ANALYSIS
   fit <- sglwqs(X, y, covariates, groups = groups)
   summary(fit)
   plot_weights(fit)

2. STABLE WEIGHTS (Bootstrap)
   fit <- sglwqs(..., bootstrap = TRUE, n_boot = 100)
   summary_bootstrap(fit)
   plot_selection_frequency(fit)

3. VALID INFERENCE (Validation)
   fit <- sglwqs(..., validation = TRUE, train_prop = 0.6)
   summary_validation(fit)
   plot_validation_results(fit)

4. PUBLICATION-READY (Both)
   library(future)
   plan(multisession, workers = 4)
   
   fit <- sglwqs(
     X = exposures,
     y = outcome,
     covariates = confounders,
     groups = list(Group1 = vars1, Group2 = vars2),
     group_structure = 'direction',
     bootstrap = TRUE,
     n_boot = 100,
     parallel = TRUE,
     validation = TRUE,
     train_prop = 0.6,
     seed = 123
   )
   
   plan(sequential)
   
   # Results
   get_inference_table(fit)
   plot_combined_results(fit)

5. EXPORT FOR PUBLICATION
   weights_df <- extract_weights(fit, include_bootstrap = TRUE)
   write.csv(weights_df, 'table_s1_weights.csv')
   
   ggsave('figure_2_selection_freq.png', 
          plot_selection_frequency(fit), 
          width = 8, height = 6)
")

# =============================================================================
# 7. Multiple Imputation with mice
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("7. MULTIPLE IMPUTATION WITH mice\n")
cat(strrep("=", 60), "\n")

if (requireNamespace("mice", quietly = TRUE)) {
  library(mice)
  
  # Create data with missing values
  cat("Creating data with missing values...\n")
  wqs_data_miss <- wqs_data[, c("y", "sex", PCBs, phthalates)]
  
  # Insert random missingness (10%)
  set.seed(123)
  for (var in c(PCBs, phthalates)) {
    miss_idx <- sample(nrow(wqs_data_miss), floor(nrow(wqs_data_miss) * 0.1))
    wqs_data_miss[miss_idx, var] <- NA
  }
  
  cat("Missing data pattern:\n")
  cat("  Total observations:", nrow(wqs_data_miss), "\n")
  cat("  Complete cases:", sum(complete.cases(wqs_data_miss)), "\n")
  cat("  Missing rate per variable: ~10%\n\n")
  
  # Multiple imputation
  cat("Running multiple imputation (m=5)...\n")
  imp <- mice(wqs_data_miss, m = 5, method = "pmm", printFlag = FALSE, seed = 123)
  
  cat("Imputation complete.\n\n")
  
  # Run sglwqs_mice
  cat("Fitting sglwqs with multiple imputation...\n")
  fit_mice <- sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c(PCBs, phthalates),
    outcome_var = "y",
    covariate_vars = "sex",
    groups = list(PCBs = PCBs, Phthalates = phthalates),
    validation = TRUE,
    train_prop = 0.6,
    seed = 123,
    verbose = TRUE
  )
  
  # Display results
  print(fit_mice)
  summary(fit_mice)
  
  # Extract pooled weights
  cat("\n--- Pooled Weights ---\n")
  pooled_weights <- extract_pooled_weights(fit_mice)
  pooled_weights <- pooled_weights[pooled_weights$weight > 0, ]
  print(head(pooled_weights[order(-pooled_weights$weight), ], 10))
  
  # Visualize pooled weights
  cat("\n--- Generating Pooled Weights Plot ---\n")
  p_pooled <- plot_pooled_weights(fit_mice, top_n = 10)
  print(p_pooled)
  
} else {
  cat("Package 'mice' is not installed.\n")
  cat("Install with: install.packages('mice')\n\n")
  cat("Example usage:\n")
  cat("
  library(mice)
  
  # Multiple imputation
  imp <- mice(data_with_missing, m = 5)
  
  # Fit sglwqs with MI
  fit <- sglwqs_mice(
    mids_obj = imp,
    exposure_vars = c('x1', 'x2', 'x3'),
    outcome_var = 'y',
    covariate_vars = c('cov1', 'cov2'),
    groups = list(GroupA = c('x1', 'x2'), GroupB = c('x3')),
    validation = TRUE,
    seed = 123
  )
  
  # Results
  summary(fit)
  plot_pooled_weights(fit)
  
  # Rubin's rules p-values
  fit$pooled$validation$wqs_pos$p_value
  fit$pooled$validation$wqs_neg$p_value
")
}

# =============================================================================
# Summary: Complete Analysis Workflow with Missing Data
# =============================================================================
cat("\n", strrep("=", 60), "\n")
cat("COMPLETE ANALYSIS WORKFLOW WITH MISSING DATA\n")
cat(strrep("=", 60), "\n")
cat("
1. HANDLE MISSING DATA
   library(mice)
   imp <- mice(data, m = 10, method = 'pmm', seed = 123)

2. FIT MODEL WITH MI
   fit <- sglwqs_mice(
     mids_obj = imp,
     exposure_vars = exposure_names,
     outcome_var = 'outcome',
     covariate_vars = covariate_names,
     groups = list(Group1 = vars1, Group2 = vars2),
     validation = TRUE,
     bootstrap = TRUE,
     n_boot = 50,
     seed = 123
   )

3. EXAMINE RESULTS
   summary(fit)                      # Full summary with Rubin's rules
   extract_pooled_weights(fit)       # Pooled weights table
   plot_pooled_weights(fit)          # Visualization

4. REPORT KEY STATISTICS
   # Pooled p-values (Rubin's rules)
   fit$pooled$validation$wqs_pos$p_value
   fit$pooled$validation$wqs_neg$p_value
   
   # Pooled estimates and 95% CI
   pos <- fit$pooled$validation$wqs_pos
   cat('Positive effect: ', pos$estimate, ' (', pos$ci_lower, ', ', pos$ci_upper, ')')
   
   # Fraction of Missing Information
   pos$fmi  # Higher values indicate more information lost to missingness

5. SENSITIVITY ANALYSIS
   # Compare complete case analysis
   fit_cc <- sglwqs(X[complete.cases(X), ], y[complete.cases(X)], ...)
   
   # Compare with different imputation methods
   imp_norm <- mice(data, method = 'norm')
   fit_norm <- sglwqs_mice(imp_norm, ...)
")
