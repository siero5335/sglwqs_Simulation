## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>", 
  eval = FALSE,
  fig.width = 10,
  fig.height = 7,
  out.width = "100%",
  dpi = 150
)

## ----workflow-----------------------------------------------------------------
# library(sglwqs)
# library(mice)
# library(future)
# 
# # 1. Multiple imputation with mice
# imp <- mice(data, m = 20, seed = 123)
# 
# # 2. Fit SGL-WQS to imputed data
# plan(multisession, workers = 4)
# 
# fit_mice <- sglwqs_mice(
#   mids_obj = imp,
#   exposure_vars = exposure_vars,
#   outcome_var = "Y",
#   covariate_vars = c("sex", "age", "bmi"),
#   groups = groups,
#   family = "binomial",
#   bootstrap = TRUE,
#   n_boot = 200,
#   validation = TRUE,
#   parallel = TRUE,
#   seed = 123
# )
# 
# plan(sequential)
# 
# # 3. View results (same functions as regular sglwqs!)
# summary(fit_mice)
# inference_table(fit_mice)

## ----workflow-plots, fig.height=8---------------------------------------------
# plot_weights(fit_mice, top_n = 5)

## ----workflow-combined, fig.width=12, fig.height=10---------------------------
# plot_combined_results(fit_mice)

## ----unified-api--------------------------------------------------------------
# # Regular analysis
# fit <- sglwqs(X, y, ...)
# summary(fit)
# plot_weights(fit)
# inference_table(fit)
# 
# # MICE analysis (same interface!)
# fit_mice <- sglwqs_mice(imp, ...)
# summary(fit_mice)
# plot_weights(fit_mice)
# inference_table(fit_mice)

## ----example------------------------------------------------------------------
# # Generate data with missing values
# set.seed(123)
# n <- 5000
# data <- data.frame(
#   x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
#   age = round(rnorm(n, 50, 12)),
#   bmi = round(rnorm(n, 25, 4), 1),
#   Y = rbinom(n, 1, 0.2)
# )
# 
# # Introduce ~10% missing in covariates
# data$age[sample(n, n*0.1)] <- NA
# data$bmi[sample(n, n*0.1)] <- NA
# 
# # Multiple imputation
# library(mice)
# imp <- mice(data, m = 20, method = "pmm", seed = 123, printFlag = FALSE)
# 
# # Fit SGL-WQS
# library(sglwqs)
# library(future)
# 
# plan(multisession, workers = 4)
# 
# fit_mice <- sglwqs_mice(
#   mids_obj = imp,
#   exposure_vars = c("x1", "x2", "x3"),
#   outcome_var = "Y",
#   covariate_vars = c("age", "bmi"),
#   family = "binomial",
#   bootstrap = TRUE,
#   n_boot = 200,
#   validation = TRUE,
#   parallel = TRUE,
#   # For large-scale analysis with many imputations:
#   checkpoint_dir = "sglwqs_mice_checkpoint",  # Save progress
#   seed = 456
# )
# 
# plan(sequential)
# 
# # Results
# summary(fit_mice)
# inference_table(fit_mice)

## ----large-scale-mice, eval=FALSE---------------------------------------------
# # Each imputation's bootstrap progress is saved separately
# fit_mice <- sglwqs_mice(
#   mids_obj = imp,  # m=20 imputations
#   exposure_vars = exposures,
#   outcome_var = "y",
#   groups = groups,
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   checkpoint_dir = "sglwqs_mice_checkpoint",
#   checkpoint_interval = 50,
#   seed = 123
# )
# 
# # Checkpoint structure:
# # sglwqs_mice_checkpoint/
# # ├── imp_1/bootstrap_checkpoint.rds
# # ├── imp_2/bootstrap_checkpoint.rds
# # ...
# # └── imp_20/bootstrap_checkpoint.rds

## ----example-plot-------------------------------------------------------------
# plot_weights(fit_mice)

