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

## ----generate-data------------------------------------------------------------
# library(MASS)
# 
# generate_sglwqs_simdata <- function(n = 100000, seed = 71) {
#   set.seed(seed)
# 
#   # Covariates
#   sex <- rbinom(n, 1, 0.5)
#   age <- pmax(pmin(round(rnorm(n, 50, 12)), 80), 20)
#   bmi <- pmax(pmin(round(22 + 0.05*age + 1.5*sex + rnorm(n, 0, 4), 1), 45), 15)
# 
#   # PCBs (positive effect, high correlation)
#   cor_pcb <- matrix(c(1, .75, .70, .65, .75, 1, .80, .70,
#                       .70, .80, 1, .75, .65, .70, .75, 1), 4)
#   pcb <- mvrnorm(n, c(3.5, 4, 4.5, 3), diag(c(1.5, 1.8, 2, 1.5)) %*% cor_pcb %*% diag(c(1.5, 1.8, 2, 1.5)))
#   colnames(pcb) <- c("pcb_138", "pcb_153", "pcb_180", "pcb_187")
# 
#   # Metals (negative effect)
#   metal <- exp(mvrnorm(n, log(c(0.5, 2, 1, 100, 20)), diag(c(0.6, 0.8, 0.5, 0.3, 0.25))))
#   colnames(metal) <- c("cd", "pb", "hg", "se", "mg")
# 
#   # PFASs (no effect)
#   pfas <- pmax(mvrnorm(n, c(3, 1.5, 5), diag(c(2, 1, 3.5))), 0.1)
#   colnames(pfas) <- c("PFOA", "PFNA", "PFOS")
# 
#   # Generate outcome (true coefficients)
#   beta_pcb <- c(0.15, 0.20, 0.08, 0.05)     # Positive
#   beta_metal <- c(-0.08, -0.15, -0.05, -0.03, -0.02)  # Negative
#   beta_pfas <- c(0, 0, 0)                   # None
# 
#   eta <- -3.0 + scale(pcb) %*% beta_pcb + scale(metal) %*% beta_metal +
#          scale(pfas) %*% beta_pfas + sex*0.3 + (age-50)/10*0.02 + (bmi-25)/5*0.05
# 
#   Y <- rbinom(n, 1, 1/(1 + exp(-eta)))
# 
#   data.frame(pcb, metal, pfas, sex = factor(sex), age = age, bmi = bmi, Y = Y)
# }
# 
# simdata <- generate_sglwqs_simdata(n = 100000, seed = 71)

## ----analysis-----------------------------------------------------------------
# library(sglwqs)
# library(future)
# 
# groups <- list(
#   PCBs = c("pcb_138", "pcb_153", "pcb_180", "pcb_187"),
#   Metals = c("cd", "pb", "hg", "se", "mg"),
#   PFASs = c("PFOA", "PFNA", "PFOS")
# )
# 
# plan(multisession, workers = 10)
# 
# # For large-scale analysis, use checkpoint_dir to enable recovery from crashes
# # and reduce memory pressure
# fit <- sglwqs(
#   X = simdata[, unlist(groups)],
#   y = simdata$Y,
#   covariates = simdata[, c("sex", "age", "bmi")],
#   groups = groups,
#   family = "binomial",
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   train_prop = 0.6,
#   checkpoint_dir = "sglwqs_checkpoint",  # Save progress every 50 iterations
#   checkpoint_interval = 50,               # Interval between saves
#   seed = 123
# )
# 
# plan(sequential)

## ----checkpoint-usage, eval=FALSE---------------------------------------------
# # Example: Resume from crashed analysis
# fit <- sglwqs(
#   X = X, y = y,
#   groups = groups,
#   bootstrap = TRUE,
#   n_boot = 500,
#   checkpoint_dir = "sglwqs_checkpoint",  # Same directory as before
#   seed = 123  # Same seed for reproducibility
# )
# # Analysis will resume from where it left off

## ----results------------------------------------------------------------------
# # Summary
# summary(fit)
# 
# # Publication-ready table
# inference_table(fit)

## ----results-plot, fig.height=9-----------------------------------------------
# # Weight visualization
# plot_weights(fit)

## ----combined-plot, fig.width=12, fig.height=10-------------------------------
# # Combined visualization
# plot_combined_results(fit, top_n = 5)

