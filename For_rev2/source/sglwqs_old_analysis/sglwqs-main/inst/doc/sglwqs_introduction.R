## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, 
  comment = "#>", 
  eval = FALSE,
  fig.width = 10,
  fig.height = 7,
  out.width = "100%",
  dpi = 150
)

## ----install------------------------------------------------------------------
# install.packages("sglwqs_0.8.10.tar.gz", repos = NULL, type = "source")

## ----workflow-----------------------------------------------------------------
# library(sglwqs)
# library(future)
# 
# # 1. Fit model
# plan(multisession, workers = 4)
# 
# fit <- sglwqs(
#   X = exposures,
#   y = outcome,
#   groups = list(PCBs = pcb_vars, Metals = metal_vars),
#   bootstrap = TRUE,
#   n_boot = 500,
#   validation = TRUE,
#   parallel = TRUE,
#   seed = 123
# )
# 
# plan(sequential)
# 
# # 2. View results
# summary(fit)
# 
# # 3. Visualize
# plot_weights(fit)
# plot_combined_results(fit)
# 
# # 4. Export for publication
# inference_table(fit)

## ----example------------------------------------------------------------------
# library(sglwqs)
# library(groupWQS)
# library(future)
# 
# data("simdata", package = "groupWQS")
# 
# # Define chemical groups
# groups <- list(
#   PCBs = c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_192"),
#   Metals = c("as", "cu", "pb", "sn"),
#   Pesticides = c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos")
# )
# 
# # Fit model
# plan(multisession, workers = 10)
# 
# fit <- sglwqs(
#   X = simdata[, 1:14],
#   y = simdata$Y,
#   family = "binomial",
#   groups = groups,
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   seed = 123
# )
# 
# plan(sequential)
# 
# # Results
# summary(fit)
# inference_table(fit)

## ----example-plots, fig.height=8----------------------------------------------
# plot_weights(fit)

## ----example-combined, fig.width=12, fig.height=10----------------------------
# plot_combined_results(fit, top_n = 5)

## ----parallel-----------------------------------------------------------------
# library(future)
# plan(multisession, workers = parallel::detectCores() - 1)
# # ... fit model with parallel = TRUE ...
# plan(sequential)

