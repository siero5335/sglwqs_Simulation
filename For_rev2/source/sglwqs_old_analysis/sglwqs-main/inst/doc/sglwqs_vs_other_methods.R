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
N_BOOT <- 500

## ----wqs-gwqs-----------------------------------------------------------------
# library(gWQS)
# data("wqs_data")
# PCBs <- names(wqs_data)[grepl("^LBX", names(wqs_data))]
# 
# # 2-index model: positive AND negative directions
# gwqs_fit <- gwqs(
#   yLBX ~ pwqs + nwqs,  # Both positive and negative WQS indices
#   mix_name = PCBs,
#   data = wqs_data,
#   q = 10,
#   validation = 0.6,
#   b = 500,
#   b1_pos = TRUE,
#   rh = 1,
#   family = "gaussian",
#   seed = 2016
# )
# 
# summary(gwqs_fit)
# gwqs_barplot(gwqs_fit)
# 
# # Extract weights for comparison
# weights_gwqs_pos <- data.frame(
#   variable = gwqs_fit$final_weights$mix_name,
#   weight = gwqs_fit$final_weights$mean_weight_p,
#   direction = "Positive",
#   method = "gWQS"
# )
# weights_gwqs_neg <- data.frame(
#   variable = gwqs_fit$final_weights$mix_name,
#   weight = gwqs_fit$final_weights$mean_weight_n,
#   direction = "Negative",
#   method = "gWQS"
# )

## ----wqs-sglwqs---------------------------------------------------------------
# library(sglwqs)
# library(future)
# 
# plan(multisession, workers = 10)
# 
# sglwqs_fit <- sglwqs(
#   X = wqs_data[, PCBs],
#   y = wqs_data$yLBX,
#   groups = list(PCBs = PCBs),
#   n_quantiles = 10,
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   train_prop = 0.6,
#   seed = 123
# )
# 
# plan(sequential)
# 
# # Results
# summary(sglwqs_fit)
# plot_weights(sglwqs_fit)
# inference_table(sglwqs_fit)

## ----qgcomp_noboot------------------------------------------------------------
# library(qgcomp)
# data("metals", package = "qgcomp")
# 
# Xnm <- c('arsenic', 'barium', 'cadmium', 'calcium', 'chromium', 'copper',
#          'iron', 'lead', 'magnesium', 'manganese', 'mercury', 'selenium',
#          'silver', 'sodium', 'zinc')
# 
# qc_noboot <- qgcomp.glm.noboot(
#   disease_state ~ .,
#   expnms = Xnm,
#   data = metals[, c(Xnm, 'disease_state')],
#   family = binomial(),
#   q = 4
# )
# 
# qc_noboot
# plot(qc_noboot)

## ----qgcomp-------------------------------------------------------------------
# qgcomp_fit <- qgcomp.boot(
#   disease_state ~ .,
#   data = metals[, c("disease_state", Xnm)],
#   family = binomial(),
#   q = 10,
#   B = 500,
#   seed = 123
# )
# 
# summary(qgcomp_fit)
# plot(qgcomp_fit)
# 
# # Extract weights
# weights_qgcomp <- data.frame(
#   variable = c(names(qc_noboot$pos.weights), names(qc_noboot$neg.weights)),
#   weight = c(qc_noboot$pos.weights, qc_noboot$neg.weights),
#   direction = c(rep("Positive", length(qc_noboot$pos.weights)),
#                 rep("Negative", length(qc_noboot$neg.weights))),
#   method = "qgcomp"
# )

## ----qgcomp-sglwqs------------------------------------------------------------
# plan(multisession, workers = 10)
# 
# sglwqs_fit2 <- sglwqs(
#   X = metals[, Xnm],
#   y = metals$disease_state,
#   family = "binomial",
#   groups = list(Metals = Xnm),
#   n_quantiles = 10,
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   train_prop = 0.6,
#   seed = 123
# )
# 
# plan(sequential)
# 
# summary(sglwqs_fit2)
# plot_weights(sglwqs_fit2)
# inference_table(sglwqs_fit2)

## ----groupwqs-----------------------------------------------------------------
# library(groupWQS)
# data("simdata", package = "groupWQS")
# 
# PCBs_g <- c("pcb_118", "pcb_138", "pcb_153", "pcb_180", "pcb_192")
# metals_g <- c("as", "cu", "pb", "sn")
# pesticides <- c("carbaryl", "propoxur", "methoxychlor", "diazinon", "chlorpyrifos")
# 
# all_vars_g <- c(PCBs_g, metals_g, pesticides)
# 
# group_list <- list(PCBs_g, metals_g, pesticides)
# X_gwqs <- make.X(simdata, num.groups = 3, groups = group_list)
# x.s <- make.x.s(simdata, num.groups = 3, groups = group_list)
# 
# n_total <- nrow(simdata)
# n_train <- 700
# Y.train <- simdata$Y[1:n_train]
# Y.valid <- simdata$Y[(n_train + 1):n_total]
# X.train <- X_gwqs[1:n_train, ]
# X.valid <- X_gwqs[(n_train + 1):n_total, ]
# 
# gwqs_group_fit <- gwqs.fit(
#   y = Y.valid,
#   y.train = Y.train,
#   x = X.valid,
#   x.train = X.train,
#   x.s = x.s,
#   B = N_BOOT,
#   n.quantiles = 4,
#   func = "binary"
# )
# 
# print(gwqs_group_fit)
# 
# # Extract weights
# # Extract weights from groupWQS
# extract_weights_groupwqs <- function(fit, var_names, group_list, method_name = "groupWQS") {
# 
#   # Get GLM coefficients (excluding Intercept)
#   glm_coefs <- coef(fit$fit)
#   gwqs_coefs <- glm_coefs[grep("^GWQS", names(glm_coefs))]
# 
#   # Get weights
#   weights_list <- fit$weights
# 
#   # Data frame to store results
#   result <- data.frame()
# 
#   for (i in seq_along(weights_list)) {
#     grp_weights <- weights_list[[i]]
#     grp_coef <- gwqs_coefs[i]
# 
#     # Determine direction based on coefficient sign
#     if (grp_coef >= 0) {
#       direction <- "Positive"
#     } else {
#       direction <- "Negative"
#     }
# 
#     df_grp <- data.frame(
#       variable = names(grp_weights),
#       weight = as.numeric(grp_weights),
#       direction = direction,
#       method = method_name,
#       stringsAsFactors = FALSE
#     )
# 
#     result <- rbind(result, df_grp)
#   }
# 
#   return(result)
# }
# 
# weights_groupwqs <- extract_weights_groupwqs(gwqs_group_fit)

## ----groupwqs-sglwqs, fig.height=9--------------------------------------------
# all_vars_g <- c(PCBs_g, metals_g, pesticides)
# 
# plan(multisession, workers = 10)
# 
# sglwqs_fit3 <- sglwqs(
#   X = simdata[, all_vars_g],
#   y = simdata$Y,
#   family = "binomial",
#   groups = list(
#     PCBs = PCBs_g,
#     Metals = metals_g,
#     Pesticides = pesticides
#   ),
#   n_quantiles = 4,
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   train_prop = 0.7,
#   seed = 123
# )
# 
# plan(sequential)
# 
# summary(sglwqs_fit3)
# plot_weights(sglwqs_fit3)
# inference_table(sglwqs_fit3)

## ----combined-plot, fig.width=12, fig.height=10-------------------------------
# plot_combined_results(sglwqs_fit3)

## ----large-scale-example, eval=FALSE------------------------------------------
# library(sglwqs)
# library(future)
# 
# plan(multisession, workers = 10)
# 
# # For large datasets, sglwqs auto-detects and adjusts memory settings
# fit <- sglwqs(
#   X = large_data[, exposures],  # e.g., N = 100,000
#   y = large_data$outcome,
#   groups = groups,
#   family = "binomial",
#   bootstrap = TRUE,
#   n_boot = 500,
#   parallel = TRUE,
#   validation = TRUE,
#   checkpoint_dir = "sglwqs_checkpoint",    # Save progress periodically
#   checkpoint_interval = 50,                 # Every 50 iterations
#   # future_globals_max_size = NULL,         # Auto-detect (default)
#   # future_globals_max_size = 2000 * 1024^2,  # Or specify manually: 2GB
#   seed = 123
# )
# # Output: "Large dataset detected (250 MB). Auto-setting future.globals.maxSize to 1000 MB"
# 
# plan(sequential)

