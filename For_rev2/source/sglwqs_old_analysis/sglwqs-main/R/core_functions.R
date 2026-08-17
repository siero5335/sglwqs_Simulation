#' Internal Function: Fit Single SGL Model
#'
#' Core fitting function used internally by sglwqs.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Processed covariate matrix (can be NULL).
#' @param var_names Exposure variable names.
#' @param cov_names Covariate names.
#' @param groups Group list.
#' @param group_by_compound Logical.
#' @param group_structure Character.
#' @param penalize_covariates Logical.
#' @param family Character.
#' @param lambda Character or numeric.
#' @param nfolds Integer.
#' @param obs_weights Optional numeric observation weights for selection loss.
#' @param ... Additional arguments to cv.sparsegl.
#'
#' @return A list with fitted model results.
#'
#' @keywords internal
fit_sgl_core <- function(X_quantile, y, cov_matrix, var_names, cov_names,
                          groups, group_by_compound, group_structure,
                          penalize_covariates, family, lambda, nfolds,
                          obs_weights = NULL, ...) {
  
  n <- nrow(X_quantile)
  p <- ncol(X_quantile)
  q <- if (is.null(cov_matrix)) 0 else ncol(cov_matrix)
  obs_weights <- .validate_obs_weights(obs_weights, n)
  sparsegl_family <- .coerce_sparsegl_family(family, obs_weights)
  
  # Create group mapping (for compound groups)
  if (!is.null(groups) && group_by_compound) {
    var_to_group <- rep(NA_integer_, p)
    for (i in seq_along(groups)) {
      var_idx <- which(var_names %in% groups[[i]])
      var_to_group[var_idx] <- i
    }
    unassigned <- which(is.na(var_to_group))
    if (length(unassigned) > 0) {
      next_group <- length(groups) + 1
      for (idx in unassigned) {
        var_to_group[idx] <- next_group
        next_group <- next_group + 1
      }
    }
    n_compound_groups <- max(var_to_group)
  } else {
    var_to_group <- NULL
    n_compound_groups <- NULL
  }
  
  # Build design matrix
  if (q > 0) {
    design_matrix <- cbind(X_quantile, -X_quantile, cov_matrix)
    
    if (!is.null(groups) && group_by_compound) {
      # group_structure = "direction": treat positive/negative as separate groups
      group_pos <- var_to_group
      group_neg <- var_to_group + n_compound_groups
      group_cov <- seq(2 * n_compound_groups + 1, 2 * n_compound_groups + q)
      group <- c(group_pos, group_neg, group_cov)
      
      n_total_groups <- 2 * n_compound_groups + q
      if (penalize_covariates) {
        pf_sparse <- rep(1, 2 * p + q)
        pf_group <- rep(1, n_total_groups)
        lower_bnd <- c(rep(0, 2 * n_compound_groups), rep(-Inf, q))
      } else {
        pf_sparse <- c(rep(1, 2 * p), rep(0, q))
        pf_group <- c(rep(1, 2 * n_compound_groups), rep(0, q))
        lower_bnd <- c(rep(0, 2 * n_compound_groups), rep(-Inf, q))
      }
    } else {
      # No user groups: G=1 direction-split (manuscript model)
      # All positive copies in group 1, all negative copies in group 2
      group <- c(rep(1L, p), rep(2L, p), seq(3L, 2L + q))
      if (penalize_covariates) {
        pf_sparse <- rep(1, 2 * p + q)
        pf_group <- rep(1, 2 + q)
        lower_bnd <- c(rep(0, 2), rep(-Inf, q))
      } else {
        pf_sparse <- c(rep(1, 2 * p), rep(0, q))
        pf_group <- c(rep(1, 2), rep(0, q))
        lower_bnd <- c(rep(0, 2), rep(-Inf, q))
      }
    }
  } else {
    design_matrix <- cbind(X_quantile, -X_quantile)

    if (!is.null(groups) && group_by_compound) {
      # group_structure = "direction": treat positive/negative as separate groups
      group_pos <- var_to_group
      group_neg <- var_to_group + n_compound_groups
      group <- c(group_pos, group_neg)
      pf_sparse <- rep(1, 2 * p)
      pf_group <- rep(1, 2 * n_compound_groups)
      lower_bnd <- rep(0, 2 * n_compound_groups)
    } else {
      # No user groups: G=1 direction-split (manuscript model)
      group <- c(rep(1L, p), rep(2L, p))
      pf_sparse <- rep(1, 2 * p)
      pf_group <- rep(1, 2)
      lower_bnd <- rep(0, 2)
    }
  }
  
  # Fit Sparse Group Lasso (with error handling)
  fit <- tryCatch({
    sparsegl::cv.sparsegl(
      x = design_matrix,
      y = y,
      group = group,
      family = sparsegl_family,
      pf_sparse = pf_sparse,
      pf_group = pf_group,
      lower_bnd = lower_bnd,
      nfolds = nfolds,
      weights = obs_weights,
      ...
    )
  }, error = function(e) {
    # Make error messages more informative
    error_msg <- e$message
    
    if (grepl("singular|rank|collinear", error_msg, ignore.case = TRUE)) {
      stop("SGL fitting failed due to multicollinearity or rank deficiency.\n",
           "  Possible causes:\n",
           "  - Too many exposures relative to sample size (p >> n)\n",
           "  - Highly correlated variables in the design matrix\n",
           "  Try: Reduce the number of exposures or increase sample size.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else if (grepl("convergence|iterate", error_msg, ignore.case = TRUE)) {
      stop("SGL fitting failed to converge.\n",
           "  Possible causes:\n",
           "  - Data may have unusual scaling or outliers\n",
           "  - Lambda sequence may be inappropriate\n",
           "  Try: Standardize your data or adjust lambda settings.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else if (grepl("fold|cv|cross", error_msg, ignore.case = TRUE)) {
      stop("Cross-validation failed.\n",
           "  Possible causes:\n",
           "  - Sample size too small for ", nfolds, "-fold CV\n",
           "  - Some folds may have too few observations\n",
           "  Try: Reduce nfolds (e.g., nfolds = 5) or increase sample size.\n",
           "  Original error: ", error_msg, call. = FALSE)
    } else {
      stop("SGL fitting failed.\n",
           "  Sample size: ", n, ", Number of exposures: ", p, "\n",
           "  This may occur when:\n",
           "  - Sample size is too small relative to number of variables\n",
           "  - Data contains NA, Inf, or unusual values\n",
           "  Original error: ", error_msg, call. = FALSE)
    }
  })
  
  # Extract coefficients
  all_coef <- coef(fit, s = lambda)
  
  intercept <- all_coef[1]
  pos_coef <- as.numeric(all_coef[2:(p + 1)])
  neg_coef <- as.numeric(all_coef[(p + 2):(2 * p + 1)])
  
  names(pos_coef) <- var_names
  names(neg_coef) <- var_names
  
  # Covariate coefficients
  if (q > 0) {
    cov_coef <- as.numeric(all_coef[(2 * p + 2):length(all_coef)])
    names(cov_coef) <- cov_names
  } else {
    cov_coef <- NULL
  }
  
  return(list(
    fit = fit,
    coefficients = all_coef,
    pos_coef = pos_coef,
    neg_coef = neg_coef,
    cov_coef = cov_coef,
    intercept = intercept,
    design_matrix = design_matrix,
    group = group
  ))
}


#' Internal Function: Bootstrap Aggregation for Weights
#'
#' Performs bootstrap aggregation to stabilize weight estimates.
#' Supports parallel processing via the future framework.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Processed covariate matrix.
#' @param var_names Exposure variable names.
#' @param cov_names Covariate names.
#' @param groups Group list.
#' @param group_by_compound Logical.
#' @param group_structure Character.
#' @param penalize_covariates Logical.
#' @param family Character.
#' @param lambda Character or numeric.
#' @param nfolds Integer.
#' @param n_boot Number of bootstrap iterations.
#' @param seed Random seed.
#' @param verbose Logical.
#' @param parallel Logical. Whether to use parallel processing (requires future.apply).
#' @param ... Additional arguments.
#'
#' @return A list with aggregated results.
#'
#' @keywords internal
bootstrap_sgl <- function(X_quantile, y, cov_matrix, var_names, cov_names,
                           groups, group_by_compound, group_structure,
                           penalize_covariates, family, lambda, nfolds,
                           n_boot, seed, verbose, parallel = FALSE,
                           checkpoint_dir = NULL, checkpoint_interval = 50,
                           cleanup_checkpoint = TRUE, 
                           stratified = TRUE, obs_weights = NULL, ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(X_quantile)
  p <- ncol(X_quantile)
  q_cov <- length(cov_names)
  obs_weights <- .validate_obs_weights(obs_weights, n)
  
  # Prepare checkpoint functionality
  use_checkpoint <- !is.null(checkpoint_dir)
  checkpoint_file <- NULL
  start_boot <- 1
  boot_pos_coef <- NULL
  boot_neg_coef <- NULL
  boot_cov_coef <- NULL
  boot_pos_index_sum <- NULL
  boot_neg_index_sum <- NULL
  boot_pos_index_sum_by_group <- NULL
  boot_neg_index_sum_by_group <- NULL
  boot_success <- NULL
  completed_boots <- integer(0)
  
  if (use_checkpoint) {
    if (!dir.exists(checkpoint_dir)) {
      dir.create(checkpoint_dir, recursive = TRUE)
    }
    checkpoint_file <- file.path(checkpoint_dir, "bootstrap_checkpoint.rds")
    
    # Load existing checkpoint if available
    if (file.exists(checkpoint_file)) {
      if (verbose) message("Loading checkpoint from ", checkpoint_file)
      checkpoint_data <- readRDS(checkpoint_file)
      
      # Validate checkpoint data
      if (!is.null(checkpoint_data$boot_pos_coef) && 
          ncol(checkpoint_data$boot_pos_coef) == p &&
          checkpoint_data$n_boot == n_boot) {
        boot_pos_coef <- checkpoint_data$boot_pos_coef
        boot_neg_coef <- checkpoint_data$boot_neg_coef
        boot_cov_coef <- checkpoint_data$boot_cov_coef
        boot_pos_index_sum <- checkpoint_data$boot_pos_index_sum
        boot_neg_index_sum <- checkpoint_data$boot_neg_index_sum
        boot_pos_index_sum_by_group <- checkpoint_data$boot_pos_index_sum_by_group
        boot_neg_index_sum_by_group <- checkpoint_data$boot_neg_index_sum_by_group
        boot_success <- checkpoint_data$boot_success
        completed_boots <- checkpoint_data$completed_boots
        start_boot <- if (length(completed_boots) > 0) max(completed_boots) + 1 else 1
        if (verbose) message("Resuming from iteration ", start_boot, 
                            " (", length(completed_boots), " completed)")
      } else {
        if (verbose) message("Checkpoint incompatible, starting fresh")
      }
    }
  }
  
  # Initialize result matrices (if not resuming from checkpoint)
  if (is.null(boot_pos_coef)) {
    boot_pos_coef <- matrix(0, nrow = n_boot, ncol = p)
    boot_neg_coef <- matrix(0, nrow = n_boot, ncol = p)
    colnames(boot_pos_coef) <- var_names
    colnames(boot_neg_coef) <- var_names
  }
  if (q_cov > 0 && is.null(boot_cov_coef)) {
    boot_cov_coef <- matrix(NA_real_, nrow = n_boot, ncol = q_cov)
    colnames(boot_cov_coef) <- cov_names
  } else if (q_cov == 0) {
    boot_cov_coef <- NULL
  }
  if (is.null(boot_pos_index_sum)) {
    boot_pos_index_sum <- rep(NA_real_, n_boot)
  }
  if (is.null(boot_neg_index_sum)) {
    boot_neg_index_sum <- rep(NA_real_, n_boot)
  }
  if (!is.null(groups)) {
    group_names <- names(groups)
    if (is.null(boot_pos_index_sum_by_group)) {
      boot_pos_index_sum_by_group <- stats::setNames(
        lapply(group_names, function(...) rep(NA_real_, n_boot)),
        group_names
      )
    }
    if (is.null(boot_neg_index_sum_by_group)) {
      boot_neg_index_sum_by_group <- stats::setNames(
        lapply(group_names, function(...) rep(NA_real_, n_boot)),
        group_names
      )
    }
  } else {
    boot_pos_index_sum_by_group <- NULL
    boot_neg_index_sum_by_group <- NULL
  }

  # Bootstrap success tracking vector (restore from checkpoint or create new)
  if (is.null(boot_success)) {
    boot_success <- logical(n_boot)
  }
  
  # Remaining bootstrap indices
  remaining_boots <- setdiff(seq_len(n_boot), completed_boots)
  
  if (length(remaining_boots) == 0) {
    if (verbose) message("All bootstrap iterations already completed")
  } else {
    
    # Single bootstrap execution function
    run_single_boot <- function(b, X_quantile, y, cov_matrix, var_names, cov_names,
                                 groups, group_by_compound, group_structure,
                                 penalize_covariates, family, lambda, nfolds,
                                 stratified = TRUE, obs_weights = NULL, ...) {
      
      n <- nrow(X_quantile)
      
      # Resampling (stratified or simple)
      if (stratified && family == "binomial") {
        # Stratified bootstrap: maintain case/control ratio
        y_levels <- unique(y)
        boot_idx <- integer(0)
        
        for (level in y_levels) {
          level_idx <- which(y == level)
          # Resample equal number from each stratum
          boot_idx <- c(boot_idx, sample(level_idx, length(level_idx), replace = TRUE))
        }
        
        # Shuffle (randomize order)
        boot_idx <- sample(boot_idx)
        
      } else {
        # Simple bootstrap
        boot_idx <- sample(n, replace = TRUE)
      }
      
      X_boot <- X_quantile[boot_idx, , drop = FALSE]
      y_boot <- y[boot_idx]
      cov_boot <- if (!is.null(cov_matrix)) cov_matrix[boot_idx, , drop = FALSE] else NULL
      weights_boot <- if (!is.null(obs_weights)) obs_weights[boot_idx] else NULL
      
      # Fitting (with error handling)
      tryCatch({
        fit_result <- fit_sgl_core(
          X_quantile = X_boot,
          y = y_boot,
          cov_matrix = cov_boot,
          var_names = var_names,
          cov_names = cov_names,
          groups = groups,
          group_by_compound = group_by_compound,
          group_structure = group_structure,
          penalize_covariates = penalize_covariates,
          family = family,
          lambda = lambda,
          nfolds = nfolds,
          obs_weights = weights_boot,
          ...
        )

        weight_result <- calculate_weights(
          pos_coef = fit_result$pos_coef,
          neg_coef = fit_result$neg_coef,
          var_names = var_names,
          groups = groups,
          handle_collinearity = "net"
        )
        
        return(list(
          success = TRUE,
          pos_coef = fit_result$pos_coef,
          neg_coef = fit_result$neg_coef,
          cov_coef = fit_result$cov_coef,
          pos_index_sum = weight_result$pos_index_sum,
          neg_index_sum = weight_result$neg_index_sum,
          pos_index_sum_by_group = weight_result$pos_index_sum_by_group,
          neg_index_sum_by_group = weight_result$neg_index_sum_by_group
        ))
        
      }, error = function(e) {
        return(list(
          success = FALSE,
          pos_coef = rep(0, length(var_names)),
          neg_coef = rep(0, length(var_names)),
          cov_coef = if (!is.null(cov_names)) {
            stats::setNames(rep(NA_real_, length(cov_names)), cov_names)
          } else {
            NULL
          },
          pos_index_sum = NA_real_,
          neg_index_sum = NA_real_,
          pos_index_sum_by_group = if (!is.null(groups)) {
            stats::setNames(as.list(rep(NA_real_, length(groups))), names(groups))
          } else {
            NULL
          },
          neg_index_sum_by_group = if (!is.null(groups)) {
            stats::setNames(as.list(rep(NA_real_, length(groups))), names(groups))
          } else {
            NULL
          },
          error_msg = e$message
        ))
      })
    }
    
    # Checkpoint save function
    append_completed_boots <- function(indices) {
      if (length(indices) == 0) {
        return(invisible(NULL))
      }
      # indices are only added once so unique is unnecessary (avoids O(n^2))
      completed_boots <<- c(completed_boots, indices)
      invisible(NULL)
    }
    
    save_checkpoint <- function() {
      if (use_checkpoint) {
        checkpoint_data <- list(
          boot_pos_coef = boot_pos_coef,
          boot_neg_coef = boot_neg_coef,
          boot_cov_coef = boot_cov_coef,
          boot_pos_index_sum = boot_pos_index_sum,
          boot_neg_index_sum = boot_neg_index_sum,
          boot_pos_index_sum_by_group = boot_pos_index_sum_by_group,
          boot_neg_index_sum_by_group = boot_neg_index_sum_by_group,
          boot_success = boot_success,
          completed_boots = completed_boots,
          n_boot = n_boot,
          var_names = var_names,
          timestamp = Sys.time()
        )
        saveRDS(checkpoint_data, checkpoint_file)
        if (verbose) {
          if (length(completed_boots) > 0) {
            message("  Checkpoint saved at iteration ", max(completed_boots))
          } else {
            message("  Checkpoint saved with no successful iterations yet")
          }
        }
      }
    }
    
    # Execute parallel processing (checkpoint per batch)
    if (parallel) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        warning("Package 'future.apply' is not installed. Falling back to sequential processing.\n",
                "Install with: install.packages('future.apply')")
        parallel <- FALSE
      } else {
        if (verbose) {
          message("Running bootstrap aggregation with ", length(remaining_boots), 
                  " remaining iterations (parallel)...")
          if (use_checkpoint) {
            message("Checkpointing every ", checkpoint_interval, " iterations to ", checkpoint_dir)
          }
        }
        
        # Checkpoint via batch processing
        batch_start <- 1
        batch_errors <- character(0)
        
        while (batch_start <= length(remaining_boots)) {
          batch_end <- min(batch_start + checkpoint_interval - 1, length(remaining_boots))
          batch_indices <- remaining_boots[batch_start:batch_end]
          
          # Execute batch (catch batch-level errors too)
          batch_result <- tryCatch({
            batch_results <- future.apply::future_lapply(
              batch_indices,
              function(b) {
                run_single_boot(b, X_quantile, y, cov_matrix, var_names, cov_names,
                                groups, group_by_compound, group_structure,
                                penalize_covariates, family, lambda, nfolds,
                                stratified = stratified,
                                obs_weights = obs_weights,
                                ...)
              },
              future.seed = TRUE
            )
            list(success = TRUE, results = batch_results)
          }, error = function(e) {
            list(success = FALSE, error_msg = e$message)
          })
          
          if (batch_result$success) {
            # Store results
            batch_completed <- integer(0)
            for (i in seq_along(batch_indices)) {
              b <- batch_indices[i]
              if (batch_result$results[[i]]$success) {
                boot_pos_coef[b, ] <- batch_result$results[[i]]$pos_coef
                boot_neg_coef[b, ] <- batch_result$results[[i]]$neg_coef
                if (!is.null(boot_cov_coef) &&
                    !is.null(batch_result$results[[i]]$cov_coef)) {
                  boot_cov_coef[b, ] <- batch_result$results[[i]]$cov_coef
                }
                boot_pos_index_sum[b] <- batch_result$results[[i]]$pos_index_sum
                boot_neg_index_sum[b] <- batch_result$results[[i]]$neg_index_sum
                if (!is.null(boot_pos_index_sum_by_group) &&
                    !is.null(batch_result$results[[i]]$pos_index_sum_by_group)) {
                  for (grp in names(boot_pos_index_sum_by_group)) {
                    boot_pos_index_sum_by_group[[grp]][b] <-
                      batch_result$results[[i]]$pos_index_sum_by_group[[grp]]
                    boot_neg_index_sum_by_group[[grp]][b] <-
                      batch_result$results[[i]]$neg_index_sum_by_group[[grp]]
                  }
                }
                boot_success[b] <- TRUE
                batch_completed <- c(batch_completed, b)
              }
            }
            append_completed_boots(batch_completed)
            
            # Free memory
            rm(batch_result)
            gc(verbose = FALSE)
          } else {
            # Entire batch failed
            batch_errors <- c(batch_errors, batch_result$error_msg)
            if (verbose) {
              message("  Warning: Batch ", batch_start, "-", batch_end, " failed: ", 
                      substr(batch_result$error_msg, 1, 100))
            }
          }
          
          # Save checkpoint
          save_checkpoint()

          batch_start <- batch_end + 1
        }
        
        if (length(batch_errors) > 0 && verbose) {
          message("Note: ", length(batch_errors), " batch(es) failed during parallel processing.")
          message("Unique errors: ", paste(unique(batch_errors)[1:min(3, length(unique(batch_errors)))], collapse = "; "))
        }
      }
    }
    
    # Sequential processing
    if (!parallel) {
      if (verbose) {
        message("Running bootstrap aggregation with ", length(remaining_boots), 
                " remaining iterations...")
        if (use_checkpoint) {
          message("Checkpointing every ", checkpoint_interval, " iterations to ", checkpoint_dir)
        }
        pb <- txtProgressBar(min = 0, max = length(remaining_boots), style = 3)
      }
      
      for (i in seq_along(remaining_boots)) {
        b <- remaining_boots[i]
        
        result <- run_single_boot(
          b, X_quantile, y, cov_matrix, var_names, cov_names,
          groups, group_by_compound, group_structure,
          penalize_covariates, family, lambda, nfolds, 
          stratified = stratified,
          obs_weights = obs_weights,
          ...
        )
        
        if (result$success) {
          boot_pos_coef[b, ] <- result$pos_coef
          boot_neg_coef[b, ] <- result$neg_coef
          if (!is.null(boot_cov_coef) && !is.null(result$cov_coef)) {
            boot_cov_coef[b, ] <- result$cov_coef
          }
          boot_pos_index_sum[b] <- result$pos_index_sum
          boot_neg_index_sum[b] <- result$neg_index_sum
          if (!is.null(boot_pos_index_sum_by_group) && !is.null(result$pos_index_sum_by_group)) {
            for (grp in names(boot_pos_index_sum_by_group)) {
              boot_pos_index_sum_by_group[[grp]][b] <- result$pos_index_sum_by_group[[grp]]
              boot_neg_index_sum_by_group[[grp]][b] <- result$neg_index_sum_by_group[[grp]]
            }
          }
          boot_success[b] <- TRUE
          append_completed_boots(b)
        }
        
        # Save checkpoint
        if (use_checkpoint && (i %% checkpoint_interval == 0 || i == length(remaining_boots))) {
          save_checkpoint()
        }
        
        if (verbose) setTxtProgressBar(pb, i)
      }
      
      if (verbose) close(pb)
    }
  }
  
  # Aggregate results (success determined by boot_success vector — all-zero coefficients included as valid)
  successful_boots <- sum(boot_success)

  if (verbose) {
    message("Completed ", successful_boots, "/", n_boot, " bootstrap iterations.")
  }

  if (successful_boots < n_boot * 0.5) {
    warning(
      sprintf(
        "Successful bootstrap iterations: %d / %d (%.1f%%).",
        successful_boots, n_boot, 100 * successful_boots / n_boot
      ),
      call. = FALSE
    )
  }

  # Mean coefficients (successful bootstraps only — including all-zero results)
  if (successful_boots > 0) {
    mean_pos_coef <- colMeans(boot_pos_coef[boot_success, , drop = FALSE])
    mean_neg_coef <- colMeans(boot_neg_coef[boot_success, , drop = FALSE])

    # Standard errors
    se_pos_coef <- apply(boot_pos_coef[boot_success, , drop = FALSE], 2, sd)
    se_neg_coef <- apply(boot_neg_coef[boot_success, , drop = FALSE], 2, sd)

    # Selection frequency (proportion of non-zero across all successful bootstraps)
    selection_freq_pos <- colMeans(boot_pos_coef[boot_success, , drop = FALSE] > 0)
    selection_freq_neg <- colMeans(boot_neg_coef[boot_success, , drop = FALSE] > 0)

    index_sum_pos_block <- boot_pos_index_sum[boot_success]
    index_sum_neg_block <- boot_neg_index_sum[boot_success]
    mean_index_sum_pos <- mean(index_sum_pos_block, na.rm = TRUE)
    mean_index_sum_neg <- mean(index_sum_neg_block, na.rm = TRUE)
    se_index_sum_pos <- stats::sd(index_sum_pos_block, na.rm = TRUE)
    se_index_sum_neg <- stats::sd(index_sum_neg_block, na.rm = TRUE)

    if (!is.null(boot_pos_index_sum_by_group)) {
      mean_index_sum_by_group_pos <- lapply(boot_pos_index_sum_by_group, function(x) {
        mean(x[boot_success], na.rm = TRUE)
      })
      mean_index_sum_by_group_neg <- lapply(boot_neg_index_sum_by_group, function(x) {
        mean(x[boot_success], na.rm = TRUE)
      })
      se_index_sum_by_group_pos <- lapply(boot_pos_index_sum_by_group, function(x) {
        stats::sd(x[boot_success], na.rm = TRUE)
      })
      se_index_sum_by_group_neg <- lapply(boot_neg_index_sum_by_group, function(x) {
        stats::sd(x[boot_success], na.rm = TRUE)
      })
    } else {
      mean_index_sum_by_group_pos <- NULL
      mean_index_sum_by_group_neg <- NULL
      se_index_sum_by_group_pos <- NULL
      se_index_sum_by_group_neg <- NULL
    }

    if (!is.null(boot_cov_coef)) {
      cov_block <- boot_cov_coef[boot_success, , drop = FALSE]
      mean_cov_coef <- colMeans(cov_block, na.rm = TRUE)
      se_cov_coef <- apply(cov_block, 2, sd, na.rm = TRUE)
      ci_lower_cov <- apply(cov_block, 2, quantile, probs = 0.025, na.rm = TRUE)
      ci_upper_cov <- apply(cov_block, 2, quantile, probs = 0.975, na.rm = TRUE)
    } else {
      mean_cov_coef <- NULL
      se_cov_coef <- NULL
      ci_lower_cov <- NULL
      ci_upper_cov <- NULL
    }
  } else {
    stop("All bootstrap iterations failed.")
  }
  
  # Remove checkpoint file
  if (use_checkpoint && cleanup_checkpoint && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    if (verbose) message("Checkpoint file removed")
  }
  
  return(list(
    mean_pos_coef = mean_pos_coef,
    mean_neg_coef = mean_neg_coef,
    se_pos_coef = se_pos_coef,
    se_neg_coef = se_neg_coef,
    selection_freq_pos = selection_freq_pos,
    selection_freq_neg = selection_freq_neg,
    mean_index_sum_pos = mean_index_sum_pos,
    mean_index_sum_neg = mean_index_sum_neg,
    se_index_sum_pos = se_index_sum_pos,
    se_index_sum_neg = se_index_sum_neg,
    mean_index_sum_by_group_pos = mean_index_sum_by_group_pos,
    mean_index_sum_by_group_neg = mean_index_sum_by_group_neg,
    se_index_sum_by_group_pos = se_index_sum_by_group_pos,
    se_index_sum_by_group_neg = se_index_sum_by_group_neg,
    boot_pos_coef = boot_pos_coef,
    boot_neg_coef = boot_neg_coef,
    mean_cov_coef = mean_cov_coef,
    se_cov_coef = se_cov_coef,
    ci_lower_cov = ci_lower_cov,
    ci_upper_cov = ci_upper_cov,
    boot_cov_coef = boot_cov_coef,
    boot_pos_index_sum = boot_pos_index_sum,
    boot_neg_index_sum = boot_neg_index_sum,
    boot_pos_index_sum_by_group = boot_pos_index_sum_by_group,
    boot_neg_index_sum_by_group = boot_neg_index_sum_by_group,
    boot_success = boot_success,
    n_successful = successful_boots
  ))
}


#' Internal Function: Refit WQS Model for Inference
#'
#' Fits a downstream inference model using fixed WQS indices.
#'
#' @param X_quantile Quantile-transformed exposure matrix.
#' @param y Outcome vector.
#' @param cov_matrix Covariate matrix.
#' @param pos_weights Fixed positive weights.
#' @param neg_weights Fixed negative weights.
#' @param family Character.
#' @param groups Group list.
#' @param group_inference Logical. Whether to include group-specific indices.
#' @param engine Character. One of \code{"glm"} or \code{"svyglm"}.
#' @param survey_design Optional pre-constructed survey design object.
#'
#' @return A list with refit results and extracted inference summaries.
#'
#' @keywords internal
refit_model <- function(X_quantile, y, cov_matrix,
                        pos_weights, neg_weights, family, groups,
                        group_inference = TRUE,
                        engine = c("glm", "svyglm"),
                        survey_design = NULL,
                        obs_weights = NULL,
                        excluded_directions = list()) {
  engine <- match.arg(engine)
  var_names <- names(pos_weights)
  n_obs <- length(y)
  
  if (!is.null(groups) && group_inference) {
    original_grp_names <- names(groups)
    safe_grp_names <- make.names(original_grp_names, unique = TRUE)
    grp_name_to_safe <- stats::setNames(safe_grp_names, original_grp_names)
    
    group_indices <- data.frame(matrix(nrow = n_obs, ncol = 0))
    group_names_pos <- character(0)
    group_names_neg <- character(0)
    
    for (grp_name in original_grp_names) {
      grp_vars <- groups[[grp_name]]
      grp_idx <- which(var_names %in% grp_vars)
      if (length(grp_idx) == 0) next

      grp_pos_w <- pos_weights[grp_idx]
      grp_neg_w <- neg_weights[grp_idx]
      if (sum(grp_pos_w) > 0) grp_pos_w <- grp_pos_w / sum(grp_pos_w)
      if (sum(grp_neg_w) > 0) grp_neg_w <- grp_neg_w / sum(grp_neg_w)

      X_grp <- X_quantile[, grp_idx, drop = FALSE]
      safe_name <- grp_name_to_safe[[grp_name]]

      pos_excluded <- identical(excluded_directions[[grp_name]], "positive")
      neg_excluded <- identical(excluded_directions[[grp_name]], "negative")

      if (!pos_excluded) {
        pos_col_name <- paste0(safe_name, "_pos")
        group_indices[[pos_col_name]] <- as.numeric(X_grp %*% grp_pos_w)
        group_names_pos <- c(group_names_pos, pos_col_name)
      }

      if (!neg_excluded) {
        neg_col_name <- paste0(safe_name, "_neg")
        group_indices[[neg_col_name]] <- as.numeric(X_grp %*% grp_neg_w)
        group_names_neg <- c(group_names_neg, neg_col_name)
      }
    }
    
    wqs_indices <- group_indices
    wqs_terms <- paste(c(group_names_pos, group_names_neg), collapse = " + ")
  } else {
    wqs_indices <- create_wqs_index(X_quantile, pos_weights, neg_weights)
    wqs_terms <- "wqs_pos + wqs_neg"
  }
  
  glm_data <- wqs_indices
  glm_data$y <- y
  
  if (!is.null(cov_matrix)) {
    cov_df <- as.data.frame(cov_matrix)
    safe_cov_names <- make.names(names(cov_df), unique = TRUE)
    names(cov_df) <- safe_cov_names
    glm_data <- cbind(glm_data, cov_df)
    cov_terms <- paste(safe_cov_names, collapse = " + ")
    formula_str <- paste("y ~", wqs_terms, "+", cov_terms)
  } else {
    formula_str <- paste("y ~", wqs_terms)
  }
  
  fit_formula <- stats::as.formula(formula_str)
  glm_family <- .glm_family_from_string(family, engine = engine)
  
  if (engine == "glm") {
    if (!is.null(obs_weights)) {
      refit_fit <- stats::glm(fit_formula, data = glm_data, family = glm_family, weights = obs_weights)
    } else {
      refit_fit <- stats::glm(fit_formula, data = glm_data, family = glm_family)
    }
  } else {
    if (!requireNamespace("survey", quietly = TRUE)) {
      stop(
        "Package 'survey' is required for `refit_engine = \"svyglm\"`. ",
        "Install it with install.packages('survey').",
        call. = FALSE
      )
    }
    if (is.null(survey_design)) {
      stop("`survey_design` must be provided when `engine = \"svyglm\"`.", call. = FALSE)
    }
    if (!identical(nrow(survey_design$variables), nrow(glm_data))) {
      stop(
        "`survey_design` must have the same number of observations as the refit data. ",
        "For NHANES analyses, build `survey_design` from the exact same final analysis dataset ",
        "used for `sglwqs()` after all filtering and row ordering are finalized.",
        call. = FALSE
      )
    }
    design_row_names <- rownames(survey_design$variables)
    data_row_names <- rownames(glm_data)
    if (!is.null(design_row_names) && !is.null(data_row_names)) {
      if (!identical(design_row_names, data_row_names)) {
        stop(
          "Row names of `survey_design$variables` do not match the refit data. ",
          "For NHANES analyses, `survey_design` must be built from the exact same final analysis dataset ",
          "in the same row order as the data passed to `sglwqs()`.",
          call. = FALSE
        )
      }
    } else {
      warning(
        "`survey_design` row names could not be checked. For NHANES analyses, ensure that ",
        "`survey_design` was built from the exact same final analysis dataset in the same row order ",
        "used by `sglwqs()`.",
        call. = FALSE
      )
    }
    
    design_refit <- survey_design
    design_refit$variables <- as.data.frame(design_refit$variables)
    for (nm in names(glm_data)) {
      design_refit$variables[[nm]] <- glm_data[[nm]]
    }
    refit_fit <- survey::svyglm(fit_formula, design = design_refit, family = glm_family)
  }
  
  refit_summary <- summary(refit_fit)
  coef_table <- refit_summary$coefficients
  p_col <- grep("^Pr\\(>\\|[tz]\\|\\)$", colnames(coef_table), value = TRUE)
  if (length(p_col) == 0) {
    p_col <- colnames(coef_table)[ncol(coef_table)]
    warning("Could not identify p-value column; using last column: ", p_col)
  } else {
    p_col <- p_col[1]
  }
  
  group_results <- list()
  if (!is.null(groups) && group_inference) {
    for (grp_name in original_grp_names) {
      safe_name <- grp_name_to_safe[[grp_name]]
      pos_term <- paste0(safe_name, "_pos")
      neg_term <- paste0(safe_name, "_neg")
      pos_excl <- identical(excluded_directions[[grp_name]], "positive")
      neg_excl <- identical(excluded_directions[[grp_name]], "negative")
      pos_exists <- !pos_excl && pos_term %in% rownames(coef_table)
      neg_exists <- !neg_excl && neg_term %in% rownames(coef_table)

      # Always insert — dropped/excluded directions get NA (not omitted)
      group_results[[grp_name]] <- list(
        pos_estimate = if (pos_exists) coef_table[pos_term, "Estimate"] else NA_real_,
        pos_se       = if (pos_exists) coef_table[pos_term, "Std. Error"] else NA_real_,
        pos_pvalue   = if (pos_exists) coef_table[pos_term, p_col] else NA_real_,
        neg_estimate = if (neg_exists) coef_table[neg_term, "Estimate"] else NA_real_,
        neg_se       = if (neg_exists) coef_table[neg_term, "Std. Error"] else NA_real_,
        neg_pvalue   = if (neg_exists) coef_table[neg_term, p_col] else NA_real_,
        pos_excluded = pos_excl,
        neg_excluded = neg_excl
      )
    }
  }
  
  if (!is.null(groups) && group_inference) {
    wqs_pos_estimate <- NA_real_
    wqs_pos_se <- NA_real_
    wqs_pos_pvalue <- NA_real_
    wqs_neg_estimate <- NA_real_
    wqs_neg_se <- NA_real_
    wqs_neg_pvalue <- NA_real_
  } else {
    wqs_pos_exists <- "wqs_pos" %in% rownames(coef_table)
    wqs_neg_exists <- "wqs_neg" %in% rownames(coef_table)
    wqs_pos_estimate <- if (wqs_pos_exists) coef_table["wqs_pos", "Estimate"] else NA_real_
    wqs_pos_se <- if (wqs_pos_exists) coef_table["wqs_pos", "Std. Error"] else NA_real_
    wqs_pos_pvalue <- if (wqs_pos_exists) coef_table["wqs_pos", p_col] else NA_real_
    wqs_neg_estimate <- if (wqs_neg_exists) coef_table["wqs_neg", "Estimate"] else NA_real_
    wqs_neg_se <- if (wqs_neg_exists) coef_table["wqs_neg", "Std. Error"] else NA_real_
    wqs_neg_pvalue <- if (wqs_neg_exists) coef_table["wqs_neg", p_col] else NA_real_
  }
  
  list(
    refit_fit = refit_fit,
    refit_summary = refit_summary,
    coef_table = coef_table,
    wqs_pos_estimate = wqs_pos_estimate,
    wqs_pos_se = wqs_pos_se,
    wqs_pos_pvalue = wqs_pos_pvalue,
    wqs_neg_estimate = wqs_neg_estimate,
    wqs_neg_se = wqs_neg_se,
    wqs_neg_pvalue = wqs_neg_pvalue,
    group_results = group_results,
    group_inference = group_inference && !is.null(groups),
    wqs_indices = wqs_indices,
    p_col = p_col,
    formula = formula_str,
    engine = engine,
    curve = "linear",
    n_obs = n_obs
  )
}


#' Internal Function: Validation Step for Inference
#'
#' Performs validation step to obtain p-values using fixed weights.
#'
#' @param X_quantile_val Quantile-transformed validation exposure matrix.
#' @param y_val Validation outcome vector.
#' @param cov_matrix_val Validation covariate matrix.
#' @param pos_weights Fixed positive weights.
#' @param neg_weights Fixed negative weights.
#' @param family Character.
#' @param groups Group list.
#' @param group_inference Logical. Whether to include group-specific indices in GLM.
#'
#' @return A list with validation GLM results.
#'
#' @keywords internal
validation_glm <- function(X_quantile_val, y_val, cov_matrix_val,
                           pos_weights, neg_weights, family, groups,
                           group_inference = TRUE,
                           obs_weights = NULL,
                           excluded_directions = list()) {
  val_result <- refit_model(
    X_quantile = X_quantile_val,
    y = y_val,
    cov_matrix = cov_matrix_val,
    pos_weights = pos_weights,
    neg_weights = neg_weights,
    family = family,
    groups = groups,
    group_inference = group_inference,
    engine = "glm",
    obs_weights = obs_weights,
    excluded_directions = excluded_directions
  )
  
  val_result$glm_fit <- val_result$refit_fit
  val_result$glm_summary <- val_result$refit_summary
  val_result$n_train <- NA
  val_result$n_val <- length(y_val)
  val_result
}


#' @keywords internal
.validate_obs_weights <- function(obs_weights, n) {
  if (is.null(obs_weights)) {
    return(NULL)
  }
  if (!is.numeric(obs_weights)) {
    stop("`obs_weights` must be a numeric vector.", call. = FALSE)
  }
  if (length(obs_weights) != n) {
    stop("`obs_weights` must have length equal to the number of observations.", call. = FALSE)
  }
  if (any(!is.finite(obs_weights))) {
    stop("`obs_weights` must contain only finite values.", call. = FALSE)
  }
  if (any(obs_weights < 0)) {
    stop("`obs_weights` must be non-negative.", call. = FALSE)
  }
  if (sum(obs_weights) <= 0) {
    stop("`obs_weights` must have positive total weight.", call. = FALSE)
  }
  obs_weights
}


#' @keywords internal
.coerce_sparsegl_family <- function(family, obs_weights = NULL) {
  if (is.null(obs_weights)) {
    return(family)
  }
  .glm_family_from_string(family)
}


#' @keywords internal
.glm_family_from_string <- function(family, engine = "glm") {
  if (inherits(family, "family")) {
    return(family)
  }
  if (identical(family, "gaussian")) {
    return(stats::gaussian())
  }
  if (identical(family, "binomial")) {
    if (identical(engine, "svyglm")) {
      return(stats::quasibinomial())
    }
    return(stats::binomial())
  }
  stop("Unsupported family: ", family, call. = FALSE)
}
