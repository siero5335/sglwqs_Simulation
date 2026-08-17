# Tests for special character handling in variable/group names
# Covers metabolomics-style names with parentheses, colons, slashes, spaces, etc.

test_that("sanitize_for_formula handles metabolomics-style names", {
  # Typical metabolomics feature names
  names <- c(
    "PC(16:0/18:1)",       # parentheses, colons, slashes
    "LysoPC 18:0",          # space, colon
    "2-Aminobutyric acid",  # starts with number, hyphen, space
    "alpha-Tocopherol",     # hyphen
    "normal_name",          # already safe
    "C18:1n-9",             # colon, hyphen
    "n-Decane/Isodecane"    # hyphen, slash
  )
  
  result <- sglwqs:::sanitize_for_formula(names)
  
  # All safe names should be valid R identifiers
  expect_equal(result$safe, make.names(names, unique = TRUE))
  
  # Original names preserved

  expect_equal(result$original, names)
  
  # Mapping is bidirectional
  for (i in seq_along(names)) {
    expect_equal(result$to_original[[result$safe[i]]], names[i])
    expect_equal(result$to_safe[[names[i]]], result$safe[i])
  }
  
  # Was modified flag
  expect_true(result$was_modified)
  
  # Safe names work in formula
  formula_str <- paste("y ~", paste(result$safe, collapse = " + "))
  expect_no_error(as.formula(formula_str))
})

test_that("sanitize_for_formula handles already-safe names", {
  names <- c("PCB_118", "PCB_153", "lead", "mercury")
  result <- sglwqs:::sanitize_for_formula(names)
  
  expect_equal(result$safe, names)
  expect_false(result$was_modified)
})

test_that("sanitize_for_formula handles NULL/empty input", {
  result <- sglwqs:::sanitize_for_formula(NULL)
  expect_null(result$safe)
  expect_false(result$was_modified)
  
  result2 <- sglwqs:::sanitize_for_formula(character(0))
  expect_equal(length(result2$safe), 0)
  expect_false(result2$was_modified)
})

test_that("sanitize_for_formula handles duplicate names after sanitization", {
  # "PC(16:0)" and "PC(16.0)" would both become "PC.16.0." without unique = TRUE
  names <- c("PC(16:0)", "PC(16.0)")
  result <- sglwqs:::sanitize_for_formula(names)
  
  # Should be unique
  expect_equal(length(unique(result$safe)), 2)
  
  # Original names preserved correctly in mapping
  expect_equal(result$to_original[[result$safe[1]]], names[1])
  expect_equal(result$to_original[[result$safe[2]]], names[2])
})

test_that("process_covariates handles special character column names", {
  # Numeric covariates with special characters
  n <- 50
  cov_df <- data.frame(
    `Body Mass Index` = rnorm(n),
    `Age (years)` = sample(20:60, n, replace = TRUE),
    `Income/$1000` = rnorm(n),
    check.names = FALSE
  )
  
  result <- sglwqs:::process_covariates(cov_df)
  
  # Column names should be formula-safe
  expect_true(all(result$names == make.names(result$names)))
  
  # Should work in formula
  formula_str <- paste("y ~", paste(result$names, collapse = " + "))
  expect_no_error(as.formula(formula_str))
  
  # Original names preserved
  expect_equal(result$original_names, c("Body Mass Index", "Age (years)", "Income/$1000"))
})

test_that("process_covariates handles factor covariates with special chars", {
  n <- 50
  cov_df <- data.frame(
    `smoking status` = factor(sample(c("never", "former", "current"), n, replace = TRUE)),
    `race/ethnicity` = factor(sample(c("White", "Black", "Hispanic/Latino"), n, replace = TRUE)),
    check.names = FALSE
  )
  
  result <- sglwqs:::process_covariates(cov_df)
  
  # Column names should be formula-safe
  expect_true(all(result$names == make.names(result$names)))
  
  # Matrix should have correct dimensions (dummy coded)
  expect_true(ncol(result$matrix) >= 2)  # at least some dummy variables
})

test_that("sglwqs works with special character variable names (no validation)", {
  set.seed(42)
  n <- 100
  p <- 6
  
  # Create data with metabolomics-style names
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c(
    "PC(16:0/18:1)", "LysoPC 18:0", "2-Aminobutyric acid",
    "alpha-Tocopherol", "C18:1n-9", "n-Decane/Isodecane"
  )
  y <- rnorm(n)
  
  # Should not error
  fit <- sglwqs(X, y, verbose = FALSE)
  
  # Original names should be preserved in output
  expect_equal(fit$var_names, colnames(X))
  expect_equal(names(fit$pos_weights), colnames(X))
  expect_equal(names(fit$neg_weights), colnames(X))
})

test_that("sglwqs works with special character variable AND group names (with validation)", {
  set.seed(42)
  n <- 200
  p <- 8
  
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- c(
    "PC(16:0/18:1)", "LysoPC 18:0", "PE(18:0/20:4)",
    "2-Aminobutyric acid", "alpha-Tocopherol",
    "TCEP", "TCIPP", "TDCIPP"
  )
  
  beta <- c(0.5, 0.3, 0, 0, -0.4, 0, 0, 0)
  y <- X %*% beta + rnorm(n, sd = 0.5)
  
  groups <- list(
    "Phospholipids (PC/PE)" = c("PC(16:0/18:1)", "LysoPC 18:0", "PE(18:0/20:4)"),
    "Amino acids" = c("2-Aminobutyric acid", "alpha-Tocopherol"),
    "OPFRs/Plasticizers" = c("TCEP", "TCIPP", "TDCIPP")
  )
  
  # Should work with validation (formula construction is the critical test)
  fit <- sglwqs(
    X, y,
    groups = groups,
    validation = TRUE,
    train_prop = 0.6,
    seed = 123,
    verbose = FALSE
  )
  
  # Original names preserved
  expect_equal(fit$var_names, colnames(X))
  expect_equal(names(fit$groups), names(groups))
  
  # Validation should have completed successfully
  expect_true(fit$validation)
  expect_true(!is.null(fit$validation_info))
  
  # Group results should use original group names
  if (!is.null(fit$validation_info$group_results)) {
    result_grp_names <- names(fit$validation_info$group_results)
    expect_true(all(result_grp_names %in% names(groups)))
  }
  
  # Summary should work
  expect_no_error(capture.output(summary(fit)))
})

test_that("sglwqs works with special character covariates + validation", {
  set.seed(42)
  n <- 200
  
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- c("PC(16:0)", "LysoPC 18:0", "TCEP", "TDCIPP")
  
  covariates <- data.frame(
    `Body Mass Index` = rnorm(n, 25, 5),
    `Age (years)` = sample(20:60, n, replace = TRUE),
    check.names = FALSE
  )
  
  y <- X[,1] * 0.5 + covariates[,1] * 0.1 + rnorm(n, sd = 0.5)
  
  fit <- sglwqs(
    X, y,
    covariates = covariates,
    validation = TRUE,
    train_prop = 0.6,
    seed = 123,
    verbose = FALSE
  )
  
  expect_true(fit$validation)
  expect_no_error(capture.output(summary(fit)))
})

test_that("extract_weights preserves original names with special characters", {
  set.seed(42)
  n <- 100
  
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- c("PC(16:0/18:1)", "LysoPC 18:0", "TCEP", "TDCIPP")
  y <- rnorm(n)
  
  fit <- sglwqs(X, y, verbose = FALSE)
  
  weights_df <- extract_weights(fit)
  
  # variable column should contain original names
  expect_true("PC(16:0/18:1)" %in% weights_df$variable)
  expect_true("LysoPC 18:0" %in% weights_df$variable)
})

test_that("plot_weights works with special character names", {
  set.seed(42)
  n <- 100
  
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- c("PC(16:0/18:1)", "LysoPC 18:0", "TCEP", "TDCIPP")
  y <- rnorm(n)
  
  groups <- list(
    "Lipids (PL)" = c("PC(16:0/18:1)", "LysoPC 18:0"),
    "OPFRs" = c("TCEP", "TDCIPP")
  )
  
  fit <- sglwqs(X, y, groups = groups, verbose = FALSE)
  
  # Plot should not error
  expect_no_error(plot_weights(fit))
})

test_that("bootstrap + validation works with special character names", {
  set.seed(42)
  n <- 200
  
  X <- matrix(rnorm(n * 4), n, 4)
  colnames(X) <- c("PC(16:0/18:1)", "LysoPC 18:0", "C18:1n-9", "2-Amino acid")
  
  y <- X[,1] * 0.5 - X[,4] * 0.3 + rnorm(n, sd = 0.5)
  
  groups <- list(
    "Phospholipids (PC)" = c("PC(16:0/18:1)", "LysoPC 18:0"),
    "Others" = c("C18:1n-9", "2-Amino acid")
  )
  
  fit <- suppressWarnings(sglwqs(
    X, y,
    groups = groups,
    bootstrap = TRUE,
    n_boot = 10,
    validation = TRUE,
    train_prop = 0.6,
    seed = 123,
    verbose = FALSE
  ))
  
  expect_true(fit$bootstrap)
  expect_true(fit$validation)
  expect_equal(fit$var_names, colnames(X))
  
  # All output functions should work
  expect_no_error(capture.output(summary(fit)))
  expect_no_error(extract_weights(fit))
  expect_no_error(plot_weights(fit))
  expect_no_error(inference_table(fit))
})

test_that("validate_names issues informative message for special characters", {
  expect_message(
    sglwqs:::validate_names(
      c("PC(16:0)", "normal", "2-Amino"),
      list("Group (A)" = c("PC(16:0)"))
    ),
    "special characters"
  )
  
  # No message for safe names
  expect_silent(
    sglwqs:::validate_names(c("PCB_118", "lead"), NULL, verbose = TRUE)
  )
})
