make_simple_data <- function(n = 180, seed = 123) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  x4 <- rnorm(n)
  y <- 1.8 * x1 - 1.6 * x2 + rnorm(n, sd = 0.4)
  data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4, y = y)
}

make_binomial_data <- function(n = 220, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  eta <- 1.2 * x1 - 1.0 * x2 + 0.3 * x3
  p <- stats::plogis(eta)
  y <- rbinom(n, 1, p)
  data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
}

make_uncertainty_data <- function(n = 140, seed = 77) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  x3 <- rnorm(n)
  y <- 1.5 * x1 - 1.0 * x2 + 0.2 * x3 + rnorm(n, sd = 0.5)
  data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
}
