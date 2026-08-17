#' Example Dataset for sglwqs
#'
#' A small reproducible dataset for package examples and tests.
#' The dataset mimics correlated chemical exposures with both
#' positive and negative mixture effects.
#'
#' @format A data frame with 1000 rows and 12 variables:
#' \describe{
#'   \item{\code{age}}{Numeric age covariate.}
#'   \item{\code{sex}}{Binary factor covariate with levels \code{"F"} and \code{"M"}.}
#'   \item{\code{metal1}}{Correlated metal exposure concentration.}
#'   \item{\code{metal2}}{Correlated metal exposure concentration.}
#'   \item{\code{metal3}}{Correlated metal exposure concentration.}
#'   \item{\code{metal4}}{Correlated metal exposure concentration.}
#'   \item{\code{pesticide1}}{Correlated pesticide exposure concentration.}
#'   \item{\code{pesticide2}}{Correlated pesticide exposure concentration.}
#'   \item{\code{pesticide3}}{Correlated pesticide exposure concentration.}
#'   \item{\code{pesticide4}}{Correlated pesticide exposure concentration.}
#'   \item{\code{outcome_cont}}{Continuous outcome with known mixture signal.}
#'   \item{\code{outcome_bin}}{Binary outcome derived from a logistic model.}
#' }
#'
#' @details
#' Recommended exposure groups:
#' \preformatted{
#' groups <- list(
#'   Metals = paste0("metal", 1:4),
#'   Pesticides = paste0("pesticide", 1:4)
#' )
#' }
#'
#' @examples
#' data("sglwqs_example")
#' str(sglwqs_example)
#'
#' exp_vars <- c(paste0("metal", 1:4), paste0("pesticide", 1:4))
#' fit <- sglwqs(
#'   X = sglwqs_example[, exp_vars],
#'   y = sglwqs_example$outcome_cont,
#'   covariates = sglwqs_example[c("age", "sex")],
#'   nfolds = 3,
#'   seed = 123,
#'   verbose = FALSE
#' )
#' print(fit)
#'
#' @name sglwqs_example
NULL
