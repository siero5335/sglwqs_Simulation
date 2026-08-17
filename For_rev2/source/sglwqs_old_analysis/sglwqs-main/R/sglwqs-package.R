#' sglwqs: Sparse Group Lasso Weighted Quantile Sum Regression
#'
#' @description
#' The sglwqs package implements Sparse Group Lasso Weighted Quantile Sum (SGL-WQS)
#' regression for analyzing the effects of chemical mixtures on health outcomes.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{sglwqs}}}{Main function for fitting SGL-WQS models}
#'   \item{\code{\link{sglwqs_mice}}}{Apply SGL-WQS to multiply imputed datasets}
#' }
#'
#' @section Unified S3 API:
#' These functions work for both sglwqs and sglwqs_mids objects:
#' \describe{
#'   \item{\code{\link{extract_weights}}}{Extract weights as data frame}
#'   \item{\code{\link{plot_weights}}}{Butterfly chart of weights}
#'   \item{\code{\link{plot_selection_frequency}}}{Selection-frequency visualization}
#'   \item{\code{\link{plot_bootstrap_ci}}}{Bootstrap coefficient CI plot}
#'   \item{\code{\link{plot_validation_results}}}{Forest plot of mixture effects}
#'   \item{\code{\link{plot_cv}}}{Cross-validation curve visualization}
#'   \item{\code{\link{calibrate}}}{Calibration assessment for binomial models}
#'   \item{\code{\link{inference_table}}}{Create publication-ready inference table}
#'   \item{\code{\link{plot_combined_results}}}{Publication-ready combined figure}
#' }
#'
#' @section Standard Model Interface:
#' \describe{
#'   \item{\code{\link{predict.sglwqs}}}{Predict WQS indices and responses on new data}
#'   \item{\code{\link{coef.sglwqs}}}{Extract coefficients or normalized weights}
#'   \item{\code{\link{confint.sglwqs}}}{Confidence intervals for validation and bootstrap results}
#'   \item{\code{\link{plot_cv}}}{Cross-validation curve visualization for lambda selection}
#' }
#'
#' @section Interoperability:
#' \describe{
#'   \item{\code{\link{tidy.sglwqs}}, \code{\link{glance.sglwqs}}, \code{\link{augment.sglwqs}}}{broom-compatible methods}
#'   \item{\code{\link{validation_metrics}}}{Predictive performance metrics}
#'   \item{\code{\link{calibrate}}}{Calibration assessment for binomial models}
#' }
#'
#' @section Additional Functions:
#' \describe{
#'   \item{\code{\link{summary_bootstrap}}}{Detailed bootstrap summary}
#'   \item{\code{\link{summary_validation}}}{Detailed validation summary}
#'   \item{\code{\link{plot_selection_frequency}}}{Bootstrap selection frequency plot}
#'   \item{\code{\link{plot_validation_results}}}{Forest plot of mixture effects}
#' }
#'
#' @importFrom ggplot2 autoplot
#' @importFrom sparsegl cv.sparsegl sparsegl
#' @importFrom stats quantile coef glm gaussian binomial as.formula model.matrix sd qnorm var pt qt pnorm plogis aggregate predict
#' @importFrom rlang .data
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate arrange desc
#' @importFrom utils txtProgressBar setTxtProgressBar head object.size
#' @importFrom generics tidy glance augment
#' @name sglwqs-package
#' @aliases sglwqs-pkg
#' @keywords internal
"_PACKAGE"
