#' Simulated data for Brilliant algorithm
#'
#' Simulated data for Brilliant algorithm with sample size = 400.
#' \itemize{
#'   \item X. X contains 200 variables within in 4 groups. Each group has 50 variables with covariance structure following first-order auto-regressive structure and rho = 0.5.
#'   \item Y. Y contains 100 variables within in 5 groups. Each group has 20 variables with covariance structure following compound symmetry structure and rho = 0.9.
#'   \item B. The coefficient matrix.
#' }
#' @docType data
#'
#' @usage data(BrilliantData)
#'
#' @format A list contains 3 elements.
#'
#' @keywords datasets
#'
#' @examples
#' data(BrilliantData)
#' X <- BrilliantData$X
#' Y <- BrilliantData$Y
#' B <- BrilliantData$B
"BrilliantData"
