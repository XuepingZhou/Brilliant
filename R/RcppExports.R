# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Brilliant algorithm
#' 
#' @param  X Predictor matrix.N by P matrix
#' @param  Y Outcome matrix. N by Q matrix.
#' @param  Omega Omega of Y. Users can use BrilliantOmega to calculate. 
#' @param  grplistP Grouping information of X. A matrix contains two columns. The number of rows corresponds to the number of groups in Y. For each row, the number in the first column is the position of the first element of the group, and the number in the second column is the position of the last element of the group. The position starts from 0.
#' @param  grplistQ Grouping information of Y. A matrix contains two columns. The number of rows corresponds to the number of groups in Y. For each row, the number in the first column is the position of the first element of the group, and the number in the second column is the position of the last element of the group. The position starts from 0.
#' @param  lambda Lasso penalty
#' @param  lambdaG Group lasso penalty
#' @param  grBnomThred 
#' @param  convThreshold. Infinity norm convergence criterion
#' @param  maxiter. Maxium iteration number.
#' @param  MaxNG = 1 . Groups are not overlapping
NULL

#' @export
MSGExt <- function(X, Y, Omega, grplistP, grplistQ, lambda, lambdaG, grBnomThred = 0.1, convThreshold = 1e-2, maxiter = 1000, MaxNG = 1L) {
    .Call(`_Brilliant_MSGExt`, X, Y, Omega, grplistP, grplistQ, lambda, lambdaG, grBnomThred, convThreshold, maxiter, MaxNG)
}

