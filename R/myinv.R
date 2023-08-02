#' myinv
#'
#' This function computes the precision matrix, when the covariance matrix is singular.
#'
#' @param X is the covriance matrix
#' @param eps is the precision when calculating the precision matrix.
#'
#' @return the precision matrix
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats cov
#' @import stats
#'
#' @examples
#' set.seed(111)
#' blocksize = 5; rho = 0.9
#' cov_block = matrix(0, blocksize ,blocksize)
#' for ( i in 1:blocksize){
#' for (j in 1:blocksize){
#' if(i==j){
#' cov_block[i,j] =  1
#' }else if( abs(i-j)>0 ){
#' cov_block[i,j] =  rho^ (abs(i-j))
#' }
#' }
#' }
#' Y = MASS::mvrnorm( n = 50,  mu = rep(0, blocksize ), Sigma = cov_block )
#' COV = cov(Y)
#' myinv(COV)
#'
myinv <-
  function(X, eps=1e-12){
    eig.X <- eigen(X, symmetric=TRUE)
    P <- eig.X[[2]]
    lambda <- eig.X[[1]]
    ind <- lambda > eps
    lambda[ind] <- 1/lambda[ind]
    lambda[!ind] <- 0
    ans <- P%*%diag(lambda, nrow=length(lambda))%*%t(P)
    return(ans)
  }
