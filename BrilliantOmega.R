#' BrilliantOmega
#'
#' The function calculates the simplified precision matrix of Y, which will be used in the Brilliant algorithm.
#'
#' @param Y Outcome matrix.
#' @param grp.list.Q A matrix contains two columns. The number of rows corresponds to the number of groups in Y. For each row, the number in the first column is the position of the first element of the group, and the number in the second column is the position of the last element of the group. The position starts from 1.
#'
#' @return A matrix of the (simplified) precision matrix of Y
#' @export
#' @importFrom MASS mvrnorm
#' @import stats
#' @examples
#' # Simulate outcomes with two groups and each group contains five variables
#' set.seed(111)
#' blocksize = 5 ; rho = .9
#' cov_block = matrix(0, blocksize , blocksize)
#' for ( i in 1:blocksize){
#'   for (j in 1:blocksize){
#'      if(i==j){ cov_block[i,j] =  1
#'   }else if( abs(i-j)>0 ){
#'      cov_block[i,j] =  rho^ (abs(i-j))
#'   }
#'  }
#' }
#' Y1 = MASS::mvrnorm( n = 50,  mu = rep(0, blocksize ), Sigma = cov_block )
#' Y2 = MASS::mvrnorm( n = 50,  mu = rep(0, blocksize ), Sigma = cov_block )
#' Y = cbind(Y1,Y2)
#' # Group structure of Y
#' Q = 10 # number variables in Y
#' RarrStarts <- seq(1, 10, by = blocksize)
#' RarrEnds <- RarrStarts +  (blocksize -1)
#' grp.list.Q <- data.frame(RarrStarts, RarrEnds)
#' # The smplified precision matrix
#' BrilliantOmega (Y = Y, grp.list.Q=grp.list.Q)
#'
#'
BrilliantOmega <- function(Y, grp.list.Q){
  R = dim(grp.list.Q)[1]
  Q = dim(Y)[2]
  Omega_thred <- matrix(0,Q, Q)
  for (r in 1:R){
    grpr <- grp.list.Q[r,1]:grp.list.Q[r,2]
    COV.grp <- cov(Y[grpr,  grpr])
    Omega_thred[grpr,  grpr] <- myinv(COV.grp)
  }
  return(Omega_thred)
}



