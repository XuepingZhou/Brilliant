// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Brilliant algorithm
//' 
//' @param  X Predictor matrix.N by P matrix
//' @param  Y Outcome matrix. N by Q matrix.
//' @param  Omega Omega of Y. Users can use BrilliantOmega to calculate. 
//' @param  grplistP Grouping information of X. A matrix contains two columns. The number of rows corresponds to the number of groups in Y. For each row, the number in the first column is the position of the first element of the group, and the number in the second column is the position of the last element of the group. The position starts from 0.
//' @param  grplistQ Grouping information of Y. A matrix contains two columns. The number of rows corresponds to the number of groups in Y. For each row, the number in the first column is the position of the first element of the group, and the number in the second column is the position of the last element of the group. The position starts from 0.
//' @param  lambda Lasso penalty
//' @param  lambdaG Group lasso penalty
//' @param  grBnomThred 
//' @param  convThreshold. Infinity norm convergence criterion
//' @param  maxiter. Maxium iteration number.
//' @param  MaxNG = 1 . Groups are not overlapping
arma::mat Position(arma::mat& grplistP , arma::mat& grplistQ, 
                   arma::uword  p, arma::uword q,   
                   arma::uword  G, arma::uword R,  arma::uword MaxNG){
  
  arma::uword gFrom, gTo, rFrom, rTo , g, r; 
  arma::uword indexg = 0 ; 
  arma::mat resmat(MaxNG, 2) ;  /// col 1st is g, 2nd is r 
  
  for(g=0; g<G; g++  ){  for(r=0; r<R; r++ ){
    gFrom = grplistP(g,0) ;  gTo = grplistP(g,1);
    rFrom = grplistQ(r,0) ;  rTo = grplistQ(r,1); 
    
    if (p >= gFrom && p <= gTo  && q >= rFrom && q <= rTo) {
      resmat(indexg, 0) = g ; 
      resmat(indexg, 1) = r ; 
      indexg++  ; 
    }
  }}
  resmat = resmat.head_rows(indexg);
  return  resmat ; 
}


double SoftShrink(double rho, double lambda){
  if (rho < -1 * lambda){
    return(rho + lambda);
  }else if(rho > lambda){
    return(rho - lambda);
  }else {
    return(0);
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::List MSGExt(arma::mat X, arma::mat Y, arma::mat Omega, 
                  arma::mat grplistP, arma::mat grplistQ , 
                  double lambda, double lambdaG,
                  double grBnomThred = 0.1,  double convThreshold =1e-2, double maxiter = 1000 , 
                  arma::uword MaxNG = 1  
){
  
  ///////////////////////////////////////////////////////////////////
  ///////////// (1) Dimensions 
  //////////////////////////////////////////////////////////////////
  int N = X.n_rows ;          /// const int N = X.nrow() ; 
  arma::uword P = X.n_cols ;  /// const int P = X.ncol() ; 
  arma::uword Q = Y.n_cols ; 
  arma::uword G = grplistP.n_rows ; 
  arma::uword R = grplistQ.n_rows ; 
  
  if(X.n_rows != Y.n_rows){ stop("Error: Check X and Y dimensions"); }
  if(Omega.n_cols != Y.n_cols){stop("Error: Check Y and Omega dimensions"); }
  
  arma::uword p, q, g, r, i; 
  arma::uword gFrom, gTo, rFrom, rTo ; 
  arma::vec grpp, grpq ; 
  arma::mat B(P, Q); 
  arma::mat BThred(P, Q); 
  
  /// double lambda; double lambdaG ; 
  int converge = 0 ; 
  int iteration = 0; 
  
  
  ///////////////////////////////////////////////////////////////////////////
  ///////////// (2) Initiation
  ////////////////////////////////////////////////////////////////////////////
  /// Initialize B matrix  
  
  /// Calculate (||x_j||_2)^2 , the column j of X   )
  arma::vec XnormSquare = arma::vec(P);
  for ( i = 0; i < P; i++) {
    XnormSquare(i) = pow(arma::norm(X.col(i), "fro") , 2) ;
  }
  
  ///  Initialize S  -- S1 = X[,j]T* (Y-XB(-ij))* Omega[.k]
  arma::mat A = X.t() * Y * Omega ; 
  arma::mat C = X.t() * X * B * Omega ;
  arma::mat omegaDiag(Q,Q); omegaDiag.diag() = diagvec(Omega); 
  arma::mat XnormDiag(P,P); XnormDiag.diag() = XnormSquare; 
  arma::mat D  = XnormDiag  * B  * omegaDiag;   
  arma::mat S  = A - C + D ; 
  
  
  ///  Initialize softShrk, apply(S, c(1, 2), function(x) {SoftShrink(x/N, lambda )})
  arma::mat softShrk(P, Q) ;     
  for(p=0;p<P;p++){
    for(q=0;q<Q;q++){
      softShrk(p,q) =  SoftShrink(S(p,q)/N, lambda ) ; 
    }
  }
  
  /// Xnormsquare_Omega
  arma::mat Xnormsquare_Omega(P,Q) ;      ////# * Element-wise multiplication
  arma::vec tmpvec = diagvec(Omega) ; 
  for( p=0;p<P;p++){ for( q=0;q<Q;q++){
    Xnormsquare_Omega(p,q) =  XnormSquare[p] *  tmpvec[q] ;
  }}
  
  ///For each group -- calculate softShrkL2Norm, grWeight  and  L2 norm (grp_Norm) 
  arma::mat softShrkL2Norm(G, R) ; 
  arma::mat grWeight(G, R) ; 
  arma::mat grpNorm(G, R) ; 
  
  for(g=0; g<G; g++  ){ for(r=0; r<R; r++ ){
    gFrom = grplistP(g,0) ;   gTo = grplistP(g,1);
    rFrom = grplistQ(r,0) ;   rTo = grplistQ(r,1); 
    grpp =  arma::linspace(gFrom, gTo,  (gTo- gFrom +1));
    grpq =  arma::linspace(rFrom, rTo,  (rTo- rFrom +1));
    softShrkL2Norm(g,r) =  arma::norm(softShrk(arma::span(gFrom, gTo) , arma::span(rFrom, rTo)), "fro") ; 
    grpNorm(g,r)   =  arma::norm(B(arma::span(gFrom, gTo) , arma::span(rFrom, rTo)), "fro") ; 
    grWeight(g,r)  = grpp.n_elem * grpq.n_elem ; /// sqrt(grpp.n_elem * grpq.n_elem ) ; 
  }}
  
  
  ///////////////////////////////////////////////////////////////////////////////
  ///////////// (3) begin update
  //////////////////////////////////////////////////////////////////////////////
  Rcpp::Rcout << " Iteration begins" << std::endl;
  converge = 0, iteration = 0; int b0flag ; 
  
  double S1, tempsoftShrk, zeroOtherNormSum, nonzeroOtherNormSum, diff1 ;
  arma::uword indexg ; 
  arma::mat Blast, M; 
  arma::mat tmpmat ;
  
  double tmpscalar ; 
  
  while(converge == 0 ){   
    iteration++ ; 
    if( (iteration % 30) == 0 ) { Rcpp::Rcout << " Iteration=" << iteration  << std::endl;}
    /// Update B 
    Blast = B ; 
    for(p=0;p<P; p++){  for(q=0;q<Q; q++){
      M = B ;  M(p,q) = 0 ; 
      S1 = as_scalar(A(p,q) - X.col(p).t() * X * M * Omega.col(q)) ;
      zeroOtherNormSum = lambda ; 
      nonzeroOtherNormSum = 0 ; 
      b0flag= 0 ; 
      
      ///  find b(p,q) belong to group  ///  # Mixed coordinate descent algorithm
      tmpmat = Position(grplistP ,  grplistQ,  p,  q, G, R,  MaxNG   ) ; 
      for (indexg=0 ; indexg < tmpmat.n_rows; indexg++) {
        g = tmpmat(indexg,0) ; r = tmpmat(indexg,1) ;
        if( softShrkL2Norm(g,r) <=  grWeight(g,r) * lambdaG   ){ // # Theorem 3.1. whole grp is 0s 
          b0flag= 1 ; 
        }else{
          tmpscalar = pow(grpNorm(g,r), 2 ) - pow(B(p,q), 2) ; 
          if( tmpscalar  <=  grBnomThred ){       // or 0 
            zeroOtherNormSum = zeroOtherNormSum + grWeight(g,r)*lambdaG ;
          }else{   
            nonzeroOtherNormSum = nonzeroOtherNormSum + grWeight(g,r) *lambdaG /grpNorm(g,r)  ;
          }
        }
      }// end all group for B(p,q) belong group 
      
      ///  Update B(p,q) 
      if (b0flag == 1) {  
        B(p,q) = 0 ; 
      }else{
        B(p,q)= SoftShrink(S1/N , zeroOtherNormSum)/( Xnormsquare_Omega(p,q) /N  + nonzeroOtherNormSum) ; 
      }
      
      
      ///  Update B(p,q) related information   ---------------------  
      if ( B(p,q) != Blast(p,q) ) {
        tempsoftShrk=SoftShrink(S1/N, zeroOtherNormSum)  ; 
        /// ## Update grp_Norm(g,r) and softShrkL2Norm(g,r)
        for (indexg=0 ; indexg < tmpmat.n_rows; indexg++) {
          g = tmpmat(indexg,0) ; r = tmpmat(indexg,1) ;
          
          /// grBnorm
          tmpscalar = pow(grpNorm(g,r), 2 ) + pow(B(p,q), 2) - pow(Blast(p,q), 2); 
          if( tmpscalar <= 0 ){ grpNorm(g,r)=0 ; 
          }else{grpNorm(g,r)= sqrt(tmpscalar) ; }
          
          /// softShrkL2Norm
          tmpscalar = pow(softShrkL2Norm(g,r) , 2 ) + pow(tempsoftShrk, 2) - pow(softShrk(p,q), 2); 
          if( tmpscalar <= 0){softShrkL2Norm(g,r) =0; 
          }else{softShrkL2Norm(g,r)=sqrt(tmpscalar); }
          
          /// softShrk(p,q)
          softShrk(p,q)  = tempsoftShrk ;
        }  
        
      } // end update b_bq related info
      
      
      
      
    } }  /// end p ,q 
    
    /*
     List result = List::create(
     Named("softShrk") = softShrk, 
     Named("grWeight") = grWeight, 
     Named("softShrkL2Norm") = softShrkL2Norm, 
     Named("S") = S,  Named("A") = A, 
     Named("B")=B)  ;
     return(result) ; 
     */
    /////////////////////////////////////////////////////////////////////////
    /// Thresholding L2 groups after all p, q
    /////////////////////////////////////////////////////////////////////////
    
    for(g=0; g<G; g++  ){  for(r=0; r<R; r++ ){
      gFrom = grplistP(g,0) ;   gTo = grplistP(g,1);
      rFrom = grplistQ(r,0) ;   rTo = grplistQ(r,1); 
      
      if( softShrkL2Norm(g,r) <=  grWeight(g,r) * lambdaG){    
        //  Rcpp::Rcout << "G" << g << "R" << r << " is a 0 group" << std::endl;
        // if (grpNorm(g,r) != 0 ){
        B(arma::span(gFrom, gTo) , arma::span(rFrom, rTo)).zeros() ; 
        grpNorm(g,r) = 0 ; 
        
        for(p = gFrom ;  p <= gTo ; p++){ for( q = rFrom ; q <= rTo ; q++){
          tmpscalar = arma::as_scalar(A(p,q) - X.col(p).t() * X * B  *  Omega.col(q)) ; 
          tempsoftShrk= SoftShrink(tmpscalar /N, lambda) ; /// softShrk(p,q) =  SoftShrink(S(p,q)/N, lambda ) ;
          tmpscalar = pow(softShrkL2Norm(g ,r),2) + pow(tempsoftShrk,2) - pow(softShrk(p, q),2); 
          if ( tmpscalar > 0 ){
            softShrkL2Norm(g,r)=  sqrt( tmpscalar) ; 
          }else{
            softShrkL2Norm(g,r) = 0 ; 
          }
          softShrk(p,q)  = tempsoftShrk  ;
        }}// end p, q 
        
        //  } // # end  grp_Norm[g,r) != 0
      } // # end  softShrkL2Norm(g,r) <= grWeight(g,r) * lambdaG
      
    } } // end g and r 
    
    /////////////////////////////////////////////////////////////
    //// Converge or not 
    /////////////////////////////////////////////////////////////
    tmpmat = B - Blast; 
    diff1 = arma::as_scalar( norm(tmpmat.as_col(), "inf")) ; 
    
    
    if( (iteration % 50) == 0 ) { Rcpp::Rcout  <<"Max difference is " << diff1 << std::endl; }
    if(  diff1  <=  convThreshold ) { 
      Rcpp::Rcout << "Converged. Max difference is " << diff1 << std::endl;
      converge = 1 ;  
    } 
    
    if(  diff1  <= convThreshold ||  iteration >= maxiter || diff1 > 1e5) {   
      BThred  = B ; 
      /// BThred.elem( find(abs(BThred) <=  betaThred ) ).zeros();
      
      for(g=0; g<G; g++  ){   for(r=0; r<R; r++ ){
        gFrom = grplistP(g,0) ;   gTo = grplistP(g,1);
        rFrom = grplistQ(r,0) ;   rTo = grplistQ(r,1); 
        grpNorm(g,r)   =  arma::norm(BThred(arma::span(gFrom, gTo) , arma::span(rFrom, rTo)), "fro") ; 
      } }
      
      List result = List::create(
        Named("beta") = BThred, Named("betaori") = B, Named("grp_Norm") = grpNorm, 
              Named("converge") = converge, Named("iteration") = iteration, 
              Named("lambda") = lambda , Named("lambdaG") = lambdaG)  ; 
      return result ; 
    }
    
    
    // result.attr("class") = "data.frame";
    // result.attr("row.names") = seq(1, N);
  } /// end of while 
  return(0); 
}  /// End of function 
