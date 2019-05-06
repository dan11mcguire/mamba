#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp:depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double  llbR1_j(NumericVector betajk_j,NumericVector sjk2_j, double tau2inv) {
  int k_j = betajk_j.size() ;
  static const double pi = 3.14159265358979323846 ;
 return -.5*sum(log(sjk2_j)) +.5*log(tau2inv) -.5*log(sum(1/sjk2_j) + tau2inv) -
   .5*(sum(pow(betajk_j,2)/sjk2_j) - 
  pow(sum(betajk_j/sjk2_j),2)/(sum(1/sjk2_j) + tau2inv)) ;

}

// [[Rcpp::export]]
double logsumexp(NumericVector x) {
  return(log(sum(exp(x - max(x)))) + max(x));
}

// [[Rcpp::export]]
NumericVector llb0_ab(double betajk, double sjk2, double lambda, double alpha) {
  NumericVector llb0_inds(2L);
  llb0_inds[0]=log(lambda) + R::dnorm(betajk, 0.0, sqrt(sjk2), true) ;
  llb0_inds[1]=log(1-lambda) + R::dnorm(betajk, 0.0, sqrt(alpha*sjk2), true);
  // cohortll0(j) = logsumexp(llb0_inds) ;
  
  return(llb0_inds) ;
} 



// [[Rcpp::export]]
double llbR0_j(NumericVector betajk_j, NumericVector sjk2_j,
               double lambda, double alpha) {
  int k_j=betajk_j.size();
NumericVector cohortll0(k_j);
for(int k=0; k < k_j; k++) {
  NumericVector llb0_inds(2L);
  llb0_inds[0]=log(lambda) + R::dnorm(betajk_j[k], 0.0, sqrt(sjk2_j[k]), true) ;
  llb0_inds[1]=log(1-lambda) + R::dnorm(betajk_j[k], 0.0, sqrt(alpha*sjk2_j[k]), true);
   cohortll0[k] = logsumexp(llb0_inds) ;
}
 return sum(cohortll0) ;
}

//[[Rcpp::export]]
NumericVector deltis(NumericVector betajk_j, NumericVector sjk2_j,
                     double lambda, double alpha) {
  int k_j=betajk_j.size();
  NumericVector cohortll0(k_j) ;
  for(int k=0; k<k_j; k++) {
    if(!Rcpp::NumericVector::is_na(betajk_j[k])){
      NumericVector llb0_inds(2L);
      llb0_inds[0]=log(lambda) + R::dnorm(betajk_j[k], 0.0, sqrt(sjk2_j[k]), true) ;
      llb0_inds[1]=log(1-lambda) + R::dnorm(betajk_j[k], 0.0, sqrt(alpha*sjk2_j[k]), true);
      cohortll0[k] = exp(llb0_inds[0]) / exp(logsumexp(llb0_inds));
      if(Rcpp::NumericVector::is_na(cohortll0[k])){
        double m=which_max(llb0_inds) ; 
        cohortll0[k]= m ;
      }
    }  else {
      cohortll0[k]=NA_REAL ;
  }
 }
 return cohortll0  ;
}


  
