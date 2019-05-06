#include <TMB.hpp>
// using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

template <class Type>
Type objective_function <Type>:: operator () ()
{
  // data ( input from R)
  DATA_VECTOR (b2s2 );
  DATA_VECTOR (bs22 );
  DATA_VECTOR (os22 );
  DATA_VECTOR (gammai );
  
  // parameters ( input from R)
  PARAMETER ( logTau2 );
  Type tau2=exp(logTau2);
  Type nll=0; 
  for(int i=0; i< b2s2.size(); i++) {
    nll -= gammai(i) * (log(1/tau2) - log(os22(i) + 1/tau2) - (b2s2(i)
                              - bs22(i)/(os22(i) + 1/tau2))) ;
  }
  return nll;
}