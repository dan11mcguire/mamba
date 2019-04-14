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
double  llbR1_i(NumericVector betaij_i,NumericVector sij2_i, double tau2inv) {
  int k_i = betaij_i.size() ;
  static const double pi = 3.14159265358979323846 ;
 return -.5*sum(log(sij2_i)) +.5*log(tau2inv) -.5*log(sum(1/sij2_i) + tau2inv) -
   .5*(sum(pow(betaij_i,2)/sij2_i) - 
  pow(sum(betaij_i/sij2_i),2)/(sum(1/sij2_i) + tau2inv)) ;

}

// [[Rcpp::export]]
double logsumexp(NumericVector x) {
  return(log(sum(exp(x - max(x)))) + max(x));
}

// [[Rcpp::export]]
NumericVector llb0_ab(double betaij, double sij2, double f, double alpha) {
  NumericVector llb0_inds(2L);
  llb0_inds[0]=log(f) + R::dnorm(betaij, 0.0, sqrt(sij2), true) ;
  llb0_inds[1]=log(1-f) + R::dnorm(betaij, 0.0, sqrt(alpha*sij2), true);
  // cohortll0(j) = logsumexp(llb0_inds) ;
  
  return(llb0_inds) ;
} 



// [[Rcpp::export]]
double llbR0_i(NumericVector betaij_i, NumericVector sij2_i,
               double f, double alpha) {
  int k_i=betaij_i.size();
NumericVector cohortll0(k_i);
for(int j=0; j < k_i; j++) {
  NumericVector llb0_inds(2L);
  llb0_inds[0]=log(f) + R::dnorm(betaij_i[j], 0.0, sqrt(sij2_i[j]), true) ;
  llb0_inds[1]=log(1-f) + R::dnorm(betaij_i[j], 0.0, sqrt(alpha*sij2_i[j]), true);
   cohortll0[j] = logsumexp(llb0_inds) ;
}
 return sum(cohortll0) ;
}

//[[Rcpp::export]]
NumericVector deltis(NumericVector betaij_i, NumericVector sij2_i,
                     double f, double alpha) {
  int k_i=betaij_i.size();
  NumericVector cohortll0(k_i) ;
  for(int j=0; j<k_i; j++) {
    if(!Rcpp::NumericVector::is_na(betaij_i[j])){
      NumericVector llb0_inds(2L);
      llb0_inds[0]=log(f) + R::dnorm(betaij_i[j], 0.0, sqrt(sij2_i[j]), true) ;
      llb0_inds[1]=log(1-f) + R::dnorm(betaij_i[j], 0.0, sqrt(alpha*sij2_i[j]), true);
      cohortll0[j] = exp(llb0_inds[0]) / exp(logsumexp(llb0_inds));
      if(Rcpp::NumericVector::is_na(cohortll0[j])){
        double m=which_max(llb0_inds) ; 
        cohortll0[j]= m ;
      }
    }  else {
      cohortll0[j]=NA_REAL ;
  }
 }
 return cohortll0  ;
}






/*** R
library(parallel)
library(matrixStats)
set.seed(100)
beta<-lapply(1:(6*10^4), function(i){rnorm(100)})
s2<-lapply(1:(6*10^4), function(i){runif(100)})
alpha=12; f=0.8;
parcores=1
mis.inds<-mclapply(1:length(beta), function(i){
  which(!is.na(beta[[i]]))
}, mc.cores = parcores)

system.time({
 deltaij.rcpp<- mclapply(1:100, function(i){
    deltis(betaij_i=beta[[i]][mis.inds[[i]]], sij2_i=s2[[i]][mis.inds[[i]]], alpha=alpha, f=f)
  })
})

system.time({
deltaij<- mclapply(1:100, function(i){
  delts<-sapply(1:length(mis.inds[[i]]), function(j){
    delt<- exp(logSumExp(log(f)  + dnorm(beta[[i]][j], 0,sqrt(s2[[i]][j]), log=TRUE))) / 
      exp(logSumExp(c(log(f) + dnorm(beta[[i]][j], 0,sqrt(s2[[i]][j]), log=TRUE),
                      log(1-f) + dnorm(beta[[i]][j], 0,sqrt(alpha*s2[[i]][j]), log=TRUE))))
    if(is.na(delt)){
      m<-which.max(c(log(1-f) + dnorm(beta[[i]][j], 0,sqrt(alpha*s2[[i]][j]), log=TRUE),
                     log(f) + dnorm(beta[[i]][j], 0,sqrt(s2[[i]], log=TRUE))))                      
      delt<-(m-1)
    }
    delt
  })
  coli<-rep(NA, length(beta[[1]]))
  coli[mis.inds[[i]]]<-delts

  coli
}, mc.cores = 1)
})
sapply(1:100, function(i){
  all.equal(deltaij[[i]], deltaij.rcpp[[i]])
})

system.time({
 vec.rcpp<-unlist( mclapply(1:length(beta), function(i){
    llbR0_i(betaij_i = beta[[i]], sij2_i = s2[[i]], alpha = 12, f=.8)
  }, mc.cores = 1))#, mc.cores = 1)
})

system.time({
  vec<-unlist(mclapply(1:length(beta), function(i){
    f<-0.8; alpha<-12;
    sum(sapply(1:length(beta[[i]]), function(j){
      logSumExp(c(log(f) + dnorm(beta[[i]][j], 0, sqrt(s2[[i]][j]), log=TRUE),
                  log(1-f) + dnorm(beta[[i]][j], 0, sqrt(alpha*s2[[i]][j]), log=TRUE)))
    }))
  }, mc.cores = 1))#, mc.cores = 1)
})
all.equal(vec.rcpp, vec)
system.time({
  vec1.rcpp<-unlist(mclapply(1:length(beta), function(i){
    llbR1_i(betaij_i = beta[[i]], sij2_i=s2[[i]], tau2inv = 1/2e-5)
  }, mc.cores=1))
})
system.time({
  vec1<-unlist(mclapply(1:length(beta), function(i){
    k_i<-length(beta[[i]])
 val<- -.5*sum(log(s2[[i]])) +.5*log(1/2e-5) -.5*log(sum(1/s2[[i]]) + 1/2e-5) -
    .5*(sum(beta[[i]]^2/s2[[i]]) - 
          sum(beta[[i]]/s2[[i]])^2/(sum(1/s2[[i]]) + 1/2e-5)) -
    (k_i/2)*log(2*pi)
 return(val)
  }))
})
all.equal(vec1, vec1.rcpp)


# # llbR1_i(rnorm(10), runif(10), 1/0.002)
# # beta<-rnorm(100)
# # s2<-runif(100)
# # 
# # 
# # 
# # logsumexp(dnorm(rnorm(20), log=TRUE))
# # llb0_ab(2, 1, .2, 10)
# # f<-0.8; alpha<-12;
# # sum(sapply(1:length(beta), function(j){
# #   logSumExp(c(log(f) + dnorm(beta[j], 0, sqrt(s2[j]), log=TRUE),
# #               log(1-f) + dnorm(beta[j], 0, sqrt(alpha*s2[j]), log=TRUE)))
# # }))
# 
#   
# llbR0_i(betaij_i = beta, sij2_i = s2, alpha = 12, f=.8)
*/
  
  
  
  
