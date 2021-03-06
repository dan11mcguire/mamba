// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// llbR1_j
double llbR1_j(NumericVector betajk_j, NumericVector sjk2_j, double tau2inv);
RcppExport SEXP _mamba_llbR1_j(SEXP betajk_jSEXP, SEXP sjk2_jSEXP, SEXP tau2invSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j(betajk_jSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j(sjk2_jSEXP);
    Rcpp::traits::input_parameter< double >::type tau2inv(tau2invSEXP);
    rcpp_result_gen = Rcpp::wrap(llbR1_j(betajk_j, sjk2_j, tau2inv));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp
double logsumexp(NumericVector x);
RcppExport SEXP _mamba_logsumexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp(x));
    return rcpp_result_gen;
END_RCPP
}
// llb0_ab
NumericVector llb0_ab(double betajk, double sjk2, double lambda, double alpha);
RcppExport SEXP _mamba_llb0_ab(SEXP betajkSEXP, SEXP sjk2SEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type betajk(betajkSEXP);
    Rcpp::traits::input_parameter< double >::type sjk2(sjk2SEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(llb0_ab(betajk, sjk2, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// llbR0_j
double llbR0_j(NumericVector betajk_j, NumericVector sjk2_j, double lambda, double alpha);
RcppExport SEXP _mamba_llbR0_j(SEXP betajk_jSEXP, SEXP sjk2_jSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j(betajk_jSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j(sjk2_jSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(llbR0_j(betajk_j, sjk2_j, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// deltis
NumericVector deltis(NumericVector betajk_j, NumericVector sjk2_j, double lambda, double alpha);
RcppExport SEXP _mamba_deltis(SEXP betajk_jSEXP, SEXP sjk2_jSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type betajk_j(betajk_jSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sjk2_j(sjk2_jSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(deltis(betajk_j, sjk2_j, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _mamba_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mamba_llbR1_j", (DL_FUNC) &_mamba_llbR1_j, 3},
    {"_mamba_logsumexp", (DL_FUNC) &_mamba_logsumexp, 1},
    {"_mamba_llb0_ab", (DL_FUNC) &_mamba_llb0_ab, 4},
    {"_mamba_llbR0_j", (DL_FUNC) &_mamba_llbR0_j, 4},
    {"_mamba_deltis", (DL_FUNC) &_mamba_deltis, 4},
    {"_mamba_rcpp_hello_world", (DL_FUNC) &_mamba_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_mamba(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
