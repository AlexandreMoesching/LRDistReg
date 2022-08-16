// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// prepare_data_cpp
List prepare_data_cpp(arma::vec X, arma::vec Y, arma::vec W);
RcppExport SEXP _LRDistReg_prepare_data_cpp(SEXP XSEXP, SEXP YSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_data_cpp(X, Y, W));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LRDistReg_prepare_data_cpp", (DL_FUNC) &_LRDistReg_prepare_data_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_LRDistReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
