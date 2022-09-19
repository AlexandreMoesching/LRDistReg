// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calibrate1_cpp
arma::mat calibrate1_cpp(arma::mat& theta, int n, arma::vec& w_jplus);
RcppExport SEXP _LRDistReg_calibrate1_cpp(SEXP thetaSEXP, SEXP nSEXP, SEXP w_jplusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_jplus(w_jplusSEXP);
    rcpp_result_gen = Rcpp::wrap(calibrate1_cpp(theta, n, w_jplus));
    return rcpp_result_gen;
END_RCPP
}
// calibrate2_cpp
arma::mat calibrate2_cpp(arma::mat& theta, int n, arma::vec& w_plusk);
RcppExport SEXP _LRDistReg_calibrate2_cpp(SEXP thetaSEXP, SEXP nSEXP, SEXP w_pluskSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_plusk(w_pluskSEXP);
    rcpp_result_gen = Rcpp::wrap(calibrate2_cpp(theta, n, w_plusk));
    return rcpp_result_gen;
END_RCPP
}
// calibrate_cpp
arma::mat calibrate_cpp(arma::mat& theta, int l, arma::imat& mM, int n, arma::mat& w, arma::vec& w_jplus, arma::vec& w_plusk, double prec);
RcppExport SEXP _LRDistReg_calibrate_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP mMSEXP, SEXP nSEXP, SEXP wSEXP, SEXP w_jplusSEXP, SEXP w_pluskSEXP, SEXP precSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_jplus(w_jplusSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type w_plusk(w_pluskSEXP);
    Rcpp::traits::input_parameter< double >::type prec(precSEXP);
    rcpp_result_gen = Rcpp::wrap(calibrate_cpp(theta, l, mM, n, w, w_jplus, w_plusk, prec));
    return rcpp_result_gen;
END_RCPP
}
// interpolate_cpp
arma::mat interpolate_cpp(arma::vec& x0, arma::vec& x, arma::mat& CDF);
RcppExport SEXP _LRDistReg_interpolate_cpp(SEXP x0SEXP, SEXP xSEXP, SEXP CDFSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type CDF(CDFSEXP);
    rcpp_result_gen = Rcpp::wrap(interpolate_cpp(x0, x, CDF));
    return rcpp_result_gen;
END_RCPP
}
// vgamma_tilde1_cpp
List vgamma_tilde1_cpp(arma::mat& theta, int l, int m, int n, arma::imat& mM, arma::mat& w_ul);
RcppExport SEXP _LRDistReg_vgamma_tilde1_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP nSEXP, SEXP mMSEXP, SEXP w_ulSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w_ul(w_ulSEXP);
    rcpp_result_gen = Rcpp::wrap(vgamma_tilde1_cpp(theta, l, m, n, mM, w_ul));
    return rcpp_result_gen;
END_RCPP
}
// vgamma_tilde2_cpp
List vgamma_tilde2_cpp(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::mat& w_ol);
RcppExport SEXP _LRDistReg_vgamma_tilde2_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP nSEXP, SEXP lLSEXP, SEXP w_olSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type lL(lLSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w_ol(w_olSEXP);
    rcpp_result_gen = Rcpp::wrap(vgamma_tilde2_cpp(theta, l, m, n, lL, w_ol));
    return rcpp_result_gen;
END_RCPP
}
// f_theta_cpp
double f_theta_cpp(arma::mat& theta, int l, int n, arma::imat& mM, arma::mat& w);
RcppExport SEXP _LRDistReg_f_theta_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP nSEXP, SEXP mMSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(f_theta_cpp(theta, l, n, mM, w));
    return rcpp_result_gen;
END_RCPP
}
// simple_step_cpp
arma::mat simple_step_cpp(arma::mat& theta, arma::mat& Psi, double delta, int l, arma::imat& mM, int n, arma::mat& w);
RcppExport SEXP _LRDistReg_simple_step_cpp(SEXP thetaSEXP, SEXP PsiSEXP, SEXP deltaSEXP, SEXP lSEXP, SEXP mMSEXP, SEXP nSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Psi(PsiSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(simple_step_cpp(theta, Psi, delta, l, mM, n, w));
    return rcpp_result_gen;
END_RCPP
}
// local_search1_cpp
List local_search1_cpp(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::imat& mM, arma::mat& w, arma::mat& w_ul);
RcppExport SEXP _LRDistReg_local_search1_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP nSEXP, SEXP lLSEXP, SEXP mMSEXP, SEXP wSEXP, SEXP w_ulSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type lL(lLSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w_ul(w_ulSEXP);
    rcpp_result_gen = Rcpp::wrap(local_search1_cpp(theta, l, m, n, lL, mM, w, w_ul));
    return rcpp_result_gen;
END_RCPP
}
// local_search2_cpp
List local_search2_cpp(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::imat& mM, arma::mat& w, arma::mat& w_ol);
RcppExport SEXP _LRDistReg_local_search2_cpp(SEXP thetaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP nSEXP, SEXP lLSEXP, SEXP mMSEXP, SEXP wSEXP, SEXP w_olSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type lL(lLSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type w_ol(w_olSEXP);
    rcpp_result_gen = Rcpp::wrap(local_search2_cpp(theta, l, m, n, lL, mM, w, w_ol));
    return rcpp_result_gen;
END_RCPP
}
// dist_reg_cpp
List dist_reg_cpp(arma::vec& X, arma::vec& Y, arma::vec& W, double delta0, arma::vec x0, bool ST);
RcppExport SEXP _LRDistReg_dist_reg_cpp(SEXP XSEXP, SEXP YSEXP, SEXP WSEXP, SEXP delta0SEXP, SEXP x0SEXP, SEXP STSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type W(WSEXP);
    Rcpp::traits::input_parameter< double >::type delta0(delta0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< bool >::type ST(STSEXP);
    rcpp_result_gen = Rcpp::wrap(dist_reg_cpp(X, Y, W, delta0, x0, ST));
    return rcpp_result_gen;
END_RCPP
}
// prepare_data_cpp
List prepare_data_cpp(arma::vec& X, arma::vec& Y, arma::vec& W);
RcppExport SEXP _LRDistReg_prepare_data_cpp(SEXP XSEXP, SEXP YSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(prepare_data_cpp(X, Y, W));
    return rcpp_result_gen;
END_RCPP
}
// lambda1_to_theta_cpp
arma::mat lambda1_to_theta_cpp(arma::mat& lambda, int l, int m, arma::imat& mM);
RcppExport SEXP _LRDistReg_lambda1_to_theta_cpp(SEXP lambdaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP mMSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type mM(mMSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda1_to_theta_cpp(lambda, l, m, mM));
    return rcpp_result_gen;
END_RCPP
}
// lambda2_to_theta_cpp
arma::mat lambda2_to_theta_cpp(arma::mat& lambda, int l, int m, arma::imat& lL);
RcppExport SEXP _LRDistReg_lambda2_to_theta_cpp(SEXP lambdaSEXP, SEXP lSEXP, SEXP mSEXP, SEXP lLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type l(lSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type lL(lLSEXP);
    rcpp_result_gen = Rcpp::wrap(lambda2_to_theta_cpp(lambda, l, m, lL));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LRDistReg_calibrate1_cpp", (DL_FUNC) &_LRDistReg_calibrate1_cpp, 3},
    {"_LRDistReg_calibrate2_cpp", (DL_FUNC) &_LRDistReg_calibrate2_cpp, 3},
    {"_LRDistReg_calibrate_cpp", (DL_FUNC) &_LRDistReg_calibrate_cpp, 8},
    {"_LRDistReg_interpolate_cpp", (DL_FUNC) &_LRDistReg_interpolate_cpp, 3},
    {"_LRDistReg_vgamma_tilde1_cpp", (DL_FUNC) &_LRDistReg_vgamma_tilde1_cpp, 6},
    {"_LRDistReg_vgamma_tilde2_cpp", (DL_FUNC) &_LRDistReg_vgamma_tilde2_cpp, 6},
    {"_LRDistReg_f_theta_cpp", (DL_FUNC) &_LRDistReg_f_theta_cpp, 5},
    {"_LRDistReg_simple_step_cpp", (DL_FUNC) &_LRDistReg_simple_step_cpp, 7},
    {"_LRDistReg_local_search1_cpp", (DL_FUNC) &_LRDistReg_local_search1_cpp, 8},
    {"_LRDistReg_local_search2_cpp", (DL_FUNC) &_LRDistReg_local_search2_cpp, 8},
    {"_LRDistReg_dist_reg_cpp", (DL_FUNC) &_LRDistReg_dist_reg_cpp, 6},
    {"_LRDistReg_prepare_data_cpp", (DL_FUNC) &_LRDistReg_prepare_data_cpp, 3},
    {"_LRDistReg_lambda1_to_theta_cpp", (DL_FUNC) &_LRDistReg_lambda1_to_theta_cpp, 4},
    {"_LRDistReg_lambda2_to_theta_cpp", (DL_FUNC) &_LRDistReg_lambda2_to_theta_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_LRDistReg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
