#include "calibrate.h"

void calibrate1_ref_cpp(arma::mat& theta, par& par) {
  theta.each_col() += - log(sum(exp(theta), 1)) + log(par.w_jplus / par.n);
}

void calibrate2_ref_cpp(arma::mat& theta, par& par) {
  theta.each_row() += - log(sum(exp(theta), 0)) + log(par.w_plusk.t() / par.n);
}

void calibrate_ref_cpp(arma::mat& theta, par& par, double& prec) {
  // Declare variable
  double f_theta, f_theta_old = f_theta_ref_cpp(theta, par);
  double delta = R_PosInf;

  // While loop
  while (delta > prec) {
    // Row-wise calibration
    calibrate1_ref_cpp(theta, par);

    // Column-wise calibration
    calibrate2_ref_cpp(theta, par);

    // Update f_theta and delta
    f_theta = f_theta_ref_cpp(theta, par);
    delta = f_theta_old - f_theta;
    f_theta_old = f_theta;
  }
}

//' Row-calibration of log-parameter
//'
//' @param theta Log-parameter
//' @param n Sample size
//' @param w_jplus Row-sums of sample weights
//'
//' @return Row-calibrated log-parameter
//'
//' @export
// [[Rcpp::export]]
arma::mat calibrate1_cpp(arma::mat theta, int n, arma::vec w_jplus) {
  // Declare variables
  par par;
  par.n = n;
  par.w_jplus = w_jplus;

  // Calibrate theta
  calibrate1_ref_cpp(theta, par);

  // Return
  return theta;
}

//' Column-calibration of log-parameter
//'
//' @param theta Log-parameter
//' @param n Sample size
//' @param w_plusk Column-sums of sample weights
//'
//' @return Row-calibrated log-parameter
//'
//' @export
// [[Rcpp::export]]
arma::mat calibrate2_cpp(arma::mat theta, int n, arma::vec w_plusk) {
  // Declare variables
  par par;
  par.n = n;
  par.w_plusk = w_plusk;

  // Calibrate theta
  calibrate2_ref_cpp(theta, par);

  // Return
  return theta;
}

//' Approximate calibration of log-parameter
//'
//' @description Computes an optimal additive adjustment of theta, i.e. a matrix
//' theta.new with values
//'     theta.new_jk = theta_jk + x_j + y_k
//' such that f(theta.new) is minimal.
//'
//' @param theta Log-parameter
//' @param ell Number of unique covariates
//' @param mM (m_j,M_j) index pairs
//' @param n Sample size
//' @param w Sample weights
//' @param w_jplus Row-sums of sample weights
//' @param w_plusk Column sums of sample weights
//' @param prec Precision for the approximate calibration
//'
//' @return Approximately calibrated log-parameter
//'
//' @export
// [[Rcpp::export]]
arma::mat calibrate_cpp(arma::mat theta, int ell, arma::imat mM, int n,
                        arma::mat w, arma::vec w_jplus, arma::vec w_plusk,
                        double prec) {
  // Declare variables
  par par;
  par.ell = ell;
  par.mM = mM;
  par.n = n;
  par.w = w;
  par.w_jplus = w_jplus;
  par.w_plusk = w_plusk;

  // Calibrate theta
  calibrate_ref_cpp(theta, par, prec);

  // Return
  return theta;
}
