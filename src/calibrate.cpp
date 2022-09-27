#include "calibrate.h"

void calibrate_ref(arma::mat& theta, const par& par, double& prec) {
  // Declare variable
  double f_new, f_old = ftheta_ref(theta, par);
  double delta = R_PosInf;

  // While loop
  while (delta > prec) {
    // Row-wise calibration
    theta.each_col() += - log(sum(exp(theta), 1)) + log(par.w_jplus / par.n);

    // Column-wise calibration
    theta.each_row() += - log(sum(exp(theta), 0)) + log(par.w_plusk.t() / par.n);

    // Update ftheta and delta
    f_new = ftheta_ref(theta, par);
    delta = f_old - f_new;
    f_old = f_new;
  }
}

//' Approximate calibration of log-parameter, C++ version
//'
//' @description Computes an optimal additive adjustment of theta, i.e. a matrix
//' theta.new with values
//'     theta.new_jk = theta_jk + x_j + y_k
//' such that f(theta.new) is minimal.
//'
//' @param theta Log-parameter
//' @param l Number of unique covariates
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
arma::mat calibrate_C(arma::mat& theta,
                      int l, arma::imat& mM, int n,
                      arma::mat& w, arma::vec& w_jplus, arma::vec& w_plusk,
                      double prec) {
  // Declare variables
  par par;
  par.l = l;
  par.mM = mM;
  par.n = n;
  par.w = w;
  par.w_jplus = w_jplus;
  par.w_plusk = w_plusk;

  // Calibrate theta
  calibrate_ref(theta, par, prec);

  // Return
  return theta;
}
