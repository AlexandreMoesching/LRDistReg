#include "reparametrize.h"

void lambda1_to_theta_ref_cpp(arma::mat& lambda, arma::mat& theta, const par& par) {
  // Reparametrize
  for (int j = 0; j < par.ell; j++) {
              theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) =
      cumsum(lambda.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)));
  }
}

void lambda2_to_theta_ref_cpp(arma::mat& lambda, arma::mat& theta, const par& par) {
  // Reparametrize
  for (int k = 0; k < par.m; k++) {
              theta.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)) =
      cumsum(lambda.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)));
  }
}

//' Transforms lambda (row) into theta, C++ version
//'
//' @param lambda Row-wise differences
//' @param ell Number of unique covariates
//' @param m Number of unique responses
//' @param mM (m_j,M_j) index pairs
//'
//' @return Transformed parameter
//'
//' @export
// [[Rcpp::export]]
arma::mat lambda1_to_theta_cpp(arma::mat lambda, int ell, int m, arma::imat mM) {
  // Declare variables
  par par;
  par.ell = ell;
  par.m = m;
  par.mM = mM;
  arma::mat theta(ell, m, arma::fill::value(R_NegInf));

  // Transform parameter
  lambda1_to_theta_ref_cpp(lambda, theta, par);

  // Return
  return theta;
}

//' Transforms lambda (column) into theta, C++ version
//'
//' @param lambda Column-wise differences
//' @param ell Number of unique covariates
//' @param m Number of unique responses
//' @param lL (l_k,L_k) index pairs
//'
//' @return Transformed parameter
//'
//' @export
// [[Rcpp::export]]
arma::mat lambda2_to_theta_cpp(arma::mat lambda, int ell, int m, arma::imat lL) {
  // Declare variables
  par par;
  par.ell = ell;
  par.m = m;
  par.lL = lL;
  arma::mat theta(ell, m, arma::fill::value(R_NegInf));

  // Transform parameter
  lambda2_to_theta_ref_cpp(lambda, theta, par);

  // Return
  return theta;
}
