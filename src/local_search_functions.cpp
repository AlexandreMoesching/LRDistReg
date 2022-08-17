#include "local_search_functions.h"

void antitonic_reg_ref_cpp(arma::vec& x, arma::vec& gamma, arma::vec& v) {

}

void simple_step_ref_cpp(arma::mat& theta, arma::mat& Psi,
                         double& delta, const par& par) {
  // Declare variables
  double rho = delta, t_star;
  double f_old = f_theta_ref_cpp(theta, par);
  double f_new = f_theta_ref_cpp(Psi, par);

  // While-loop
  while (f_new > f_old) {
    Psi = (theta + Psi) / 2.0;
    rho = rho / 2.0;
    f_new = f_theta_ref_cpp(Psi, par);
  }

  // Compute t_star
  t_star = rho / (2 * (rho - f_old + f_new));
  if (t_star > 1.0) {
    t_star = 1.0;
  }

  // Move theta towards Psi by an amount of t_star
  for (int j = 0; j < par.ell; j++) {
    theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) =
      (1.0 - t_star) * theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) +
      t_star * Psi.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1));
  }
}

void local_search1_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           double& delta, const par& par) {

}

void local_search2_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           double& delta, const par& par) {

}

//' Computes the antitonic regression of a vector gamma with weights v
//'
//' @param gamma Vector to be regressed
//' @param v Weights
//' @param size Length of each vectors
//'
//' @return Antitonic regression of a vector gamma with weights v
//'
//' @export
//[[Rcpp::export]]
arma::vec antitonic_reg_cpp(arma::vec gamma, arma::vec v, int size) {
  // Declare variables
  arma::vec x(size);

  // Compute antitonic regression of gamma with weights v
  antitonic_reg_ref_cpp(x, gamma, v);

  // Return
  return x;
}

//' Function to take a simple step
//'
//' @param theta Log-parameter
//' @param Psi Proposal
//' @param delta Delta
//' @param ell Number of unique covariates
//' @param m Number of unique responses
//' @param n Sample size
//' @param PP Reduced index space
//' @param w Sample weights
//'
//' @return Updated theta parameter
//'
//' @export
//[[Rcpp::export]]
arma::mat simple_step_cpp(arma::mat theta, arma::mat Psi, double delta,
                          int ell, int m, arma::imat mM, int n, int w) {
  // Declare variables
  par par;
  par.ell = ell;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.w = w;

  // Update theta
  simple_step_ref_cpp(theta, Psi, delta, par);

  // Return
  return theta;
}

//' Local search (row)
//'
//' @param theta Log-parameter
//' @param ell Number of unique covariates
//' @param m Number of unique responses
//' @param n Sample size
//' @param lL (l_k,L_k) index pairs
//' @param mM (m_j,M_j) index pairs
//' @param PP Reduced index space
//' @param w Sample weights
//' @param w_ul Cumulative row-sums of sample weights
//'
//' @return New proposal Psi and step-size delta
//'
//' @export
//[[Rcpp::export]]
List local_search1_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM, arma::imat PP,
                       arma::mat w, arma::mat w_ul) {
  // Declare variables
  par par;
  par.ell = ell;
  par.lL = lL;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.PP = PP;
  par.w = w;
  par.w_ul = w_ul;
  arma::mat Psi(ell, m, arma::fill::zeros);
  double delta;

  // Compute Psi and delta
  local_search1_ref_cpp(theta, Psi, delta, par);

  // Return
  return List::create(Named("Psi") = Psi,
                      Named("delta") = delta);
}

//' Local search (column)
//'
//' @param theta Log-parameter
//' @param ell Number of unique covariates
//' @param m Number of unique responses
//' @param n Sample size
//' @param lL (l_k,L_k) index pairs
//' @param mM (m_j,M_j) index pairs
//' @param PP Reduced index space
//' @param w Sample weights
//' @param w_ol Cumulative column-sums of sample weights
//'
//' @return New proposal Psi and step-size delta
//'
//' @export
//[[Rcpp::export]]
List local_search2_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM, arma::imat PP,
                       arma::mat w, arma::mat w_ol) {
  // Declare variables
  par par;
  par.ell = ell;
  par.lL = lL;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.PP = PP;
  par.w = w;
  par.w_ol = w_ol;
  arma::mat Psi(ell, m, arma::fill::zeros);
  double delta;

  // Compute Psi and delta
  local_search2_ref_cpp(theta, Psi, delta, par);

  // Return
  return List::create(Named("Psi") = Psi,
                      Named("delta") = delta);
}
