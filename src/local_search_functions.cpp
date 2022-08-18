#include "local_search_functions.h"

void simple_step_ref_cpp(arma::mat& theta, arma::mat& Psi,
                         double& delta, const par& par) {
  // CHECK IF WE NEED THETA_NEW LIKE IN THE R CODE!

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
             t_star  *   Psi.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1));
  }
}

void local_search1_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           arma::mat& v, arma::mat& gamma,
                           arma::mat& lambda_star, pava_par& par1,
                           double& delta, const par& par) {
  // Declare variables
  int lk1, Lk0, d;

  // Compute v-tilde and gamma-tilde
  vgamma_tilde1_ref_cpp(theta, v, gamma, par);

  // Compute lambda_star
  // (i) Baseline
  for (int j = 0; j < par.ell; j++) {
    lambda_star.at(j, par.mM.at(j, 0)) = gamma.at(j, par.mM.at(j, 0));
  }
  // (ii) Isotonic regression
  if (par.m > 1) {
    for (int k = 1; k < par.m; k++) {
      lk1 = par.lL.at(k, 0);
      Lk0 = par.lL.at(k - 1, 1);
      if (lk1 <= Lk0) {
        // Initialize
        d = 1;
        par1.PP[lk1 + d - 1] = lk1 - 1;
        par1.PP[lk1 + d] = lk1;
        par1.WW[lk1 + d] = v.at(lk1, k);
        par1.MM[lk1 + d - 1] = R_NegInf;
        par1.MM[lk1 + d] = gamma.at(lk1, k);

        // PAVA
        for (int j = lk1 + 1; j < Lk0 + 1; j++) {
          // Increase partition
          d++;
          par1.PP[lk1 + d] = j;
          par1.WW[lk1 + d] = v.at(j, k);
          par1.MM[lk1 + d] = gamma.at(j, k);

          // Pool adjacent violators
          while (par1.MM[lk1 + d - 1] >= par1.MM[lk1 + d]) {
            d--;
            par1.MM[lk1 + d] =
              par1.WW[lk1 + d]     * par1.MM[lk1 + d] +
              par1.WW[lk1 + d + 1] * par1.MM[lk1 + d + 1];
            par1.WW[lk1 + d] = par1.WW[lk1 + d] + par1.WW[lk1 + d + 1];
            par1.MM[lk1 + d] = par1.MM[lk1 + d] / par1.WW[lk1 + d];
            par1.PP[lk1 + d] = par1.PP[lk1 + d + 1];
          }
        }

        // Reconstruct solution
        for (int u = 0; u < d; u++) {
          lambda_star.col(
            k
          ).subvec(
              par1.PP[lk1 + u] + 1, par1.PP[lk1 + u + 1]
          ).fill(
              par1.MM[lk1 + u + 1]
          );
        }
      }
    }
  }

  // Transform lambda_star into Psi
  lambda1_to_theta_ref_cpp(lambda_star, Psi, par);

  // Update delta
  delta = 0;
  for (int j = 0; j < par.ell; j++) {
    delta += sum(
      (
                      - par.w.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) +
            par.n * exp(theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)))
      ) % (
          theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) -
            Psi.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
      )
    );
  }
}

void local_search2_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           arma::mat& v, arma::mat& gamma,
                           arma::mat& lambda_star, pava_par& par2,
                           double& delta, const par& par) {
  // Declare variables
  int mj1, Mj0, d;

  // Compute v-tilde and gamma-tilde
  vgamma_tilde2_ref_cpp(theta, v, gamma, par);

  // Compute lambda_star
  // (i) Baseline
  for (int k = 0; k < par.m; k++) {
    lambda_star.at(par.lL.at(k, 0), k) = gamma.at(par.lL.at(k, 0), k);
  }
  // (ii) Isotonic regression
  if (par.ell > 1) {
    for (int j = 1; j < par.ell; j++) {
      mj1 = par.mM.at(j, 0);
      Mj0 = par.mM.at(j - 1, 1);
      if (mj1 <= Mj0) {
        // Initialize
        d = 1;
        par2.PP[mj1 + d - 1] = mj1 - 1;
        par2.PP[mj1 + d] = mj1;
        par2.WW[mj1 + d] = v.at(j, mj1);
        par2.MM[mj1 + d - 1] = R_NegInf;
        par2.MM[mj1 + d] = gamma.at(j, mj1);

        // PAVA
        for (int k = mj1 + 1; k < Mj0 + 1; k++) {
          // Increase partition
          d++;
          par2.PP[mj1 + d] = k;
          par2.WW[mj1 + d] = v.at(j, k);
          par2.MM[mj1 + d] = gamma.at(j, k);

          // Pool adjacent violators
          while (par2.MM[mj1 + d - 1] >= par2.MM[mj1 + d]) {
            d--;
            par2.MM[mj1 + d] =
              par2.WW[mj1 + d]     * par2.MM[mj1 + d] +
              par2.WW[mj1 + d + 1] * par2.MM[mj1 + d + 1];
            par2.WW[mj1 + d] = par2.WW[mj1 + d] + par2.WW[mj1 + d + 1];
            par2.MM[mj1 + d] = par2.MM[mj1 + d] / par2.WW[mj1 + d];
            par2.PP[mj1 + d] = par2.PP[mj1 + d + 1];
          }
        }

        // Reconstruct solution
        for (int u = 0; u < d; u++) {
          lambda_star.row(
            j
          ).subvec(
              par2.PP[mj1 + u] + 1, par2.PP[mj1 + u + 1]
          ).fill(
              par2.MM[mj1 + u + 1]
          );
        }
      }
    }
  }

  // Transform lambda_star into Psi
  lambda2_to_theta_ref_cpp(lambda_star, Psi, par);

  // Update delta
  delta = 0;
  for (int j = 0; j < par.ell; j++) {
    delta += sum(
      (
          - par.w.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) +
            par.n * exp(theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)))
      ) % (
          theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) -
            Psi.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
      )
    );
  }
}

//' Function to take a simple step
//'
//' @param theta Log-parameter
//' @param Psi Proposal
//' @param delta Delta
//' @param ell Number of unique covariates
//' @param mM (m_j,M_j) index pairs
//' @param n Sample size
//' @param w Sample weights
//'
//' @return Updated theta parameter
//'
//' @export
//[[Rcpp::export]]
arma::mat simple_step_cpp(arma::mat theta, arma::mat Psi, double delta,
                          int ell, arma::imat mM, int n, arma::mat w) {
  // Declare variables
  par par;
  par.ell = ell;
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
//' @param w Sample weights
//' @param w_ul Cumulative row-sums of sample weights
//'
//' @return New proposal Psi and step-size delta
//'
//' @export
//[[Rcpp::export]]
List local_search1_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM,
                       arma::mat w, arma::mat w_ul) {
  // Declare variables
  arma::mat Psi(ell, m, arma::fill::value(R_NegInf));
  arma::mat v(ell, m), gamma(ell, m), lambda_star(ell, m);
  arma::ivec PP1(ell + 1);
  arma::vec MM1(ell + 1);
  arma::vec WW1(ell + 1);
  pava_par par1 = {PP1, MM1, WW1};
  par par;
  par.ell = ell;
  par.lL = lL;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.w = w;
  par.w_ul = w_ul;
  double delta;

  // Compute Psi and delta
  local_search1_ref_cpp(theta, Psi, v, gamma, lambda_star, par1, delta, par);

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
//' @param w Sample weights
//' @param w_ol Cumulative column-sums of sample weights
//'
//' @return New proposal Psi and step-size delta
//'
//' @export
//[[Rcpp::export]]
List local_search2_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM,
                       arma::mat w, arma::mat w_ol) {
  // Declare variables
  arma::mat Psi(ell, m, arma::fill::value(R_NegInf));
  arma::mat v(ell, m), gamma(ell, m), lambda_star(ell, m);
  arma::ivec PP2(m + 1);
  arma::vec MM2(m + 1);
  arma::vec WW2(m + 1);
  pava_par par2 = {PP2, MM2, WW2};
  par par;
  par.ell = ell;
  par.lL = lL;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.w = w;
  par.w_ol = w_ol;
  double delta;

  // Compute Psi and delta
  local_search2_ref_cpp(theta, Psi, v, gamma, lambda_star, par2, delta, par);

  // Return
  return List::create(Named("Psi") = Psi,
                      Named("delta") = delta);
}
