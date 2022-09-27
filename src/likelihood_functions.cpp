#include "likelihood_functions.h"

void vg1_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par) {
  for (int j = 0; j < par.l; j++) {
    // Compute v
    v.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) =
      par.n *
      reverse(
        cumsum(
          reverse(
            exp(
              theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
            )
          )
        )
      );

    // Compute gamma using v
    g.at(j, par.mM.at(j, 0)) = theta.at(j, par.mM.at(j, 0));
    if (par.mM.at(j, 0) < par.mM.at(j, 1)) {
      g.row(j).subvec(par.mM.at(j, 0) + 1, par.mM.at(j, 1)) =
        diff(theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)));
    }
    g.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) +=
      - 1.0 + (
          par.w_ul.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) /
            v.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
      );
  }
}

void vg2_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par) {
  // Compute v and gamma
  for (int k = 0; k < par.m; k++) {
    // Compute v
    v.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)) =
      par.n *
      reverse(
        cumsum(
          reverse(
            exp(
              theta.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1))
            )
          )
        )
      );

    // Compute gamma using v
    g.at(par.lL.at(k, 0), k) = theta.at(par.lL.at(k, 0), k);
    if (par.lL.at(k, 0) < par.lL.at(k, 1)) {
      g.col(k).subvec(par.lL.at(k, 0) + 1, par.lL.at(k, 1)) =
        diff(theta.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)));
    }
    g.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)) +=
      - 1.0 + (
          par.w_ol.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1)) /
            v.col(k).subvec(par.lL.at(k, 0), par.lL.at(k, 1))
      );
  }
}

double ftheta_ref(arma::mat& theta, const par& par) {
  // Declare variable
  double ftheta = 0.0;

  // Update ftheta
  for (int j = 0; j < par.l; j++) {
    ftheta += sum(
      - (
          par.w.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)) %
          theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
      ) + (
          par.n * exp(
              theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1))
          )
      )
    );
  }

  // Return
  return ftheta;
}

//' v-tilde and gamma-tilde functions (row), C++ version
//'
//' @param theta Log-parameter
//' @param l Number of unique covariates
//' @param m Number of unique responses
//' @param n Sample size
//' @param mM (m_j,M_j) index pairs
//' @param w_ul Cumulative row-sums of sample weights
//'
//' @return v-tilde and gamma-tilde functions
//'
//' @export
//[[Rcpp::export]]
List vg1_C(arma::mat& theta, int l, int m, int n, arma::imat& mM, arma::mat& w_ul) {
  // Declare variables
  par par;
  par.l = l;
  par.m = m;
  par.mM = mM;
  par.n = n;
  par.w_ul = w_ul;
  arma::mat v(l, m, arma::fill::zeros), g(l, m, arma::fill::zeros);

  // Compute v and gamma
  vg1_ref(theta, v, g, par);

  // Return
  return List::create(Named("v") = v, Named("g") = g);
}

//' v-tilde and gamma-tilde functions (column), C++ version
//'
//' @param theta Log-parameter
//' @param l Number of unique covariates
//' @param m Number of unique responses
//' @param n Sample size
//' @param lL (l_k,L_k) index pairs
//' @param w_ol Cumulative column-sums of sample weights
//'
//' @return v-tilde and gamma-tilde functions
//'
//' @export
//[[Rcpp::export]]
List vg2_C(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::mat& w_ol) {
  // Declare variables
  par par;
  par.l = l;
  par.lL = lL;
  par.m = m;
  par.n = n;
  par.w_ol = w_ol;
  arma::mat v(l, m, arma::fill::zeros);
  arma::mat g(l, m, arma::fill::zeros);

  // Compute v and gamma
  vg2_ref(theta, v, g, par);

  // Return
  return List::create(Named("v") = v, Named("g") = g);
}

//' Negative log-likelihood in terms of log-parameter, C++ version
//'
//' @param theta Log-parameter
//' @param l Number of unique covariates
//' @param n Sample size
//' @param mM (m_j,M_j) index pairs
//' @param w Sample weights
//'
//' @return Negative log-likelihood in terms of log-parameter
//'
//' @export
//[[Rcpp::export]]
double ftheta_C(arma::mat& theta, int l, int n, arma::imat& mM, arma::mat& w) {
  // Declare variables
  par par;
  par.l = l;
  par.mM = mM;
  par.n = n;
  par.w = w;

  // Return ftheta
  return ftheta_ref(theta, par);
}
