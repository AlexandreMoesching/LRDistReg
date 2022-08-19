#include "main.h"

void TP2_fit_ref_cpp(arma::mat& h_TP2, arma::mat& theta, arma::mat& Psi,
                     arma::mat& v, arma::mat& gamma, arma::mat& lambda_star,
                     pava_par& par1, pava_par& par2,
                     double& delta, double delta0, par& par) {
  // Declare variables
  double tmp_dbl = -log(accu(par.PP)); // SHOULD WE ALSO REPLACE SUM BY ACCU ELSEWHERE ??
  double prec = 1e-5;
  int s = 0;

  // Initialize delta and theta
  delta = R_PosInf;
  for (int j = 0; j < par.ell; j++) {
    theta.row(j).subvec(par.mM.at(j, 0), par.mM.at(j, 1)).fill(tmp_dbl);
  }

  // Main while-loop
  while (delta > delta0) {
    // Calibrate theta
    calibrate_ref_cpp(theta, par, prec);
    std::cout << "Iteration = " << s << "\n";

    // Find a new proposal
    if (s % 2 == 0) {
      local_search1_ref_cpp(theta, Psi, v, gamma, lambda_star, par1, delta, par);
    } else {
      local_search2_ref_cpp(theta, Psi, v, gamma, lambda_star, par2, delta, par);
    }

    // Perform real improvement
    simple_step_ref_cpp(theta, Psi, delta, par);
    std::cout << "delta = " << delta << "\n\n";

    // Change parity
    s++;
  }

  std::cout << "delta is smaller than delta0 = " << delta0 << "\n";

  // Return probability weights
  h_TP2 = exp(theta);
}


//' TP2 fit function
//'
//' @param X Covariates
//' @param Y Responses
//' @param W User-specified weights
//' @param delta0 Threshhold
//'
//' @return h matrix, delta and parameters (estimation time)
//'
//' @export
// [[Rcpp::export]]
List TP2_fit_cpp(arma::vec X, arma::vec Y, arma::vec W, double delta0) {
  // Compute parameters
  par par = prepare_data_par_cpp(X, Y, W);
  List par_list = List::create(Named("ell") = par.ell,
                               Named("lL") = par.lL,
                               Named("m") = par.m,
                               Named("mM") = par.mM,
                               Named("n") = par.n,
                               Named("PP") = par.PP,
                               Named("w") = par.w,
                               Named("w_jplus") = par.w_jplus,
                               Named("w_plusk") = par.w_plusk,
                               Named("w_ul") = par.w_ul,
                               Named("w_ol") = par.w_ol,
                               Named("W") = par.W,
                               Named("x") = par.x,
                               Named("X") = par.X,
                               Named("y") = par.y,
                               Named("Y") = par.Y);

  // Declare variables
  arma::mat h_TP2(par.ell, par.m);
  arma::mat theta(par.ell, par.m, arma::fill::value(R_NegInf));
  arma::mat Psi(par.ell, par.m);
  arma::mat v(par.ell, par.m);
  arma::mat gamma(par.ell, par.m);
  arma::mat lambda_star(par.ell, par.m);

  arma::ivec PP1(par.ell + 1);
  arma::vec MM1(par.ell + 1);
  arma::vec WW1(par.ell + 1);
  pava_par par1 = {PP1, MM1, WW1};

  arma::ivec PP2(par.m + 1);
  arma::vec MM2(par.m + 1);
  arma::vec WW2(par.m + 1);
  pava_par par2 = {PP2, MM2, WW2};

  double delta;

  // Estimate
  TP2_fit_ref_cpp(h_TP2, theta, Psi, v, gamma, lambda_star,
                  par1, par2, delta, delta0, par);

  // Return
  return List::create(Named("h_TP2") = h_TP2,
                      Named("delta") = delta,
                      Named("par") = par_list);
}
