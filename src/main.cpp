#include "main.h"

void TP2_fit_ref_cpp(arma::mat& h_TP2, arma::mat& q_LR, arma::mat& CDF_LR,
                     arma::mat& theta, arma::mat& Psi,
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
    // std::cout << "Iteration = " << s << "\n";

    // Find a new proposal
    if (s % 2 == 0) {
      local_search1_ref_cpp(theta, Psi, v, gamma, lambda_star, par1, delta, par);
    } else {
      local_search2_ref_cpp(theta, Psi, v, gamma, lambda_star, par2, delta, par);
    }

    // Perform real improvement
    simple_step_ref_cpp(theta, Psi, delta, par);
    // std::cout << "delta = " << delta << "\n\n";

    // Change parity
    s++;
  }

  // std::cout << "delta is smaller than delta0 = " << delta0 << "\n";

  // Compute probability weights
  h_TP2 = exp(theta);

  // Compute q_LR
  for (int j = 0; j < par.ell; j++) {
    q_LR.row(j) = h_TP2.row(j) / sum(h_TP2.row(j));
  }

  // Compute CDF_LR
  CDF_LR = cumsum(q_LR, 1);
}

void ST_fit_ref_cpp(arma::mat& CDF_EMP, arma::mat& CDF_ST,
                    pava_par& par3, par& par) {
  // Compute empirical distribution function
  CDF_EMP = cumsum(par.w, 1);
  CDF_EMP.each_col() %= 1.0 / par.w_jplus;

  // Compute isotonic distributional regression (under ST-order)
  // That is, for each k = 0, ..., m-1, determine the antitonic vector f
  // minimizing:
  //      sum_{j = 0}^{m - 1} w_jplus[j] (f[j] - CDP_EMP.at(j, k))^2.
  // In other words, f is the antitonic regression of CDP_EMP.col(k) with
  // weights w_jplus.
  //
  // This could be done more efficiently, see
  //    Henzi, A., M\"osching, A. & D\"umbgen, L.
  //    Accelerating the Pool-Adjacent-Violators Algorithm for Isotonic
  //    Distributional Regression.
  //    Methodol Comput Appl Probab (2022).
  //    https://doi.org/10.1007/s11009-022-09937-2
  int d;
  for (int k = 0; k < par.m; k++) {
    // Initialize
    d = 1;
    par3.PP[0] = - 1;
    par3.PP[d] =   0;
    par3.WW[d] = par.w_jplus[0];
    par3.MM[0] = R_PosInf;
    par3.MM[d] = CDF_EMP.at(0, k);

    // PAVA (antitonic regression)
    for (int j = 1; j < par.ell; j++) {
      // Increase partition
      d++;
      par3.PP[d] = j;
      par3.WW[d] = par.w_jplus[j];
      par3.MM[d] = CDF_EMP.at(j, k);

      // Pool adjacent violators
      while (par3.MM[d - 1] <= par3.MM[d]) {
        d--;
        par3.MM[d] = par3.WW[d] * par3.MM[d] + par3.WW[d + 1] * par3.MM[d + 1];
        par3.WW[d] = par3.WW[d] + par3.WW[d + 1];
        par3.MM[d] = par3.MM[d] / par3.WW[d];
        par3.PP[d] = par3.PP[d + 1];
      }
    }

    // Reconstruct solution
    for (int u = 0; u < d; u++) {
      CDF_ST.col(
        k
      ).subvec(
          par3.PP[u] + 1, par3.PP[u + 1]
      ).fill(
          par3.MM[u + 1]
      );
    }
  }
}

//' Isotonic distributional regression (LR, ST, EMP)
//'
//' @param X Covariates
//' @param Y Responses
//' @param W User-specified sample weights
//' @param delta0 Threshhold
//' @param ST Boolean indicating whether or not the classical isotonic
//' distributional regression will also be computed
//'
//' @return Isotonic distributional regression(s) and estimation parameters
//'
//' @export
// [[Rcpp::export]]
List dist_reg_cpp(arma::vec X, arma::vec Y, arma::vec W, double delta0,
                  bool ST = false) {
  // Compute parameters
  par par = prepare_data_par_cpp(X, Y, W);

  // Prepare parameter list to return
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

  // Declare variables for TP2 estimation
  arma::mat  h_TP2(par.ell, par.m);
  arma::mat   q_LR(par.ell, par.m);
  arma::mat CDF_LR(par.ell, par.m);

  arma::mat       theta(par.ell, par.m, arma::fill::value(R_NegInf));
  arma::mat         Psi(par.ell, par.m);
  arma::mat           v(par.ell, par.m);
  arma::mat       gamma(par.ell, par.m);
  arma::mat lambda_star(par.ell, par.m);

  arma::ivec PP1(par.ell + 1);
  arma::vec  MM1(par.ell + 1);
  arma::vec  WW1(par.ell + 1);
  pava_par par1 = {PP1, MM1, WW1};

  arma::ivec PP2(par.m + 1);
  arma::vec  MM2(par.m + 1);
  arma::vec  WW2(par.m + 1);
  pava_par par2 = {PP2, MM2, WW2};

  double delta;

  // Estimate TP2 distribution & LR-ordered family of distributions
  TP2_fit_ref_cpp(h_TP2, q_LR, CDF_LR,
                  theta, Psi, v, gamma, lambda_star,
                  par1, par2, delta, delta0, par);


  if (ST) {
    // Declare variables for isotonic distributional regression
    arma::mat CDF_EMP(par.ell, par.m);
    arma::mat  CDF_ST(par.ell, par.m);

    arma::ivec PP3(par.ell + 1);
    arma::vec  MM3(par.ell + 1);
    arma::vec  WW3(par.ell + 1);
    pava_par par3 = {PP3, MM3, WW3};

    // Estimate ST-ordered family of distributions
    ST_fit_ref_cpp(CDF_EMP, CDF_ST, par3, par);

    // Return
    return List::create(Named("h_TP2") = h_TP2,
                        Named("q_LR") = q_LR,
                        Named("CDF_LR") = CDF_LR,
                        Named("CDF_ST") = CDF_ST,
                        Named("CDF_EMP") = CDF_EMP,
                        Named("delta") = delta,
                        Named("par") = par_list);
  } else {
    // Return
    return List::create(Named("h_TP2") = h_TP2,
                        Named("q_LR") = q_LR,
                        Named("CDF_LR") = CDF_LR,
                        Named("delta") = delta,
                        Named("par") = par_list);
  }
}
