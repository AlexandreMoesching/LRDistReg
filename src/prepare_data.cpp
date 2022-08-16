#include "prepare_data.h"

void prepare_data_void_cpp(par& par) {
  // Number of observations
  par.n = par.X.n_elem;

  // Sorted values of X and Y
  arma::vec sorted_X = sort(par.X);
  arma::vec sorted_Y = sort(par.Y);

  // New values of from the sorted X and Y
  arma::uvec new_X(par.n, arma::fill::zeros);
  arma::uvec new_Y(par.n, arma::fill::zeros);
  new_X[0] = 1;
  new_Y[0] = 1;
  for (int i = 1; i < par.n; i++) {
    if (sorted_X[i-1] < sorted_X[i]) {
      new_X[i] = 1;
    }
    if (sorted_Y[i-1] < sorted_Y[i]) {
      new_Y[i] = 1;
    }
  }

  // Number of unique X and Y
  par.ell = sum(new_X);
  par.m = sum(new_Y);

  // Unique values of X and Y
  arma::vec x(par.ell);
  arma::vec y(par.m);
  int j = 0, k = 0;
  for (int i = 0; i < par.n; i++) {
    if (new_X[i] == 1) {
      x[j] = sorted_X[i];
      j++;
    }
    if (new_Y[i] == 1) {
      y[k] = sorted_Y[i];
      k++;
    }
  }
  par.x = x;
  par.y = y;

  // Index of x/y corresponding to a given X_i/Y_i
  // Function 'find' requires uvec type
  arma::uvec xk_Xi;
  arma::uvec yj_Yi;

  // Weight matrix
  arma::mat w(par.ell, par.m, arma::fill::zeros);
  for (int i = 0; i < par.n; i++) {
    xk_Xi = find(par.X[i] == x);
    yj_Yi = find(par.Y[i] == y);
    w.at(xk_Xi[0], yj_Yi[0]) += par.W[i];
  }
  par.w = w;

  // Weights w_j+ = sum_k wjk and w_+k = sum_j wjk
  par.w_jplus = sum(w, 1);
  par.w_plusk = sum(w, 0).t();

  // Reverted cumulative row and column sums of weight matrix
  par.w_ul = fliplr(cumsum(fliplr(w), 1));
  par.w_ol = flipud(cumsum(flipud(w), 0));

  // lL
  arma::imat lL(par.m, 2, arma::fill::zeros);
  int tmp, lk = par.ell - 1, Lk = 0;
  for (int k = 0; k < par.m; k++) {
    tmp = min(find(w.col(par.m - 1 - k) > 0));
    if (tmp < lk) {
      lk = tmp;
    }
    lL.at(par.m - 1 - k, 0) = lk;

    tmp = max(find(w.col(k) > 0));
    if (tmp > Lk) {
      Lk = tmp;
    }
    lL.at(k, 1) = Lk;
  }
  par.lL = lL;

  // mM
  arma::imat mM(par.ell, 2, arma::fill::zeros);
  int mj = par.m - 1, Mj = 0;
  for (int j = 0; j < par.ell; j++) {
    tmp = min(find(w.row(par.ell - 1 - j) > 0));
    if (tmp < mj) {
      mj = tmp;
    }
    mM.at(par.ell - 1 - j, 0) = mj;

    tmp = max(find(w.row(j) > 0));
    if (tmp > Mj) {
      Mj = tmp;
    }
    mM.at(j, 1) = Mj;
  }
  par.mM = mM;

  // PP
  arma::imat PP(par.ell, par.m, arma::fill::zeros);
  for (int j = 0; j < par.ell; j++) {
    PP.row(j).subvec(mM.at(j, 0), mM.at(j, 1)).fill(1);
  }
  par.PP = PP;
}

//' Prepare the data, C++ version
//'
//' @param X Covariates
//' @param Y Responses
//' @param W User-specified weights
//'
//' @return A list of pre-computed parameters necessary for the estimation
//'
//' @export
// [[Rcpp::export]]
List prepare_data_cpp(arma::vec X, arma::vec Y, arma::vec W) {
  // Declare variables
  par par;
  par.X = X;
  par.Y = Y;
  par.W = W;

  // Compute parameters
  prepare_data_void_cpp(par);

  // Return
  return List::create(Named("ell") = par.ell,
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
}
