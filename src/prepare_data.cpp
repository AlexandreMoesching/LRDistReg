#include "prepare_data.h"

par prepare_data_par(arma::vec& X, arma::vec& Y, arma::vec& W) {
  // Declare variables
  par par;
  par.X = X;
  par.Y = Y;
  par.W = W;
  par.n = X.n_elem;

  // Preprocess X and Y:
  arma::uvec ordX = stable_sort_index(X); // Sorted indices of elements in X
  arma::uvec ordY = stable_sort_index(Y); // Sorted indices of elements in Y

  arma::vec Xs(par.n); // Sorted X's with repetition
  arma::vec Ys(par.n); // Sorted Y's with repetition
  for (int i = 0; i < par.n; i++) {
    Xs[i] = X[ordX[i]];
    Ys[i] = Y[ordY[i]];
  }

  int l = 0;
  int m   = 0;

  // Xj_tmp[i] gives the index j of x such that x[j] = Xs[i]
  // Yk_tmp[i] gives the index k of y such that y[k] = Ys[i]
  arma::vec Xj_tmp(par.n, arma::fill::zeros);
  arma::vec Yk_tmp(par.n, arma::fill::zeros);

  arma::vec x_tmp(par.n, arma::fill::value(Xs[0])); // Unique sorted X's
  arma::vec y_tmp(par.n, arma::fill::value(Ys[0])); // Unique sorted Y's

  for (int i = 1; i < par.n; i++) {
    if (Xs[i] > x_tmp[l]) {
      l++;
      x_tmp[l] = Xs[i];
    }
    Xj_tmp[i] = l;
    if (Ys[i] > y_tmp[m]) {
      m++;
      y_tmp[m] = Ys[i];
    }
    Yk_tmp[i] = m;
  }

  par.x = x_tmp.subvec(0, l); // Sorted unique covariates
  par.y = y_tmp.subvec(0, m); // Sorted unique responses

  par.l = l + 1; // Number of unique covariates
  par.m   = m   + 1; // Number of unique responses

  // Xj[i] gives the index j of x such that x[j] = X[i]
  // Yk[i] gives the index k of y such that y[k] = Y[i]
  arma::uvec Xj(par.n);
  arma::uvec Yk(par.n);
  for (int i = 0; i < par.n; i++) {
    Xj[ordX[i]] = Xj_tmp[i];
    Yk[ordY[i]] = Yk_tmp[i];
  }

  // Weight matrix
  arma::mat w(par.l, par.m, arma::fill::zeros);
  for (int i = 0; i < par.n; i++) {
    w.at(Xj[i], Yk[i]) += W[i];
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
  int tmp, lk = par.l - 1, Lk = 0;
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
  arma::imat mM(par.l, 2, arma::fill::zeros);
  int mj = par.m - 1, Mj = 0;
  for (int j = 0; j < par.l; j++) {
    tmp = min(find(w.row(par.l - 1 - j) > 0));
    if (tmp < mj) {
      mj = tmp;
    }
    mM.at(par.l - 1 - j, 0) = mj;

    tmp = max(find(w.row(j) > 0));
    if (tmp > Mj) {
      Mj = tmp;
    }
    mM.at(j, 1) = Mj;
  }
  par.mM = mM;

  // PP
  arma::imat PP(par.l, par.m, arma::fill::zeros);
  for (int j = 0; j < par.l; j++) {
    PP.row(j).subvec(mM.at(j, 0), mM.at(j, 1)).fill(1);
  }
  par.PP = PP;

  // Return
  return par;
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
List prepare_data_C(arma::vec& X, arma::vec& Y, arma::vec& W) {
  // Declare variables
  par par = prepare_data_par(X, Y, W);

  // Return
  return List::create(Named("l") = par.l,
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
