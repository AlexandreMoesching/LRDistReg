#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

struct par {
  arma::vec x;
  arma::vec X;
  arma::vec y;
  arma::vec Y;
  int l;
  arma::imat lL;
  int m;
  arma::imat mM;
  int n;
  arma::imat PP;
  arma::mat w;
  arma::vec w_jplus;
  arma::vec w_plusk;
  arma::mat w_ul;
  arma::mat w_ol;
  arma::vec W;
};

par  prepare_data_par(arma::vec& X, arma::vec& Y, arma::vec& W);
List prepare_data_C(  arma::vec& X, arma::vec& Y, arma::vec& W);
