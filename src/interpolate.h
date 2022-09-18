#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp17")]]

arma::mat interpolate_cpp(arma::vec& x0, arma::vec& x, arma::mat& CDF);
