#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' rcpparma_hello_world
//'
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);

    return m1 + 3 * (m1 + m2);
}

//' rcpparma_outerproduct
//' @param x A vector
//'
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

//' rcpparma_innerproduct
//' @param x A vector
//'
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


//' rcpparma_bothproducts
//' @param x A vector
//'
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}
