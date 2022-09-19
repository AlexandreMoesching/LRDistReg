#include "interpolate.h"

//' Linear interpolation of the cdf's, C++ version
//'
//' @param x0 Set of covariates on which to extend CDF
//' @param x Set of covariates of CDF
//' @param CDF Step-function matrix, its j-th row contains the cdf for X = x_j
//'
//' @return Linear interpolation of the cdf's on the new set x0
//'
//' @export
//[[Rcpp::export]]
arma::mat interpolate_cpp(arma::vec& x0, arma::vec& x, arma::mat& CDF) {
  // Declare variables
  int l0 = x0.n_elem, l = x.n_elem, m = CDF.n_cols;
  arma::mat CDF0(l0, m);
  double x_min = min(x), x_max = max(x);

  // Determine r and s such that min(x) <= x0[r] and x0[s] <= max(x)
  int r = 0;
  int s = l0 - 1;
  while (x0[r] < x_min) {
    CDF0.row(r) = CDF.row(0);
    r++;
  }
  while (x_max < x0[s]) {
    CDF0.row(s) = CDF.row(l - 1);
    s--;
  }

  // Interpolate if number of points of x0 in [min(x), max(x)] is positive
  if (s - r + 1 > 0) {
    arma::vec x_int = x0.subvec(r, s);
    arma::vec y_int(s - r + 1);
    for (int k = 0; k < m; k++) {
      interp1(x, CDF.col(k), x_int, y_int, "*linear");
      CDF0.col(k).subvec(r, s) = y_int;
    }
  }

  // Return
  return(CDF0);
}
