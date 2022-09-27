#include "likelihood_functions.h"

void calibrate_ref(arma::mat& theta, const par& par, double& prec);

arma::mat calibrate_C(arma::mat& theta, int l, arma::imat& mM, int n,
                      arma::mat& w, arma::vec& w_jplus, arma::vec& w_plusk,
                      double prec);
