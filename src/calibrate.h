#include "likelihood_functions.h"

void calibrate1_ref_cpp(arma::mat& theta, const par& par);
void calibrate2_ref_cpp(arma::mat& theta, const par& par);
void  calibrate_ref_cpp(arma::mat& theta, const par& par, double& prec);

arma::mat calibrate1_cpp(arma::mat theta, int n, arma::vec w_jplus);
arma::mat calibrate2_cpp(arma::mat theta, int n, arma::vec w_plusk);
arma::mat  calibrate_cpp(arma::mat theta, int ell, arma::imat mM, int n,
                         arma::mat w, arma::vec w_jplus, arma::vec w_plusk,
                         double prec);
