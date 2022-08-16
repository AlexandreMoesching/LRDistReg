#include "prepare_data.h"

void lambda1_to_theta_void_cpp(arma::mat& lambda, arma::mat& theta, const par& par);
void lambda2_to_theta_void_cpp(arma::mat& lambda, arma::mat& theta, const par& par);

arma::mat lambda1_to_theta_cpp(arma::mat lambda, int ell, int m, arma::imat mM);
arma::mat lambda2_to_theta_cpp(arma::mat lambda, int ell, int m, arma::imat lL);
