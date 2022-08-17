#include "prepare_data.h"

void vgamma_tilde1_ref_cpp(arma::mat& theta, arma::mat& v, arma::mat& gamma, const par& par);
void vgamma_tilde2_ref_cpp(arma::mat& theta, arma::mat& v, arma::mat& gamma, const par& par);

List vgamma_tilde1_cpp(arma::mat theta, int ell, int m, int n, arma::imat mM, arma::mat w_ul);
List vgamma_tilde2_cpp(arma::mat theta, int ell, int m, int n, arma::imat lL, arma::mat w_ol);

double f_theta_ref_cpp(arma::mat& theta, const par& par);
double f_theta_cpp(arma::mat theta, int ell, int n, arma::imat mM, arma::mat w);
