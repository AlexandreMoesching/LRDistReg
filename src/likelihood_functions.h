#include "reparametrize.h"

void vgamma_tilde1_ref(arma::mat& theta, arma::mat& v, arma::mat& gamma, const par& par);
void vgamma_tilde2_ref(arma::mat& theta, arma::mat& v, arma::mat& gamma, const par& par);

List vgamma_tilde1_C(arma::mat& theta, int l, int m, int n, arma::imat& mM, arma::mat& w_ul);
List vgamma_tilde2_C(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::mat& w_ol);

double f_theta_ref(arma::mat& theta, const par& par);
double f_theta_C(arma::mat& theta, int l, int n, arma::imat& mM, arma::mat& w);
