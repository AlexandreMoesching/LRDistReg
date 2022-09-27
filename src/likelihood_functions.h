#include "reparametrize.h"

void vg_tilde1_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par);
void vg_tilde2_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par);

List vg_tilde1_C(arma::mat& theta, int l, int m, int n, arma::imat& mM, arma::mat& w_ul);
List vg_tilde2_C(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::mat& w_ol);

double ftheta_ref(arma::mat& theta, const par& par);
double ftheta_C(arma::mat& theta, int l, int n, arma::imat& mM, arma::mat& w);
