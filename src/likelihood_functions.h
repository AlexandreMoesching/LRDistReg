#include "reparametrize.h"

void vg1_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par);
void vg2_ref(arma::mat& theta, arma::mat& v, arma::mat& g, const par& par);

List vg1_C(arma::mat& theta, int l, int m, int n, arma::imat& mM, arma::mat& w_ul);
List vg2_C(arma::mat& theta, int l, int m, int n, arma::imat& lL, arma::mat& w_ol);

double ftheta_ref(arma::mat& theta, const par& par);
double ftheta_C(arma::mat& theta, int l, int n, arma::imat& mM, arma::mat& w);
