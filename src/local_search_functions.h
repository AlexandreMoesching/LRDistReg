#include "likelihood_functions.h"

void antitonic_reg_ref_cpp(arma::vec& x, arma::vec& gamma, arma::vec& v);
void simple_step_ref_cpp(  arma::mat& theta, arma::mat& Psi, double& delta, const par& par);
void local_search1_ref_cpp(arma::mat& theta, arma::mat& Psi, double& delta, const par& par);
void local_search2_ref_cpp(arma::mat& theta, arma::mat& Psi, double& delta, const par& par);

arma::vec antitonic_reg_cpp(arma::vec gamma, arma::vec v, int size);
arma::mat simple_step_cpp(arma::mat theta, arma::mat Psi, double delta,
                          int ell, int m, arma::imat mM, int n, int w);
List local_search1_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM, arma::imat PP,
                       arma::mat w, arma::mat w_ul);
List local_search2_cpp(arma::mat theta, int ell, int m, int n,
                       arma::imat lL, arma::imat mM, arma::imat PP,
                       arma::mat w, arma::mat w_ol);
