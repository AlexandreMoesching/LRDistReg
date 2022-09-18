#include "calibrate.h"

struct pava_par {
  arma::ivec PP;
  arma::vec WW;
  arma::vec MM;
};

void simple_step_ref_cpp(  arma::mat& theta, arma::mat& Psi,
                           double& delta, const par& par);
void local_search1_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           arma::mat& v, arma::mat& gamma,
                           arma::mat& lambda_star, pava_par& par1,
                           double& delta, const par& par);
void local_search2_ref_cpp(arma::mat& theta, arma::mat& Psi,
                           arma::mat& v, arma::mat& gamma,
                           arma::mat& lambda_star, pava_par& par2,
                           double& delta, const par& par);

arma::mat simple_step_cpp(arma::mat& theta, arma::mat& Psi, double delta,
                          int ell, arma::imat& mM, int n, arma::mat& w);
List local_search1_cpp(arma::mat& theta, int ell, int m, int n,
                       arma::imat& lL, arma::imat& mM,
                       arma::mat& w, arma::mat& w_ul);
List local_search2_cpp(arma::mat& theta, int ell, int m, int n,
                       arma::imat& lL, arma::imat& mM,
                       arma::mat& w, arma::mat& w_ol);
