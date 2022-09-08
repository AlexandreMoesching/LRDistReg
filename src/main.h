#include "local_search_functions.h"

void TP2_fit_ref_cpp(arma::mat& h_TP2, arma::mat& theta, arma::mat& Psi,
                     arma::mat& v, arma::mat& gamma, arma::mat& lambda_star,
                     pava_par& par1, pava_par& par2,
                     double& delta, double delta0, par& par);
List TP2_fit_cpp(arma::vec X, arma::vec Y, arma::vec W, double delta0);
