#include "local_search_functions.h"
#include "interpolate.h"

void TP2_fit_ref(arma::mat& h_TP2, arma::mat& q_LR, arma::mat& CDF_LR,
                 arma::mat& theta, arma::mat& Psi,
                 arma::mat& v, arma::mat& gamma, arma::mat& lambda_star,
                 pava_par& par1, pava_par& par2,
                 double& delta, double delta0, par& par);

void ST_fit_ref(arma::mat& CDF_EMP, arma::mat& CDF_ST,
                pava_par& par3, par& par);

List dist_reg_C(arma::vec& X, arma::vec& Y, arma::vec& W,
                double delta0, arma::vec x0, bool ST);

