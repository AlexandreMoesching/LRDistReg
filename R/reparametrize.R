#' Transforms lambda (row) into theta, R version
#'
#' @param lambda Row-wise differences
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param mM (m_j,M_j) index pairs
#'
#' @return Transformed parameter
#' @export
lambda1_to_theta_R <- function(lambda, l, m, mM) {
  theta <- matrix(-Inf, nrow = l, ncol = m)
  for (j in 1:l) {
    kk <- mM[j, 1]:mM[j, 2]
    theta[j, kk] <- cumsum(lambda[j, kk])
  }
  return(theta)
}

#' Transforms lambda (column) into theta, R version
#'
#' @param lambda Column-wise differences
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param lL (l_k,L_k) index pairs
#'
#' @return Transformed parameter
#' @export
lambda2_to_theta_R <- function(lambda, l, m, lL) {
  theta <- matrix(-Inf, nrow = l, ncol = m)
  for (k in 1:m) {
    jj <- lL[k, 1]:lL[k, 2]
    theta[jj, k] <- cumsum(lambda[jj, k])
  }
  return(theta)
}
