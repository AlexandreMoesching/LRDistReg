#' Transforms lambda (row) into theta
#'
#' @param lambda Row-wise differences
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param mM (m_j,M_j) index pairs
#'
#' @return Transformed parameter
#' @export
#'
#' @examples
lambda1.to.theta <- function(lambda, ell, m, mM) {
  theta <- matrix(-Inf, nrow = ell, ncol = m)
  for (j in 1:ell) {
    kk <- mM[j, 1]:mM[j, 2]
    theta[j, kk] <- cumsum(lambda[j, kk])
  }
  return(theta)
}

#' Transforms lambda (column) into theta
#'
#' @param lambda Column-wise differences
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param lL (l_k,L_k) index pairs
#'
#' @return Transformed parameter
#' @export
#'
#' @examples
lambda2.to.theta <- function(lambda, ell, m, lL) {
  theta <- matrix(-Inf, nrow = ell, ncol = m)
  for (k in 1:m) {
    jj <- lL[k, 1]:lL[k, 2]
    theta[jj, k] <- cumsum(lambda[jj, k])
  }
  return(theta)
}
