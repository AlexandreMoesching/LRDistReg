#' v-tilde and gamma-tilde functions (row)
#'
#' @param theta Log-parameter
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param w_cumul.1 Cumulative row-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#' @export
vgamma.tilde1 <- function(theta, ell, m, n, mM, w_cumul.1) {
  v <- matrix(0, nrow = ell, ncol = m)
  gamma <- matrix(0, nrow = ell, ncol = m)
  for (j in 1:ell) {
    kk <- mM[j, 1]:mM[j, 2]
    tmp1 <- n * rev(cumsum(rev(exp(theta[j, kk]))))
    tmp2 <- c(theta[j, mM[j, 1]], diff(theta[j, kk])) # theta^*
    v[j, kk] <- tmp1
    gamma[j, kk] <- tmp2 + w_cumul.1[j, kk] / tmp1 - 1
  }
  return(list(v = v, gamma = gamma))
}

#' v-tilde and gamma-tilde functions (column)
#'
#' @param theta Log-parameter
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param lL (l_k,L_k) index pairs
#' @param w_cumul.2 Cumulative column-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#' @export
vgamma.tilde2 <- function(theta, ell, m, n, lL, w_cumul.2) {
  v <- matrix(0, nrow = ell, ncol = m)
  gamma <- matrix(0, nrow = ell, ncol = m)
  for (k in 1:m) {
    jj <- lL[k, 1]:lL[k, 2]
    tmp1 <- n * rev(cumsum(rev(exp(theta[jj, k]))))
    tmp2 <- c(theta[lL[k, 1], k], diff(theta[jj, k])) # theta^*
    v[jj, k] <- tmp1
    gamma[jj, k] <- tmp2 + w_cumul.2[jj, k] / tmp1 - 1
  }
  return(list(v = v, gamma = gamma))
}

#' Negative log-likelihood in terms of log-parameter
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w Sample weights
#' @param PP Reduced index space
#'
#' @return Negative log-likelihood in terms of log-parameter
#' @export
f.theta <- function(theta, n, w, PP) {
  return(sum(-w[PP] * theta[PP] + n * exp(theta[PP])))
}
