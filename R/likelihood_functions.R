#' v-tilde and gamma-tilde functions (row), R version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param w_ul Cumulative row-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#' @export
vg1_R <- function(theta, l, m, n, mM, w_ul) {
  v <- g <- matrix(0, nrow = l, ncol = m)
  for (j in 1:l) {
    kk <- mM[j, 1]:mM[j, 2]
    tmp1 <- n * rev(cumsum(rev(exp(theta[j, kk]))))
    tmp2 <- c(theta[j, mM[j, 1]], diff(theta[j, kk])) # theta^*
    v[j, kk] <- tmp1
    g[j, kk] <- tmp2 + w_ul[j, kk] / tmp1 - 1
  }
  return(list(v = v, g = g))
}

#' v-tilde and gamma-tilde functions (column), R version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param lL (l_k,L_k) index pairs
#' @param w_ol Cumulative column-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#' @export
vg2_R <- function(theta, l, m, n, lL, w_ol) {
  v <- g <- matrix(0, nrow = l, ncol = m)
  for (k in 1:m) {
    jj <- lL[k, 1]:lL[k, 2]
    tmp1 <- n * rev(cumsum(rev(exp(theta[jj, k]))))
    tmp2 <- c(theta[lL[k, 1], k], diff(theta[jj, k])) # theta^*
    v[jj, k] <- tmp1
    g[jj, k] <- tmp2 + w_ol[jj, k] / tmp1 - 1
  }
  return(list(v = v, g = g))
}

#' Negative log-likelihood in terms of log-parameter, R version
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w Sample weights
#' @param PP Reduced index space
#'
#' @return Negative log-likelihood in terms of log-parameter
#' @export
ftheta_R <- function(theta, n, w, PP) {
  return(sum(-w[PP] * theta[PP] + n * exp(theta[PP])))
}
