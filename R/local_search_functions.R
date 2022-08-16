#' Function to take a simple step
#'
#' @param theta Log-parameter
#' @param Psi Proposal
#' @param delta Delta
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param w Sample weights
#' @param PP Reduced index space
#'
#' @return Updated theta parameter
#' @export
simple.step <- function(theta, Psi, delta, ell, m, n, w, PP) {
  rho <- delta
  f.old <- f.theta(theta, n, w, PP)
  f.new <- f.theta(Psi, n, w, PP)
  while (f.new > f.old) {
    Psi <- (theta + Psi) / 2
    rho <- rho / 2
    f.new <- f.theta(Psi, n, w, PP)
  }
  c0 <- rho - f.old + f.new
  if (c0 <= 0) {
    stop("c0 is not strictly positive!")
  }
  t.star <- min(1, rho / (2 * c0))
  theta.new <- matrix(-Inf, nrow = ell, ncol = m)
  theta.new[PP] <- (1 - t.star) * theta[PP] + t.star * Psi[PP]
  return(theta.new)
}

#' Local search (row)
#'
#' @param theta Log-parameter
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param lL (l_k,L_k) index pairs
#' @param PP Reduced index space
#' @param w Sample weights
#' @param w_ul Cumulative row-sums of sample weights
#'
#' @return New proposal Psi and step-size delta
#' @export
local.search1 <- function(theta, ell, m, n, mM, lL, PP, w, w_ul) {
  tmp <- vgamma.tilde1(theta, ell, m, n, mM, w_ul)
  v.tilde <- tmp$v
  gamma.tilde <- tmp$gamma

  lambda.star <- matrix(0, nrow = ell, ncol = m)
  lambda.star[cbind(1:ell, mM[, 1])] <- gamma.tilde[cbind(1:ell, mM[, 1])]
  if (m >= 2) {
    for (k in 2:m) {
      jj <- lL[k, 1]:lL[k - 1, 2]
      if (length(jj) >= 1) {
        lambda.star[jj, k] <- Iso::pava(gamma.tilde[jj, k], v.tilde[jj, k])
      }
    }
  }
  Psi <- lambda1.to.theta(lambda.star, ell, m, mM)
  delta <- sum((-w[PP] + n * exp(theta[PP])) * (theta[PP] - Psi[PP]))
  return(list(Psi = Psi, delta = delta))
}

#' Local search (column)
#'
#' @param theta Log-parameter
#' @param ell Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param lL (l_k,L_k) index pairs
#' @param PP Reduced index space
#' @param w Sample weights
#' @param w_ol Cumulative column-sums of sample weights
#'
#' @return New proposal Psi and step-size delta
#' @export
local.search2 <- function(theta, ell, m, n, mM, lL, PP, w, w_ol) {
  tmp <- vgamma.tilde2(theta, ell, m, n, lL, w_ol)
  v.tilde <- tmp$v
  gamma.tilde <- tmp$gamma
  lambda.star <- matrix(0, nrow = ell, ncol = m)

  lambda.star[cbind(lL[, 1], 1:m)] <- gamma.tilde[cbind(lL[, 1], 1:m)]
  if (ell >= 2) {
    for (j in 2:ell) {
      kk <- mM[j, 1]:mM[j - 1, 2]
      if (length(kk) >= 1) {
        lambda.star[j, kk] <- Iso::pava(gamma.tilde[j, kk], v.tilde[j, kk])
      }
    }
  }
  Psi <- lambda2.to.theta(lambda.star, ell, m, lL)
  delta <- sum((-w[PP] + n * exp(theta[PP])) * (theta[PP] - Psi[PP]))
  return(list(Psi = Psi, delta = delta))
}
