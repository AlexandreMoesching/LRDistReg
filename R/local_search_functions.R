#' Function to take a simple step, R version
#'
#' @param theta Log-parameter
#' @param Psi Proposal
#' @param delta Delta
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param w Sample weights
#' @param PP Reduced index space
#'
#' @return Updated theta parameter
#' @export
simple_step_R <- function(theta, Psi, delta, l, m, n, w, PP) {
  rho <- delta
  f_old <- ftheta_R(theta, n, w, PP)
  f_new <- ftheta_R(Psi, n, w, PP)
  while (f_new > f_old) {
    Psi <- (theta + Psi) / 2
    rho <- rho / 2
    f_new <- ftheta_R(Psi, n, w, PP)
  }
  t_star <- min(1, rho / (2 * (rho - f_old + f_new)))
  theta_new <- matrix(-Inf, nrow = l, ncol = m)
  theta_new[PP] <- (1 - t_star) * theta[PP] + t_star * Psi[PP]
  return(theta_new)
}

#' Local search (row), R version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
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
local_search1_R <- function(theta, l, m, n, mM, lL, PP, w, w_ul) {
  tmp <- vg_tilde1_R(theta, l, m, n, mM, w_ul)
  v_tilde <- tmp$v
  g_tilde <- tmp$g
  lambda_star <- matrix(0, nrow = l, ncol = m)
  lambda_star[cbind(1:l, mM[, 1])] <- g_tilde[cbind(1:l, mM[, 1])]
  if (m >= 2) {
    for (k in 2:m) {
      jj <- lL[k, 1]:lL[k - 1, 2]
      if (length(jj) >= 1) {
        lambda_star[jj, k] <- Iso::pava(g_tilde[jj, k], v_tilde[jj, k])
      }
    }
  }
  Psi <- lambda1_to_theta_R(lambda_star, l, m, mM)
  # We take max(0,-) to avoid numerical instability in case theta ~ Psi
  delta <- max(0, sum((-w[PP] + n * exp(theta[PP])) * (theta[PP] - Psi[PP])))
  return(list(Psi = Psi, delta = delta))
}

#' Local search (column), R version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
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
local_search2_R <- function(theta, l, m, n, mM, lL, PP, w, w_ol) {
  tmp <- vg_tilde2_R(theta, l, m, n, lL, w_ol)
  v_tilde <- tmp$v
  g_tilde <- tmp$g
  lambda_star <- matrix(0, nrow = l, ncol = m)

  lambda_star[cbind(lL[, 1], 1:m)] <- g_tilde[cbind(lL[, 1], 1:m)]
  if (l >= 2) {
    for (j in 2:l) {
      kk <- mM[j, 1]:mM[j - 1, 2]
      if (length(kk) >= 1) {
        lambda_star[j, kk] <- Iso::pava(g_tilde[j, kk], v_tilde[j, kk])
      }
    }
  }
  Psi <- lambda2_to_theta_R(lambda_star, l, m, lL)
  # We take max(0,-) to avoid numerical instability in case theta ~ Psi
  delta <- max(0, sum((-w[PP] + n * exp(theta[PP])) * (theta[PP] - Psi[PP])))
  return(list(Psi = Psi, delta = delta))
}
