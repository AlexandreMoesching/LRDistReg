#' Row-calibration of log-parameter
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w_j.plus Row-sums of sample weights
#'
#' @return Row-calibrated log-parameter
#' @export
calibrate1 <- function(theta, n, w_j.plus) {
  theta <- theta - log(rowSums(exp(theta))) + log(w_j.plus / n)
  return(theta)
}

#' Column-calibration of log-parameter
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w_plus.k Column-sums of sample weights
#'
#' @return Row-calibrated log-parameter
#' @export
calibrate2 <- function(theta, n, w_plus.k) {
  theta <- t(t(theta) - log(colSums(exp(theta))) + log(w_plus.k / n))
  return(theta)
}

#' Approximate calibration of log-parameter
#'
#' @description Computes an optimal additive adjustment of theta, i.e. a matrix
#' theta.new with values
#'     theta.new_jk = theta_jk + x_j + y_k
#' such that f(theta.new) is minimal.
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w Sample weights
#' @param w_j.plus Row-sums of sample weights
#' @param w_plus.k Column sums of sample weights
#' @param PP Reduced index space
#' @param prec Precision for the approximate calibration
#'
#' @return Approximately calibrated log-parameter
#' @export
calibrate <- function(theta, n, w, w_j.plus, w_plus.k, PP, prec = 1e-7)
{
  ftheta.old <- f.theta(theta, n, w, PP)
  delta <- Inf
  while (delta > prec) {
    # Row-wise calibration
    theta <- theta - log(rowSums(exp(theta))) + log(w_j.plus / n)

    # Column-wise calibration
    theta <- t(t(theta) - log(colSums(exp(theta))) + log(w_plus.k / n))

    ftheta <- f.theta(theta, n, w, PP)
    delta <- ftheta.old - ftheta
    ftheta.old <- ftheta
  }
  return(theta)
}
