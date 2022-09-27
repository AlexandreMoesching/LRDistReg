#' Approximate calibration of log-parameter, R version
#'
#' @description Computes an optimal additive adjustment of theta, i.e. a matrix
#' theta_new with values
#'     theta_new_jk = theta_jk + x_j + y_k
#' such that f(theta_new) is minimal.
#'
#' @param theta Log-parameter
#' @param n Sample size
#' @param w Sample weights
#' @param w_jplus Row-sums of sample weights
#' @param w_plusk Column sums of sample weights
#' @param PP Reduced index space
#' @param prec Precision for the approximate calibration
#'
#' @return Approximately calibrated log-parameter
#' @export
calibrate_R <- function(theta, n, w, w_jplus, w_plusk, PP, prec = 1e-10) {
  # Initialize f(theta) and delta
  f_old <- ftheta_R(theta, n, w, PP)
  delta <- Inf

  # Loop until improvement is too small
  while (delta > prec) {
    # Row-wise calibration
    theta <- theta - log(rowSums(exp(theta))) + log(w_jplus / n)

    # Column-wise calibration
    theta <- t(t(theta) - log(colSums(exp(theta))) + log(w_plusk / n))

    # Update f(theta) and delta
    f_new <- ftheta_R(theta, n, w, PP)
    delta <- f_old - f_new
    f_old <- f_new
  }

  # Return parameter
  return(theta)
}
