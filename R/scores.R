#' Computes the 'simple score'
#'
#' @param x0 Extended set of covariates
#' @param l0 Cardinality of extended set of covariates `x0`
#' @param res A list containing `CDF_LR`, `CDF_ST`, `CDF_EMP` and `par`. I.e.,
#' the estimated conditional cumulative distribution functions under
#' likelihood-ratio order and usual stochastic order, the conditional empirical
#' distribution functions, all of which having `x0` for covariates, and a list
#' of parameters, of which the unique observations `y` and their number `m` are
#' relevant
#' @param a,b Shape and scale functions of the true Gamma conditional
#' distribution functions
#' @param rel.tol Relative tolerance for numerical integration
#'
#' @return A l0-by-3 matrix of conditional simple scores
#' @export
SS <- function(x0, l0, res, a, b, rel.tol = 1e-9) {
  y.ext <- c(res$par$y, Inf)
  SS <- matrix(0, nrow = l0, ncol = 3)
  colnames(SS) <- c("LR", "ST", "EMP")
  for (j in 1:l0) {
    # Integral on (-\infty,y_1)
    fun <- function(z) {
      stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])) *
      stats::dgamma(z, shape = a(x0[j]), scale = b(x0[j]))
    }
    SS[j, ] <-
      stats::integrate(fun, -Inf, res$par$y[1], rel.tol = rel.tol)$value

    # Integral on [y_k,y_k+1)
    for (k in 1:res$par$m) {
      # LR
      fun <- function(z) {
        abs(res$CDF_LR[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j]))) *
          stats::dgamma(z, shape = a(x0[j]), scale = b(x0[j]))
      }
      SS[j, 1] <- SS[j, 1] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value

      # ST
      fun <- function(z) {
        abs(res$CDF_ST[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j]))) *
          stats::dgamma(z, shape = a(x0[j]), scale = b(x0[j]))
      }
      SS[j, 2] <- SS[j, 2] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value

      # EMP
      fun <- function(z) {
        abs(res$CDF_EMP[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j]))) *
          stats::dgamma(z, shape = a(x0[j]), scale = b(x0[j]))
      }
      SS[j, 3] <- SS[j, 3] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value
    }
  }
  return(SS)
}

#' Computes the 'CRPS'
#'
#' @param x0 Extended set of covariates
#' @param l0 Cardinality of extended set of covariates `x0`
#' @param res A list containing `CDF_LR`, `CDF_ST`, `CDF_EMP` and `par`. I.e.,
#' the estimated conditional cumulative distribution functions under
#' likelihood-ratio order and usual stochastic order, the conditional empirical
#' distribution functions, all of which having `x0` for covariates, and a list
#' of parameters, of which the unique observations `y` and their number `m` are
#' relevant
#' @param a,b Shape and scale functions of the true Gamma conditional
#' distribution functions
#' @param rel.tol Relative tolerance for numerical integration
#'
#' @return A l0-by-3 matrix of conditional CRPS
#' @export
CRPS <- function(x0, l0, res, a, b, rel.tol = 1e-9) {
  y.ext <- c(res$par$y, Inf)
  CRPS <- matrix(b(x0) / beta(1 / 2, a(x0)), nrow = l0, ncol = 3)
  colnames(CRPS) <- c("LR", "ST", "EMP")
  for (j in 1:l0) {
    # Integral on (-\infty,y_1)
    fun <- function(z) {
      (stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])))^2
    }
    CRPS[j, ] <- CRPS[j, ] +
      stats::integrate(fun, -Inf, res$par$y[1], rel.tol = rel.tol)$value

    # Integral on [y_k,y_k+1)
    for (k in 1:res$par$m) {
      # LR
      fun <- function(z) {
        (res$CDF_LR[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])))^2
      }
      CRPS[j, 1] <- CRPS[j, 1] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value

      # ST
      fun <- function(z) {
        (res$CDF_ST[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])))^2
      }
      CRPS[j, 2] <- CRPS[j, 2] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value

      # EMP
      fun <- function(z) {
        (res$CDF_EMP[j, k] -
          stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])))^2
      }
      CRPS[j, 3] <- CRPS[j, 3] +
        stats::integrate(fun, y.ext[k], y.ext[k + 1], rel.tol = rel.tol)$value
    }
  }
  return(CRPS)
}
