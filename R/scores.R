#' Special transformation
#'
#' @param d,z  Scalars
#' @noRd
rho <- function(d, z) {
  d * z - z^2 / 2
}

#' Computes the 'simple score' and 'CRPS' for the gamma model
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
#' @param rel_tol Relative tolerance for numerical integration
#'
#' @return A list of two l0-by-3 matrix of conditional simple scores and CRPS's
#' @export
SS_CRPS_gamma <- function(x0, l0, res, a, b, rel_tol = 1e-9) {
  y <- c(res$par$y)
  m <- res$par$m

  # Pre-compute matrix of Gamma CDF's
  G0_jk <- outer(x0, y,
                      function(xj, yk) stats::pgamma(yk,
                                                     shape = a(xj),
                                                     scale = b(xj)))
  G1_jk <- outer(x0, y,
                      function(xj, yk) stats::pgamma(yk,
                                                     shape = a(xj) + 1,
                                                     scale = b(xj)))

  # Pre-compute some constants
  # Predefine some constants
  cc <- b(x0) * gamma(a(x0) + 1) / gamma(a(x0))

  # The three columns of SS are identical in the beginning. They correspond to
  # the sum of:
  # - The integral of      G_xj  g_xj on (0,   y_1)
  # - The integral of (1 - G_xj) g_xj on (y_m, Inf)
  SS <- matrix(
    (1/2) * (1 + G0_jk[, 1]^2 + G0_jk[, m]^2) - G0_jk[, m],
    nrow = l0, ncol = 3
  )
  colnames(SS) <- c("LR", "ST", "EMP")

  # The three columns of CRPS are identical in the beginning. They correspond to
  # the last term in the expression of the CRPS involving the Beta function
  CRPS <- matrix(b(x0) / beta(1 / 2, a(x0)), nrow = l0, ncol = 3)
  colnames(CRPS) <- c("LR", "ST", "EMP")

  if (m > 1) {
    # For each x_j in x0
    for (j in 1:l0) {
      # First term corresponds to the sum of:
      # - The integral of      G_xj ^2 on (0,   y_m)
      tmp1 <- stats::integrate(f = function(z)
        ( stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])) )^2,
        lower = 0,
        upper = y[m],
        subdivisions = 1e4,
        rel.tol = rel_tol
      )$value
      # - The integral of (1 - G_xj)^2 on (y_m, Inf)
      tmp2 <- stats::integrate(f = function(z)
        ( 1 - stats::pgamma(z, shape = a(x0[j]), scale = b(x0[j])) )^2,
        lower = y[m],
        upper = Inf,
        subdivisions = 1e4,
        rel.tol = rel_tol
      )$value

      # Add the results to the current CRPS matrix
      CRPS[j, ] <- CRPS[j, ] + tmp1 + tmp2

      # Compute remaining terms
      for (k in 1:(m - 1)) { # (y_k, y_k+1)
        # For each orders
        for (ord in 1:3) {
          # Define appropriate \tilde{G}_jk
          if (ord == 1) {
            tG_jk <- res$CDF_LR[j, k]
          } else if (ord == 2) {
            tG_jk <- res$CDF_ST[j, k]
          } else {
            tG_jk <- res$CDF_EMP[j, k]
          }

          # Computations related to SS
          if (tG_jk >= G0_jk[j, k + 1]) {
            SS[j, ord] <- SS[j, ord] +
              rho(tG_jk, G0_jk[j, k + 1]) - rho(tG_jk, G0_jk[j, k])
          } else if (tG_jk <= G0_jk[j, k]) {
            SS[j, ord] <- SS[j, ord] +
              rho(tG_jk, G0_jk[j, k]) - rho(tG_jk, G0_jk[j, k + 1])
          } else {
            SS[j, ord] <- SS[j, ord] +
              tG_jk^2 - rho(tG_jk, G0_jk[j, k]) - rho(tG_jk, G0_jk[j, k + 1])
          }

          # Computations related to CRPS
          CRPS[j, ord] <- CRPS[j, ord] +
            tG_jk^2 * (y[k + 1] - y[k]) - 2 * tG_jk * (
              y[k + 1] * G0_jk[j, k + 1] -
                  y[k] * G0_jk[j, k] -
                cc[j] * (G1_jk[j, k + 1] - G1_jk[j, k])
            )
        }
      }
    }
  }
  return(list(SS = SS, CRPS = CRPS))
}
