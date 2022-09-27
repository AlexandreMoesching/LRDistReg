## usethis namespace: start
#' @useDynLib LRDistReg, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' TP2 fit function, R version
#'
#' @param par Model parameters
#' @param delta0 Threshold
#'
#' @return h matrix and estimation time
#' @export
TP2_fit_R <- function(par, delta0 = 1e-8) {
  # Extract parameters
  l <- par$l
  lL <- par$lL
  m <- par$m
  mM <- par$mM
  w <- par$w
  w_jplus <- par$w_jplus
  w_plusk <- par$w_plusk
  w_ul <- par$w_ul
  w_ol <- par$w_ol
  n <- par$n
  PP <- par$PP

  # Initialize
  s <- 0
  theta <- matrix(-Inf, nrow = l, ncol = m)
  theta[PP] <- -log(sum(PP))

  # Calibrate
  theta <- calibrate_R(theta, n, w, w_jplus, w_plusk, PP, prec = delta0)

  # New candidate
  tmp <- local_search1_R(theta, l, m, n, mM, lL, PP, w, w_ul)
  s <- s + 1

  # Main while-loop
  while (tmp$delta > delta0) {
    # Real improvement
    theta <- simple_step_R(theta, tmp$Psi, tmp$delta, l, m, n, w, PP)

    # Calibrate
    theta <- calibrate_R(theta, n, w, w_jplus, w_plusk, PP, prec = delta0)

    # New candidate
    if (s %% 2 == 0) {
      tmp <- local_search1_R(theta, l, m, n, mM, lL, PP, w, w_ul)
    } else {
      tmp <- local_search2_R(theta, l, m, n, mM, lL, PP, w, w_ol)
    }

    # Change parity
    s <- s + 1
  }

  # Return probability weights
  h_TP2 <- exp(theta)

  # Compute q
  q_LR <- h_TP2 / rowSums(h_TP2)

  # Compute CDF
  CDF_LR <- t(apply(q_LR, 1, cumsum))

  # Return
  return(list(theta = theta,
              h_TP2 = h_TP2,
              q_LR = q_LR,
              CDF_LR = CDF_LR,
              delta = tmp$delta))
}

#' Isotonic distributional regression (LR, ST, EMP), R version
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified sample weights
#' @param delta0 Threshold
#' @param x0 Set of covariates on which to estimate the distributions
#' @param ST Boolean indicating whether or not the classical isotonic
#' distributional regression will also be computed
#'
#' @return A list of results which depends on the option chosen
#' @export
#'
#' @examples # To be done
dist_reg_R <- function(X, Y, W = rep(1, length(X)),
                     delta0 = 1e-8, x0 = NULL, ST = FALSE) {
  # Compute model parameters
  par <- prepare_data_R(X, Y, W)

  # Fit under LR constraint
  res <- TP2_fit_R(par, delta0)
  res$par <- par

  # Interpolate to x0 if x0 != NULL
  if (!is.null(x0)) {
    res$CDF_LR <- interpolate_R(x0, par$x, res$CDF_LR)
  }
  res$x0 <- x0

  # Fit under usual stochastic ordering constraint
  if (ST) {
    CDF_EMP <- matrix(0, par$l, par$m)
    for (j in 1:par$l) {
      CDF_EMP[j, ] <- cumsum(par$w[j, ]) / par$w_jplus[j]
    }
    res$CDF_ST <- apply(
      CDF_EMP, 2,
      function(z) Iso::pava(z, w = par$w_jplus, decreasing = TRUE)
    )
    res$CDF_EMP <- CDF_EMP

    # Interpolate to x0 if x0 != NULL
    if (any(!is.null(x0))) {
      res$CDF_ST <- interpolate_R(x0, par$x, res$CDF_ST)
      res$CDF_EMP <- interpolate_R(x0, par$x, res$CDF_EMP)
    }
  }

  # Return results
  return(res)
}
