## usethis namespace: start
#' @useDynLib LRDistReg, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' TP2 fit function
#'
#' @param par Model parameters
#' @param delta0 Threshold
#'
#' @return h matrix and estimation time
#' @export
TP2.fit <- function(par, delta0 = 1e-8) {
  # Start timer
  start.time <- Sys.time()

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
  theta <- matrix(-Inf, nrow = l, ncol = m)
  theta[PP] <- -log(sum(PP))
  delta <- Inf

  # Main while-loop
  s <- 0
  while (delta > delta0) {
    # Calibrate
    theta <- calibrate(theta, n, w, w_jplus, w_plusk, PP, prec = 1e-10)

    # New candidate
    if (s %% 2 == 0) {
      tmp <- local.search1(theta, l, m, n, mM, lL, PP, w, w_ul)
    } else {
      tmp <- local.search2(theta, l, m, n, mM, lL, PP, w, w_ol)
    }

    Psi <- tmp$Psi
    delta <- tmp$delta

    # Real improvement
    theta <- simple.step(theta, Psi, delta, l, m, n, w, PP)

    # Change parity
    s <- s + 1
  }

  # Return probability weights
  h_TP2 <- exp(theta)

  # End timer
  tot.time <- Sys.time() - start.time

  # Compute q
  q_LR <- h_TP2 / rowSums(h_TP2)

  # Compute CDF
  CDF_LR <- t(apply(q_LR, 1, cumsum))

  # Return
  return(list(h_TP2 = h_TP2, q_LR = q_LR, CDF_LR = CDF_LR, tot.time = tot.time))
}

#' Isotonic distributional regression (LR, ST, EMP)
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified sample weights
#' @param show.design Whether or not to plot the design
#' @param indices Whether or not to plot indices instead of true values
#' @param suggest.delta0 Whether or not to suggest the threshold delta0
#' @param delta0 Threshold
#' @param x0 Set of covariates on which to estimate the distributions
#' @param ST Boolean indicating whether or not the classical isotonic
#' distributional regression will also be computed
#'
#' @return A list of results which depends on the option chosen
#' @export
#'
#' @examples # To be done
dist.reg <- function(X, Y, W = rep(1, length(X)),
                     show.design = FALSE, indices = FALSE,
                     suggest.delta0 = FALSE,
                     delta0 = 1e-8, x0 = NULL, ST = FALSE) {
  # Compute model parameters
  par <- prepare.data(X, Y, W)

  # Plot design
  if (show.design == TRUE) {
    plotD(par, indices = indices)
  }

  # Suggest delta0
  if (suggest.delta0 == TRUE) {
    theta <- matrix(-Inf, nrow = par$l, ncol = par$m)
    theta[par$PP] <- -log(sum(par$PP))
    calibrate(
      theta, par$n, par$w, par$w_jplus, par$w_plusk, par$PP
    )
    f.tmp <- f.theta(theta, par$n, par$w, par$PP)
    delta0 <- 10^floor(log10(f.tmp) * 1.1 - 3)
  }

  # Fit under LR constraint
  res <- TP2.fit(par, delta0)

  # Interpolate to x0 if x0 != NULL
  if (!is.null(x0)) {
    res$CDF_LR <- interpolate(x0, par$x, res$CDF_LR)
  }

  # Fit under usual stochastic ordering constraint
  if (ST) {
    CDF_EMP <- matrix(0, par$l, par$m)
    for (j in 1:par$l) {
      CDF_EMP[j, ] <- cumsum(par$w[j, ]) / par$w_jplus[j]
    }
    res$CDF_EMP <- CDF_EMP
    res$CDF_ST <- apply(
      CDF_EMP, 2,
      function(z) Iso::pava(z, w = par$w_jplus, decreasing = TRUE)
    )

    # Interpolate to x0 if x0 != NULL
    if (any(!is.null(x0))) {
      res$CDF_EMP <- interpolate(x0, par$x, res$CDF_EMP)
      res$CDF_ST <- interpolate(x0, par$x, res$CDF_ST)
    }
  }

  # Return results
  res$par <- par
  res$delta0 <- delta0
  res$x0 <- x0

  return(res)
}
