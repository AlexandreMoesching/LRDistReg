## usethis namespace: start
#' @useDynLib LRDistReg, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

#' Isotonic distributional regression (LR, US, EMP)
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified sample weights
#' @param show.design Whether or not to plot the design
#' @param indices Whether or not to plot indices instead of true values
#' @param suggest.delta0 Whether or not to suggest the threshold delta0
#' @param delta0 Threshhold
#' @param echo Whether or not to print intermediate steps
#' @param out.file Whether or not to wwite intermediate steps in a separate file
#' @param IDR Whether or not to compute Isotonic Distributional Regression
#' @param x0 Set of covariates on which to estimate the distributions
#'
#' @return A list of results which depends on the option chosen
#' @export
#'
#' @examples #To be done
dist.reg <- function(X, Y, W = rep(1, length(X)),
                     show.design = FALSE, indices = FALSE,
                     suggest.delta0 = FALSE, delta0 = 1e-2,
                     echo = FALSE, out.file = FALSE,
                     IDR = FALSE, x0 = NULL) {
  # Setup of the design
  D.Setup <- prepare.data(X, Y, W)

  # Plot design
  if (show.design == TRUE) {
    plotD(D.Setup, indices = indices)
  }

  # Suggest delta0
  if (suggest.delta0 == TRUE) {
    theta <- matrix(-Inf, nrow = D.Setup$ell, ncol = D.Setup$m)
    theta[D.Setup$PP] <- -log(sum(D.Setup$PP))
    theta <- calibrate1(theta, D.Setup$n, D.Setup$w_jplus)
    theta <- calibrate2(theta, D.Setup$n, D.Setup$w_plusk)
    f.tmp <- f.theta(theta, D.Setup$n, D.Setup$w, D.Setup$PP)
    delta0 <- 10^floor(log10(f.tmp) * 1.1 - 3)
  }

  # Fit under LR constraint
  res <- TP2.fit(D.Setup, delta0, echo, out.file)
  res$q.LR <- res$h.TP2 / rowSums(res$h.TP2)
  CDF.LR <- t(apply(res$q.LR, 1, cumsum))

  # Interpolate to x0 if x0 != NULL
  if (!is.null(x0)) {
    CDF.LR <- interpolate(x0, D.Setup$x, CDF.LR)
  }
  res$CDF.LR <- CDF.LR

  # Fit under usual stochastic ordering constraint
  if (IDR == TRUE) {
    CDF.EMP <- matrix(0, D.Setup$ell, D.Setup$m)
    for (j in 1:D.Setup$ell) {
      CDF.EMP.fun <- stats::ecdf(Y[X == D.Setup$x[j]])
      CDF.EMP[j, ] <- CDF.EMP.fun(D.Setup$y)
    }
    CDF.ST <- apply(CDF.EMP, 2, function(z) Iso::pava(z, decreasing = TRUE))

    # Interpolate to x0 if x0 != NULL
    if (any(!is.null(x0))) {
      CDF.EMP <- interpolate(x0, D.Setup$x, CDF.EMP)
      CDF.ST <- interpolate(x0, D.Setup$x, CDF.ST)
    }
    res$CDF.EMP <- CDF.EMP
    res$CDF.ST <- CDF.ST
  }

  # Return results
  res$D.Setup <- D.Setup
  res$delta0 <- delta0
  res$x0 <- x0

  return(res)
}

#' TP2 fit function
#'
#' @param D.Setup Parameters
#' @param delta0 Threshhold
#' @param echo Whether or not to print intermediate steps
#' @param out.file Whether or not to wwite intermediate steps in a separate file
#'
#' @return h matrix and estimation time
#' @export
TP2.fit <- function(D.Setup, delta0 = 1e-1, echo = FALSE, out.file = FALSE) {
  # Start timer
  start.time <- Sys.time()

  # Extract parameters
  ell <- D.Setup$ell
  lL <- D.Setup$lL
  m <- D.Setup$m
  mM <- D.Setup$mM
  w <- D.Setup$w
  w_jplus <- D.Setup$w_jplus
  w_plusk <- D.Setup$w_plusk
  w_ul <- D.Setup$w_ul
  w_ol <- D.Setup$w_ol
  n <- D.Setup$n
  PP <- D.Setup$PP

  # Initialize
  theta <- matrix(-Inf, nrow = ell, ncol = m)
  theta[PP] <- -log(sum(PP))
  delta <- Inf

  # Start writing an output file
  s <- 0
  if (out.file == TRUE) {
    sink("outfile.txt")
  }
  if (echo == TRUE | out.file == TRUE) {
    cat("iteration ", s,
      ", delta = ", delta,
      ", f(theta) = ", f.theta(theta, n, w, PP), "\n",
      sep = ""
    )
  }

  # Main while-loop
  while (delta > delta0) {
    # Calibrate
    theta <- calibrate(theta, n, w, w_jplus, w_plusk, PP, prec = 1e-5)

    if (s %% 2 == 0) {
      # Old calibration
      # theta <- calibrate1(theta, n, w_jplus)

      # New candidate
      tmp <- local.search1(theta, ell, m, n, mM, lL, PP, w, w_ul)
    } else {
      # Old calibration
      # theta <- calibrate2(theta, n, w_plusk)

      # New candidate
      tmp <- local.search2(theta, ell, m, n, mM, lL, PP, w, w_ol)
    }

    Psi <- tmp$Psi
    delta <- tmp$delta

    # Real improvement
    theta <- simple.step(theta, Psi, delta, ell, m, n, w, PP)

    # Change parity
    s <- s + 1

    # Write into a file
    if (echo == TRUE | out.file == TRUE) {
      cat("iteration ", s,
        ", delta = ", delta,
        ", f(theta) = ", f.theta(theta, n, w, PP), "\n",
        sep = ""
      )
    }
  }

  # End writing an output file
  if (out.file == TRUE) {
    sink()
  }

  # Return probability weights
  h.TP2 <- exp(theta)

  # End timer
  tot.time <- Sys.time() - start.time

  return(list(h.TP2 = h.TP2, tot.time = tot.time))
}
