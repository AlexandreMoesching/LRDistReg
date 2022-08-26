#' Plot the design
#'
#' @param par Model parameters
#' @param indices Whether or not to plot indices instead of true values
#'
#' @return Plot
#' @export
plotD <- function(par, indices = FALSE) {
  ell <- par$ell
  m <- par$m
  P <- which(par$PP, arr.ind = TRUE)

  if (indices == TRUE) {
    D <- which(par$w >= 1, arr.ind = TRUE)
    D0 <- D
    x <- 1:ell
    y <- 1:m
    G <- P
    graphics::plot(D,
      type = "n", mgp = c(1, 1, 0),
      xaxt = "n", yaxt = "n",
      xlab = expression(italic(j)), ylab = expression(italic(k))
    )
  } else {
    D <- cbind(par$X, par$Y)
    D0 <- which(par$w >= 1, arr.ind = TRUE)
    x <- par$x
    y <- par$y
    G <- cbind(x[P[, 1]], y[P[, 2]])
    graphics::plot(D,
      type = "n", mgp = c(1, 1, 0),
      xaxt = "n", yaxt = "n",
      xlab = expression(italic(x)), ylab = expression(italic(y))
    )
  }

  # Vertical segments
  mM <- par$mM
  for (j in 1:ell) {
    graphics::segments(
      x0 = x[j], y0 = y[mM[j, 1]], x1 = x[j], y1 = y[mM[j, 2]],
      col = "lightgrey"
    )
  }

  # Horizontal segments
  lL <- par$lL
  for (k in 1:m) {
    graphics::segments(
      x0 = x[lL[k, 1]], y0 = y[k], x1 = x[lL[k, 2]], y1 = y[k],
      col = "lightgrey"
    )
  }

  # Color
  col <- par$w[D0] + 1

  # Points
  graphics::points(G, pch = 16, cex = 0.5)
  graphics::points(D, pch = 16, cex = 0.8, col = col)
}
