#' Linear interpolation of the cdf's
#'
#' @param x0 Set of covariates on which to extend CDF
#' @param x Set of covariates of CDF, {X_1, X_2, ..., X_N}
#' @param CDF Step-function matrix, its j-th row contains the cdf for X = x_j
#'
#' @return Linear interpolation of the cdf's on the new set x0
#' @export
interpolate <- function(x0, x, CDF) {
  my.approx <- function(z) stats::approx(x, z, xout = x0, rule = 2)$y
  res <- apply(CDF, 2, my.approx)
  return(res)
}
