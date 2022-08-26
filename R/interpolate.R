#' Linear interpolation of the cdf's
#'
#' @param x0 Set of covariates on which to extend CDF.MAT
#' @param x Set of covariates of CDF.MAT, {X_1, X_2, ..., X_N}
#' @param CDF.MAT Step-function matrix, its j-th row contains the cdf for X = x_j
#'
#' @return Linear interpolation of the cdf's on the new set x0
#' @export
interpolate <- function(x0, x, CDF.MAT) {
  my.approx <- function(z) stats::approx(x, z, xout = x0, rule = 2)$y
  res <- apply(CDF.MAT, 2, my.approx)
  return(res)
}
