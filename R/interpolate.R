#' Linear interpolation of CDF's, R version
#'
#' @param x0 Set of covariates on which to extend CDF
#' @param x Set of covariates of CDF
#' @param CDF Step-function matrix, its j-th row contains the CDF for X = x_j
#'
#' @return Linear interpolation of CDF's on the new set x0
#' @export
interpolate_R <- function(x0, x, CDF) {
  return(apply(CDF, 2, function(z) stats::approx(x, z, xout = x0, rule = 2)$y))
}
