## usethis namespace: start
#' @useDynLib LRDistReg, .registration = TRUE
## usethis namespace: end
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

strsplit1 <- function(x, split) {
  strsplit(x, split = split)[[1]]
}
