#' Prepare the data
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified weights
#'
#' @return
#' @export
#'
#' @examples
prepare.data <- function(X, Y, W = rep(1, length(X))) {
  # Number of observation pairs
  n <- length(X)

  # Sorted X/Y
  sX <- sort(X)
  sY <- sort(Y)

  # Which of the sorted X/Y's are unique values
  newx <- (c(-Inf, sX[1:(n - 1)]) < sX)
  newy <- (c(-Inf, sY[1:(n - 1)]) < sY)

  # Unique X/Y's
  x <- sX[newx]
  y <- sY[newy]

  # Number of unique X/Y values
  ell <- length(x)
  m <- length(y)

  # For each i = 1,...,n, compute k and j such that X_i = x_k and Y_i = y_j
  rsX <- cumsum(newx)
  iX <- rsX[rank(X, ties.method = "f")]
  rsY <- cumsum(newy)
  iY <- rsY[rank(Y, ties.method = "f")]

  # Compute w
  w <- matrix(0, nrow = ell, ncol = m)
  for (i in 1:n) {
    w[iX[i], iY[i]] <- w[iX[i], iY[i]] + W[i]
  }

  # Compute the weights w_j+ = sum_k wjk and w_+k = sum_j wjk
  w_j.plus <- rowSums(w)
  w_plus.k <- colSums(w)

  # Compute cumulative weights
  #     \underline{w}_jk := sum_{k'=k}^M_j w_jk'
  # for m_j <= k <= M_j and 1 <= j <= ell and
  #     \overline{w}_jk := sum_{j'=j}^L_k w_j'k
  # for ell_k <= j <= L_k and 1 <= k <= m and
  w_cumul.1 <- t(apply(w, 1, function(v) rev(cumsum(rev(v)))))
  w_cumul.2 <- apply(w, 2, function(v) rev(cumsum(rev(v))))

  # Determine mM
  mM <- matrix(0, ell, 2)
  tmp <- 1
  for (j in 1:ell) {
    tmp <- max(tmp, which(w[j, ] > 0))
    mM[j, 2] <- tmp
  }
  tmp <- m
  for (j in ell:1) {
    tmp <- min(tmp, which(w[j, ] > 0))
    mM[j, 1] <- tmp
  }

  # Determine lL
  lL <- matrix(0, m, 2)
  tmp <- 1
  for (k in 1:m) {
    tmp <- max(tmp, which(w[, k] > 0))
    lL[k, 2] <- tmp
  }
  tmp <- ell
  for (k in m:1) {
    tmp <- min(tmp, which(w[, k] > 0))
    lL[k, 1] <- tmp
  }

  # Determine PP
  PP <- matrix(FALSE, ell, m)
  for (j in 1:ell) {
    PP[j, mM[j, 1]:mM[j, 2]] <- TRUE
  }

  return(list(
    x = x, X = X, ell = ell, lL = lL,
    y = y, Y = Y, m = m, mM = mM,
    w = w, w_j.plus = w_j.plus, w_plus.k = w_plus.k,
    w_cumul.1 = w_cumul.1, w_cumul.2 = w_cumul.2,
    n = n, PP = PP
  ))
}