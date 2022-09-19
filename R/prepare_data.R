#' Prepare the data, R version
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified weights
#'
#' @return A list of pre-computed parameters necessary for the estimation
#'
#' @export
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
  l <- length(x)
  m <- length(y)

  # For each i = 1,...,n, compute k and j such that X_i = x_k and Y_i = y_j
  rsX <- cumsum(newx)
  iX <- rsX[rank(X, ties.method = "f")]
  rsY <- cumsum(newy)
  iY <- rsY[rank(Y, ties.method = "f")]

  # Compute w
  w <- matrix(0, nrow = l, ncol = m)
  for (i in 1:n) {
    w[iX[i], iY[i]] <- w[iX[i], iY[i]] + W[i]
  }

  # Compute the weights w_j+ = sum_k wjk and w_+k = sum_j wjk
  w_jplus <- rowSums(w)
  w_plusk <- colSums(w)

  # Compute cumulative weights
  #     \underline{w}_jk := sum_{k'=k}^M_j w_jk'
  # for m_j <= k <= M_j and 1 <= j <= l and
  #     \overline{w}_jk := sum_{j'=j}^L_k w_j'k
  # for l_k <= j <= L_k and 1 <= k <= m and
  w_ul <- t(apply(w, 1, function(v) rev(cumsum(rev(v)))))
  w_ol <- apply(w, 2, function(v) rev(cumsum(rev(v))))

  # Determine mM
  mM <- matrix(0, l, 2)
  tmp <- 1
  for (j in 1:l) {
    tmp <- max(tmp, which(w[j, ] > 0))
    mM[j, 2] <- tmp
  }
  tmp <- m
  for (j in l:1) {
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
  tmp <- l
  for (k in m:1) {
    tmp <- min(tmp, which(w[, k] > 0))
    lL[k, 1] <- tmp
  }

  # Determine PP
  PP <- matrix(FALSE, l, m)
  for (j in 1:l) {
    PP[j, mM[j, 1]:mM[j, 2]] <- TRUE
  }

  return(list(
    l = l,
    lL = lL,
    m = m,
    mM = mM,
    n = n,
    PP = PP,
    w = w,
    w_jplus = w_jplus,
    w_plusk = w_plusk,
    w_ul = w_ul,
    w_ol = w_ol,
    W = W,
    x = x,
    X = X,
    y = y,
    Y = Y
  ))
}
