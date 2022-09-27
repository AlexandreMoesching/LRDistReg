#' Taut string
#'
#' @description The three vectors Y.low, Y.upp and x of length n > 2 are such
#' that
#'     Y.low.i <= Y.upp.i  for all i
#' with equality for i in {1,n}, and
#'     x1 < x2 < ... < xn.
#' The function outputs the corresponding taut string Y. That means, Y is a
#' vector of length n such that
#'     Y.low.i <= Y.i <= Y.upp.i  for all i,
#' for 1 < i < n,
#'     Y.i = Y.low.i  if  dY.i-1/dx.i-1 > dY.i/dx.i ,
#'     Y.i = Y.upp.i  if  dY.i-1/dx.i-1 < dY.i/dx.i ,
#' where dz.j := z.j+1 - z.j for any vector z and 1 <= j < length(z).
#'
#' @param Y.low Lower vector
#' @param Y.upp Upper vector
#' @param x Vector of support points
#'
#' @return Taut string between Y.low and Y.upp
#' @export
taut_string <- function(Y.low, Y.upp, x = 1:length(Y.low)) {
  n <- length(Y.low)
  a <- 1
  # Largest index up to which lower and upper string coincide.
  k.low <- c(1, 2, rep(NA, n - 2))
  # indices of knots of lower string
  slope.low <- c(Inf, (Y.low[2] - Y.low[1]) / (x[2] - x[1]), rep(NA, n - 2))
  #  left slopes of the lower string at knots
  b.low <- 2
  # index of last knot of lower string
  k.upp <- c(1, 2, rep(NA, n - 2))
  # indices of knots of upper string
  slope.upp <- c(-Inf, (Y.upp[2] - Y.upp[1]) / (x[2] - x[1]), rep(NA, n - 2))
  # left slopes of the upper string at knots
  b.upp <- 2
  # index of last know of upper string
  for (c in 3:n) {
    # Extend lower string:
    b.low <- b.low + 1
    k.low[b.low] <- c
    slope.low[b.low] <- (Y.low[c] - Y.low[c - 1]) / (x[c] - x[c - 1])
    # Pull lower string:
    while (b.low > a + 1 && slope.low[b.low] >= slope.low[b.low - 1]) {
      slope.low[b.low - 1] <-
        (Y.low[c] - Y.low[k.low[b.low - 2]]) / (x[c] - x[k.low[b.low - 2]])
      k.low[b.low - 1] <- c
      b.low <- b.low - 1
    }

    # Extend upper string:
    b.upp <- b.upp + 1
    k.upp[b.upp] <- c
    slope.upp[b.upp] <- (Y.upp[c] - Y.upp[c - 1]) / (x[c] - x[c - 1])
    # Pull upper string:
    while (b.upp > a + 1 && slope.upp[b.upp] <= slope.upp[b.upp - 1]) {
      slope.upp[b.upp - 1] <-
        (Y.upp[c] - Y.upp[k.upp[b.upp - 2]]) / (x[c] - x[k.upp[b.upp - 2]])
      k.upp[b.upp - 1] <- c
      b.upp <- b.upp - 1
    }

    # Check whether lower and upper strings are still ordered and modifiy them:
    if (b.low == a + 1) {
      # Bend lower string:
      while (b.upp > a + 1 && slope.low[b.low] >= slope.upp[a + 1]) {
        a <- a + 1
        b.low <- b.low + 1
        k.low[a] <- k.upp[a]
        slope.low[a] <- slope.upp[a]
        k.low[b.low] <- c
        slope.low[b.low] <-
          (Y.low[c] - Y.upp[k.upp[a]]) / (x[c] - x[k.upp[a]])
        Y.low[k.low[a]] <- Y.upp[k.upp[a]]
      }
    }
    if (b.upp == a + 1) {
      # Bend upper string:
      while (b.low > a + 1 && slope.upp[b.upp] <= slope.low[a + 1]) {
        a <- a + 1
        b.upp <- b.upp + 1
        k.upp[a] <- k.low[a]
        slope.upp[a] <- slope.low[a]
        k.upp[b.upp] <- c
        slope.upp[b.upp] <-
          (Y.upp[c] - Y.low[k.low[a]]) / (x[c] - x[k.low[a]])
        Y.upp[k.upp[a]] <- Y.low[k.low[a]]
      }
    }
  }

  # Transform slope.low (== slope.upp) into a proper vector of length n:
  slope <- c(NA, rep(0, n - 1))
  a <- 1
  for (b in 2:b.low) {
    slope[(a + 1):k.low[b]] <- slope.low[b]
    a <- k.low[b]
  }

  # Compute vector Y from its slopes:
  Y <- Y.low
  dx <- x[2:n] - x[1:(n - 1)]
  Y[2:n] <- Y[1] + cumsum(dx * slope[2:n])

  return(list(knot.ind = k.low, Y = Y))
}
