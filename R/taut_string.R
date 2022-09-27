#' Taut string
#'
#' @description The three vectors Y_low, Y_upp and x of length n > 2 are such
#' that
#'     Y_low_i <= Y_upp_i  for all i
#' with equality for i in {1,n}, and
#'     x1 < x2 < ... < xn.
#' The function outputs the corresponding taut string Y. That means, Y is a
#' vector of length n such that
#'     Y_low_i <= Y_i <= Y_upp_i  for all i,
#' for 1 < i < n,
#'     Y_i = Y_low_i  if  dY_i-1/dx_i-1 > dY_i/dx_i ,
#'     Y_i = Y_upp_i  if  dY_i-1/dx_i-1 < dY_i/dx_i ,
#' where dz_j := z_j+1 - z_j for any vector z and 1 <= j < length(z).
#'
#' @param Y_low Lower vector
#' @param Y_upp Upper vector
#' @param x Vector of support points
#'
#' @return Taut string between Y_low and Y_upp
#' @export
taut_string <- function(Y_low, Y_upp, x = 1:length(Y_low)) {
  n <- length(Y_low)
  a <- 1
  # Largest index up to which lower and upper string coincide.
  k_low <- c(1, 2, rep(NA, n - 2))
  # indices of knots of lower string
  slope_low <- c(Inf, (Y_low[2] - Y_low[1]) / (x[2] - x[1]), rep(NA, n - 2))
  #  left slopes of the lower string at knots
  b_low <- 2
  # index of last knot of lower string
  k_upp <- c(1, 2, rep(NA, n - 2))
  # indices of knots of upper string
  slope_upp <- c(-Inf, (Y_upp[2] - Y_upp[1]) / (x[2] - x[1]), rep(NA, n - 2))
  # left slopes of the upper string at knots
  b_upp <- 2
  # index of last know of upper string
  for (c in 3:n) {
    # Extend lower string:
    b_low <- b_low + 1
    k_low[b_low] <- c
    slope_low[b_low] <- (Y_low[c] - Y_low[c - 1]) / (x[c] - x[c - 1])
    # Pull lower string:
    while (b_low > a + 1 && slope_low[b_low] >= slope_low[b_low - 1]) {
      slope_low[b_low - 1] <-
        (Y_low[c] - Y_low[k_low[b_low - 2]]) / (x[c] - x[k_low[b_low - 2]])
      k_low[b_low - 1] <- c
      b_low <- b_low - 1
    }

    # Extend upper string:
    b_upp <- b_upp + 1
    k_upp[b_upp] <- c
    slope_upp[b_upp] <- (Y_upp[c] - Y_upp[c - 1]) / (x[c] - x[c - 1])
    # Pull upper string:
    while (b_upp > a + 1 && slope_upp[b_upp] <= slope_upp[b_upp - 1]) {
      slope_upp[b_upp - 1] <-
        (Y_upp[c] - Y_upp[k_upp[b_upp - 2]]) / (x[c] - x[k_upp[b_upp - 2]])
      k_upp[b_upp - 1] <- c
      b_upp <- b_upp - 1
    }

    # Check whether lower and upper strings are still ordered and modifiy them:
    if (b_low == a + 1) {
      # Bend lower string:
      while (b_upp > a + 1 && slope_low[b_low] >= slope_upp[a + 1]) {
        a <- a + 1
        b_low <- b_low + 1
        k_low[a] <- k_upp[a]
        slope_low[a] <- slope_upp[a]
        k_low[b_low] <- c
        slope_low[b_low] <-
          (Y_low[c] - Y_upp[k_upp[a]]) / (x[c] - x[k_upp[a]])
        Y_low[k_low[a]] <- Y_upp[k_upp[a]]
      }
    }
    if (b_upp == a + 1) {
      # Bend upper string:
      while (b_low > a + 1 && slope_upp[b_upp] <= slope_low[a + 1]) {
        a <- a + 1
        b_upp <- b_upp + 1
        k_upp[a] <- k_low[a]
        slope_upp[a] <- slope_low[a]
        k_upp[b_upp] <- c
        slope_upp[b_upp] <-
          (Y_upp[c] - Y_low[k_low[a]]) / (x[c] - x[k_low[a]])
        Y_upp[k_upp[a]] <- Y_low[k_low[a]]
      }
    }
  }

  # Transform slope_low (== slope_upp) into a proper vector of length n:
  slope <- c(NA, rep(0, n - 1))
  a <- 1
  for (b in 2:b_low) {
    slope[(a + 1):k_low[b]] <- slope_low[b]
    a <- k_low[b]
  }

  # Compute vector Y from its slopes:
  Y <- Y_low
  dx <- x[2:n] - x[1:(n - 1)]
  Y[2:n] <- Y[1] + cumsum(dx * slope[2:n])

  return(list(knot_ind = k_low, Y = Y))
}
