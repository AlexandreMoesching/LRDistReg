rm(list = ls())
library(LRDistReg)
library(doSNOW)

set.seed(2)

# Parameters
l0 <- 1e1
n <- 1e1
delta0 <- 1e-5
rel.tol <- 1e-5
n.sim <- 1e3

# Shape and scale functions
a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)

# True family of distributions
r.cond.dist <- function(x)        rgamma(1,     shape = a(x), scale = b(x))
d.cond.dist <- function(x, alpha) dgamma(alpha, shape = a(x), scale = b(x))
p.cond.dist <- function(x, alpha) pgamma(alpha, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))

# Potential covariates
x0 <- 1 + 3 * (1:l0) / l0

result <- matrix(0, nrow = n.sim, ncol = 6*l0)
for (it in 1:n.sim) {
  cat("it = ", it, "\n", sep = "")
  # Generate a sample
  X <- sort(sample(x0, size = n, replace = TRUE))
  Y <- rep(0, n)
  for (i in 1:n) Y[i] <- r.cond.dist(X[i])
  Y <- round(Y, 2)
  W <- rep(1, n)

  # Estimate CDF's
  cat("Start estimating\n", sep = "")
  res <- dist_reg_C(X, Y, W, delta0, x0, TRUE)
  cat("Done\n", sep = "")
  cat("delta = ", res$delta, "\n", sep = "")
  # cat("CDF_LR = ", res$CDF_LR, "\n", sep = "")

  # Compute SS and CRPS
  cat("Start computing scores\n", sep = "")
  res_SS   <- SS(  x0, l0, res, a, b, rel.tol)
  res_CRPS <- CRPS(x0, l0, res, a, b, rel.tol)
  cat("Done\n", sep = "")

  # Return
  result[i,] <- c(res_SS[  ,1], res_SS[  ,2], res_SS[  ,3],
                  res_CRPS[,1], res_CRPS[,2], res_CRPS[,3])
  cat("\n")
}

################################################################################
# DEBUG

# it <- 44
# cat("it = ", it, "\n", sep = "")
# # Generate a sample
# X <- sort(sample(x0, size = n, replace = TRUE))
# Y <- rep(0, n)
# for (i in 1:n) Y[i] <- r.cond.dist(X[i])
# Y <- round(Y, 2)
# W <- rep(1, n)
#
# plot(X, Y)
#
# # Estimate CDF's
# cat("Start estimating\n", sep = "")
# res <- dist_reg_cpp(X, Y, W, delta0, x0, TRUE)
# cat("Done\n", sep = "")
# cat("delta = ", res$delta, "\n", sep = "")
# # cat("CDF_LR = ", res$CDF_LR, "\n", sep = "")
#
# res$CDF_LR
#
# res <- dist.reg(X, Y, W, FALSE, FALSE, FALSE, delta0, x0, TRUE)
#
# # Compute SS and CRPS
# cat("Start computing scores\n", sep = "")
# res_SS   <- SS(  x0, l0, res, a, b, rel.tol)
# res_CRPS <- CRPS(x0, l0, res, a, b, rel.tol)
# cat("Done\n", sep = "")
#
# # Return
# result[i,] <- c(res_SS[  ,1], res_SS[  ,2], res_SS[  ,3],
#                 res_CRPS[,1], res_CRPS[,2], res_CRPS[,3])
# cat("\n")
