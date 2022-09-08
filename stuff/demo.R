# Clean environment and load package
rm(list = ls())
library(LRDistReg)

#_______________________________________________________________________________
# TRUE FAMILY OF DISTRIBUTIONS
#
# Define true family of distributions
a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)
r.cond.dist <- function(x)        rgamma(1,     shape = a(x), scale = b(x))
d.cond.dist <- function(x, alpha) dgamma(alpha, shape = a(x), scale = b(x))
p.cond.dist <- function(x, alpha) pgamma(alpha, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))

# Plot of true family
xx <- seq(1, 4, length.out = 2e2)
yy <- seq(q.cond.dist(1, 0.05), q.cond.dist(4, 0.95), length.out = 2e2)
lattice::levelplot(outer(xx, yy, FUN = "d.cond.dist"),
                   nlevels = 1e2,
                   aspect = "fill",
                   col.regions = hcl.colors(1e2),
                   row.values = xx,
                   column.values = yy,
                   xlab = expression(italic(x)),
                   ylab = expression(italic(y)),
                   xlim = range(xx),
                   ylim = range(yy))

#_______________________________________________________________________________
# GENERATE A SAMPLE
#
ell0 <- 5e0                   # Number of potential covariates
x0 <- 1 + 3 * (1:ell0) / ell0 # Potential covariates

n <- 1e3                     # Sample size
X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
Y <- round(Y, 1)              # Creates some ties
W <- rep(1, n)                # Sample weights

delta0 <- 1e-8                # Threshold for estimation precision

#_______________________________________________________________________________
# ESTIMATE
#
# (1) Estimation with R
par_R <- prepare.data(X, Y, W)
res_R <- TP2.fit(par_R, delta0, echo = FALSE, out.file = FALSE)

# (2) Estimation with C++
res_cpp <- TP2_fit_cpp(X, Y, W, delta0)

# Compare results
sum(abs(res_R$h.TP2 - res_cpp$h_TP2))
sum(abs(res_R$q.LR - res_cpp$q_LR))
sum(abs(res_R$CDF.LR - res_cpp$CDF_LR))
