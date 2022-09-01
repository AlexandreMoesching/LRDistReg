library(LRDistReg)
rm(list = ls())

a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)
r.cond.dist <- function(x) rgamma(1, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))

xx <- seq(1, 4, length.out = 2e2)
yy <- seq(q.cond.dist(1, 0.05), q.cond.dist(4, 0.95), length.out = 2e2)
lattice::levelplot(outer(xx, yy, FUN = "d.cond.dist"), nlevels = 1e2,
                   aspect = "fill", col.regions = hcl.colors(1e2),
                   row.values = xx, column.values = yy,
                   xlab = expression(italic(x)), ylab = expression(italic(y)),
                   xlim = range(xx), ylim = range(yy))

n <- 1e4
ell0 <- 1e2
x0 <- 1+(1:ell0)/ell0*3

X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
Y <- round(Y, 1) # Should create some ties
W <- rep(1, n)

delta0 <- 1e-6

par_R <- prepare.data(X, Y, W)
res_R <- TP2.fit(par_R, delta0, echo = FALSE, out.file = FALSE)

res_cpp <- TP2_fit_cpp(X, Y, W, delta0)

max(abs(res_R$h.TP2 - res_cpp$h_TP2))

# microbenchmark::microbenchmark(TP2.fit(prepare.data(X, Y, W), delta0, echo = FALSE, out.file = FALSE),
#                                TP2_fit_cpp(X, Y, W, delta0),
#                                times = 10)
