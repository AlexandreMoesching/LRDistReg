library(LRDistReg)
rm(list = ls())

a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)
r.cond.dist <- function(x) rgamma(1, shape = a(x), scale = b(x))
d.cond.dist <- function(x, alpha) dgamma(alpha, shape = a(x), scale = b(x))
p.cond.dist <- function(x, alpha) pgamma(alpha, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))

xx <- seq(1, 4, length.out = 2e2)
yy <- seq(q.cond.dist(1, 0.05), q.cond.dist(4, 0.95), length.out = 2e2)
lattice::levelplot(outer(xx, yy, FUN = "d.cond.dist"), nlevels = 1e2,
                   aspect = "fill", col.regions = hcl.colors(1e2),
                   row.values = xx, column.values = yy,
                   xlab = expression(italic(x)), ylab = expression(italic(y)),
                   xlim = range(xx), ylim = range(yy))

n <- 2e2
ell0 <- 1e2
x0 <- 1+(1:ell0)/ell0*3

X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
Y <- round(Y, 1) # Should create some ties
W <- rep(1, n)

delta0 <- 1e-8

par_R <- prepare.data(X, Y, W)
res_R <- TP2.fit(par_R, delta0, echo = FALSE, out.file = FALSE)

res_cpp <- TP2_fit_cpp(X, Y, W, delta0)

sum(abs(res_R$h.TP2 - res_cpp$h_TP2))
sum(abs(res_R$q.LR - res_cpp$q_LR))
sum(abs(res_R$CDF.LR - res_cpp$CDF_LR))

n_boot <- 1e2

h_boot_manual <- matrix(0, nrow = nrow(res_R$h.TP2), ncol = ncol(res_R$h.TP2))
CDF_boot_manual <- matrix(0, nrow = nrow(res_R$h.TP2), ncol = ncol(res_R$h.TP2))

for (i in 1:n_boot) {
  W <- rexp(n, rate = 0.5)
  W <- W * (n/sum(W))
  tmp <- TP2_fit_cpp(X, Y, W, delta0)
  h_boot_manual <- h_boot_manual + tmp$h_TP2 / n_boot
  CDF_boot_manual <- CDF_boot_manual + tmp$CDF_LR / n_boot
}

res_bag <- TP2_fit_bag_cpp(X, Y, W, delta0, n_boot)
max(abs(res_bag$h_TP2 - res_cpp$h_TP2))

beta <- (1:9)/10
plot(0, xlim = range(yy), ylim = c(0,1), xlab = expression(italic(y)), ylab = "CDF")
for (s in seq_along(beta)) {
  j0 <- floor(par_R$ell * beta[s])
  x0 <- par_R$x[j0]
  lines(yy, p.cond.dist(x0, yy))
  lines(par_R$y, res_cpp$CDF_LR[j0,], col = "blue")
  lines(par_R$y, res_bag$CDF_LR[j0,], col = "red")
  lines(par_R$y, CDF_boot_manual[j0,], col = "green")
}



# microbenchmark::microbenchmark(TP2.fit(prepare.data(X, Y, W), delta0, echo = FALSE, out.file = FALSE),
#                                TP2_fit_cpp(X, Y, W, delta0),
#                                times = 10)


