n <- 1e5; ell0 <- 1e2; m0 <- 1e2
# set.seed(1)
X <- sample(rnorm(ell0), n, replace = TRUE)
Y <- sample(rnorm(m0), n, replace = TRUE)
# plot(X, Y, pch = 16, cex = 0.2)
W <- rep(1, n)

res_R <- prepare.data(X, Y, W)
res_cpp <- prepare_data_cpp(X, Y, W)

all(res_R$x == res_cpp$x)
all(res_R$X == res_cpp$X)
all(res_R$y == res_cpp$y)
all(res_R$Y == res_cpp$Y)
all(res_R$ell == res_cpp$ell)
all(res_R$m == res_cpp$m)
all(res_R$w == res_cpp$w)
all(res_R$w_jplus == res_cpp$w_jplus)
all(res_R$w_plusk == res_cpp$w_plusk)
all(res_R$w_ul == res_cpp$w_ul)
all(res_R$w_ol == res_cpp$w_ol)

all(res_R$lL == res_cpp$lL + 1)
all(res_R$mM == res_cpp$mM + 1)

# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,1]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,2]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[res_R$lL[,1]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)
# points(res_R$x[res_R$lL[,2]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)

all(res_R$PP == res_cpp$PP)

# microbenchmark::microbenchmark(prepare.data(X, Y, W),
#                                prepare_data_cpp(X, Y, W),
#                                times = 1e2)

ell <- res_cpp$ell
m <- res_cpp$m
lL <- res_cpp$lL
mM <- res_cpp$mM

lambda <- matrix(rnorm(ell * m), nrow = ell, ncol = m)

all.equal(lambda1.to.theta(lambda, ell, m, mM + 1),
          lambda1_to_theta_cpp(lambda, ell, m, mM), tolerance = 1e-14)
all.equal(lambda2.to.theta(lambda, ell, m, lL + 1),
          lambda2_to_theta_cpp(lambda, ell, m, lL), tolerance = 1e-14)

# microbenchmark::microbenchmark(lambda1.to.theta(lambda, ell, m, mM + 1),
#                                lambda1_to_theta_cpp(lambda, ell, m, mM),
#                                times = 1e2)
