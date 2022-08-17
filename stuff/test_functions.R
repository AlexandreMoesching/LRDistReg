####____________________________________________________________________________
#### SETUP                                                                  ####

n <- 1e5; ell0 <- 1e2; m0 <- 1e2
# set.seed(1)
X <- sample(rnorm(ell0), n, replace = TRUE)
Y <- sample(rnorm(m0), n, replace = TRUE)
W <- rep(1, n)

####____________________________________________________________________________
#### TEST prepare_data                                                      ####

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

# plot(X, Y, pch = 16, cex = 0.2)
# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,1]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,2]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[res_R$lL[,1]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)
# points(res_R$x[res_R$lL[,2]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)

all(res_R$PP == res_cpp$PP)

# microbenchmark::microbenchmark(prepare.data(X, Y, W),
#                                prepare_data_cpp(X, Y, W),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST reparametrize                                                     ####

ell <- res_cpp$ell
m <- res_cpp$m

lambda <- matrix(rnorm(ell * m), nrow = ell, ncol = m)

lL <- res_cpp$lL
mM <- res_cpp$mM

all.equal(lambda1.to.theta(lambda, ell, m, mM + 1),
          lambda1_to_theta_cpp(lambda, ell, m, mM), tolerance = 1e-14)
all.equal(lambda2.to.theta(lambda, ell, m, lL + 1),
          lambda2_to_theta_cpp(lambda, ell, m, lL), tolerance = 1e-14)

# microbenchmark::microbenchmark(lambda1.to.theta(lambda, ell, m, mM + 1),
#                                lambda1_to_theta_cpp(lambda, ell, m, mM),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST likelihood_functions                                              ####

theta <- matrix(rnorm(ell * m), nrow = ell, ncol = m)

w_ul <- res_cpp$w_ul

tmp_R <- vgamma.tilde1(theta, ell, m, n, mM + 1, w_ul)
tmp_cpp <- vgamma_tilde1_cpp(theta, ell, m, n, mM, w_ul)
all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14)
all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14)

w_ol <- res_cpp$w_ol

tmp_R <- vgamma.tilde2(theta, ell, m, n, lL + 1, w_ol)
tmp_cpp <- vgamma_tilde2_cpp(theta, ell, m, n, lL, w_ol)
all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14)
all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14)

w <- res_cpp$w
PP <- res_R$PP

f.theta(theta, n, w, PP) == f_theta_cpp(theta, ell, n, mM, w)

# microbenchmark::microbenchmark(vgamma.tilde1(theta, ell, m, n, mM + 1, w_ul),
#                                vgamma_tilde1_cpp(theta, ell, m, n, mM, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(vgamma.tilde2(theta, ell, m, n, lL + 1, w_ol),
#                                vgamma_tilde2_cpp(theta, ell, m, n, lL, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(f.theta(theta, n, w, PP),
#                                f_theta_cpp(theta, ell, n, mM, w),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST calibrate                                                         ####

w_jplus <- c(res_cpp$w_jplus)
w_plusk <- c(res_cpp$w_plusk)

all.equal(calibrate1(theta, n, w_jplus),
          calibrate1_cpp(theta, n, w_jplus),
          tolerance = 1e-14)
all.equal(calibrate2(theta, n, w_plusk),
          calibrate2_cpp(theta, n, w_plusk),
          tolerance = 1e-14)

prec <- 1e-10

all.equal(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
          calibrate_cpp(theta, ell, mM, n, w, w_jplus, w_plusk, prec),
          tolerance = 1e-14)

# microbenchmark::microbenchmark(calibrate1(theta, n, w_jplus),
#                                calibrate1_cpp(theta, n, w_jplus),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate2(theta, n, w_plusk),
#                                calibrate2_cpp(theta, n, w_plusk),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
#                                calibrate_cpp(theta, ell, mM, n, w, w_jplus, w_plusk, prec),
#                                times = 1e2)
