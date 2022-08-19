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

cat(all(res_R$x == res_cpp$x), "\n")
cat(all(res_R$X == res_cpp$X), "\n")
cat(all(res_R$y == res_cpp$y), "\n")
cat(all(res_R$Y == res_cpp$Y), "\n")
cat(all(res_R$ell == res_cpp$ell), "\n")
cat(all(res_R$m == res_cpp$m), "\n")
cat(all(res_R$w == res_cpp$w), "\n")
cat(all(res_R$w_jplus == res_cpp$w_jplus), "\n")
cat(all(res_R$w_plusk == res_cpp$w_plusk), "\n")
cat(all(res_R$w_ul == res_cpp$w_ul), "\n")
cat(all(res_R$w_ol == res_cpp$w_ol), "\n")

cat(all(res_R$lL == res_cpp$lL + 1), "\n")
cat(all(res_R$mM == res_cpp$mM + 1), "\n")

# plot(X, Y, pch = 16, cex = 0.2)
# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,1]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[1:res_R$ell], res_R$y[res_R$mM[,2]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[res_R$lL[,1]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)
# points(res_R$x[res_R$lL[,2]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)

cat(all(res_R$PP == res_cpp$PP), "\n")

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

cat(all.equal(lambda1.to.theta(lambda, ell, m, mM + 1),
              lambda1_to_theta_cpp(lambda, ell, m, mM),
              tolerance = 1e-14), "\n")
cat(all.equal(lambda2.to.theta(lambda, ell, m, lL + 1),
              lambda2_to_theta_cpp(lambda, ell, m, lL),
              tolerance = 1e-14), "\n")

# microbenchmark::microbenchmark(lambda1.to.theta(lambda, ell, m, mM + 1),
#                                lambda1_to_theta_cpp(lambda, ell, m, mM),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST likelihood_functions                                              ####

theta <- matrix(rnorm(ell * m), nrow = ell, ncol = m)

w_ul <- res_cpp$w_ul

tmp_R <- vgamma.tilde1(theta, ell, m, n, mM + 1, w_ul)
tmp_cpp <- vgamma_tilde1_cpp(theta, ell, m, n, mM, w_ul)
cat(all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14), "\n")
cat(all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14), "\n")

w_ol <- res_cpp$w_ol

tmp_R <- vgamma.tilde2(theta, ell, m, n, lL + 1, w_ol)
tmp_cpp <- vgamma_tilde2_cpp(theta, ell, m, n, lL, w_ol)
cat(all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14), "\n")
cat(all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14), "\n")

w <- res_cpp$w
PP <- res_R$PP

cat(all.equal(f.theta(theta, n, w, PP),
              f_theta_cpp(theta, ell, n, mM, w),
              tolerance = 1e-12), "\n")

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

cat(all.equal(calibrate1(theta, n, w_jplus),
              calibrate1_cpp(theta, n, w_jplus),
              tolerance = 1e-14), "\n")
cat(all.equal(calibrate2(theta, n, w_plusk),
              calibrate2_cpp(theta, n, w_plusk),
              tolerance = 1e-14), "\n")

prec <- 1e-12

cat(all.equal(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
              calibrate_cpp(theta, ell, mM, n, w, w_jplus, w_plusk, prec),
              tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(calibrate1(theta, n, w_jplus),
#                                calibrate1_cpp(theta, n, w_jplus),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate2(theta, n, w_plusk),
#                                calibrate2_cpp(theta, n, w_plusk),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
#                                calibrate_cpp(theta, ell, mM, n, w, w_jplus, w_plusk, prec),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST local_search_functions                                            ####

res_R <- local.search1(theta, ell, m, n, mM + 1, lL + 1, PP, w, w_ul)
res_cpp <- local_search1_cpp(theta, ell, m, n, lL, mM, w, w_ul)
cat(all.equal(res_R$Psi, res_cpp$Psi, tolerance = 1e-12), "\n")
cat(all.equal(res_R$delta, res_cpp$delta, tolerance = 1e-12), "\n")

res_R <- local.search2(theta, ell, m, n, mM + 1, lL + 1, PP, w, w_ol)
res_cpp <- local_search2_cpp(theta, ell, m, n, lL, mM, w, w_ol)
cat(all.equal(res_R$Psi, res_cpp$Psi, tolerance = 1e-12), "\n")
cat(all.equal(res_R$delta, res_cpp$delta, tolerance = 1e-12), "\n")

Psi <- res_cpp$Psi
delta <- res_cpp$delta

res_R <- simple.step(theta, Psi, delta, ell, m, n, w, PP)
res_cpp <- simple_step_cpp(theta, Psi, delta, ell, mM, n, w)
cat(all.equal(res_R[PP], res_cpp[PP], tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(local.search1(theta, ell, m, n, mM + 1, lL + 1, PP, w, w_ul),
#                                local_search1_cpp(theta, ell, m, n, lL, mM, w, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(local.search2(theta, ell, m, n, mM + 1, lL + 1, PP, w, w_ol),
#                                local_search2_cpp(theta, ell, m, n, lL, mM, w, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(simple.step(theta, Psi, delta, ell, m, n, w, PP),
#                                simple_step_cpp(theta, Psi, delta, ell, mM, n, w),
#                                times = 1e2)


####____________________________________________________________________________
#### TEST main                                                              ####
delta0 <- 1e-8

par_R <- prepare.data(X, Y, W)
res_R <- TP2.fit(par_R, delta0, echo = TRUE, out.file = FALSE)

res_cpp <- TP2_fit_cpp(X, Y, W, delta0)

all.equal(res_R$h.TP2, res_cpp$h_TP2, tolerance = 1e-10)

