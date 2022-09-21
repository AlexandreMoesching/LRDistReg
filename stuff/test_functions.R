library(LRDistReg)

####____________________________________________________________________________
#### SETUP                                                                  ####

n <- 1e4; l0 <- 1e2; m0 <- 1e2
# set.seed(1)
X <- sample(rnorm(l0), n, replace = TRUE)
Y <- sample(rnorm(m0), n, replace = TRUE)
W <- rep(1, n)

####____________________________________________________________________________
#### TEST prepare_data                                                      ####

res_R <- prepare.data(X, Y, W)
res_cpp <- prepare_data_cpp(X, Y, W)

cat("x's are the same:", all(res_R$x == res_cpp$x), "\n")
cat("X's are the same:", all(res_R$X == res_cpp$X), "\n")
cat("y's are the same:", all(res_R$y == res_cpp$y), "\n")
cat("Y's are the same:", all(res_R$Y == res_cpp$Y), "\n")
cat("l's are the same:", all(res_R$l == res_cpp$l), "\n")
cat("m's are the same:", all(res_R$m == res_cpp$m), "\n")
cat("w's are the same:", all(res_R$w == res_cpp$w), "\n")
cat("w_jplus's are the same:", all(res_R$w_jplus == res_cpp$w_jplus), "\n")
cat("w_plusk's are the same:", all(res_R$w_plusk == res_cpp$w_plusk), "\n")
cat("w_ul's are the same:", all(res_R$w_ul == res_cpp$w_ul), "\n")
cat("w_ol's are the same:", all(res_R$w_ol == res_cpp$w_ol), "\n")
cat("lL's are the same:", all(res_R$lL == res_cpp$lL + 1), "\n")
cat("mM's are the same:", all(res_R$mM == res_cpp$mM + 1), "\n")

# plot(X, Y, pch = 16, cex = 0.2)
# points(res_R$x[1:res_R$l], res_R$y[res_R$mM[,1]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[1:res_R$l], res_R$y[res_R$mM[,2]], col = "red", pch = 1, cex = 0.7)
# points(res_R$x[res_R$lL[,1]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)
# points(res_R$x[res_R$lL[,2]], res_R$y[1:res_R$m], col = "blue", pch = 16, cex = 0.5)

cat("PP's are the same:", all(res_R$PP == res_cpp$PP), "\n")

# microbenchmark::microbenchmark(prepare.data(X, Y, W),
#                                prepare_data_cpp(X, Y, W),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST reparametrize                                                     ####

l <- res_cpp$l
m <- res_cpp$m

lambda <- matrix(rnorm(l * m), nrow = l, ncol = m)

lL <- res_cpp$lL
mM <- res_cpp$mM

cat("lambda1.to.theta's are the same:",
    all.equal(lambda1.to.theta(lambda, l, m, mM + 1),
              lambda1_to_theta_cpp(lambda, l, m, mM),
              tolerance = 1e-14), "\n")
cat("lambda2.to.theta's are the same:",
    all.equal(lambda2.to.theta(lambda, l, m, lL + 1),
              lambda2_to_theta_cpp(lambda, l, m, lL),
              tolerance = 1e-14), "\n")

# microbenchmark::microbenchmark(lambda1.to.theta(lambda, l, m, mM + 1),
#                                lambda1_to_theta_cpp(lambda, l, m, mM),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST likelihood_functions                                              ####

theta <- matrix(rnorm(l * m), nrow = l, ncol = m)

w_ul <- res_cpp$w_ul

tmp_R <- vgamma.tilde1(theta, l, m, n, mM + 1, w_ul)
tmp_cpp <- vgamma_tilde1_cpp(theta, l, m, n, mM, w_ul)
cat("v1's are the same:",
    all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14), "\n")
cat("gamma1's are the same:",
    all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14), "\n")

w_ol <- res_cpp$w_ol

tmp_R <- vgamma.tilde2(theta, l, m, n, lL + 1, w_ol)
tmp_cpp <- vgamma_tilde2_cpp(theta, l, m, n, lL, w_ol)
cat("v2's are the same:",
    all.equal(tmp_R$v, tmp_cpp$v, tolerance = 1e-14), "\n")
cat("gamma2's are the same:",
    all.equal(tmp_R$gamma, tmp_cpp$gamma, tolerance = 1e-14), "\n")

w <- res_cpp$w
PP <- res_R$PP

cat("f.theta's are the same:",
    all.equal(f.theta(theta, n, w, PP),
              f_theta_cpp(theta, l, n, mM, w),
              tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(vgamma.tilde1(theta, l, m, n, mM + 1, w_ul),
#                                vgamma_tilde1_cpp(theta, l, m, n, mM, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(vgamma.tilde2(theta, l, m, n, lL + 1, w_ol),
#                                vgamma_tilde2_cpp(theta, l, m, n, lL, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(f.theta(theta, n, w, PP),
#                                f_theta_cpp(theta, l, n, mM, w),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST calibrate                                                         ####

w_jplus <- c(res_cpp$w_jplus)
w_plusk <- c(res_cpp$w_plusk)

cat("calibrate1's are the same:",
    all.equal(calibrate1(theta, n, w_jplus),
              calibrate1_cpp(theta, n, w_jplus),
              tolerance = 1e-14), "\n")
cat("calibrate2's are the same:",
    all.equal(calibrate2(theta, n, w_plusk),
              calibrate2_cpp(theta, n, w_plusk),
              tolerance = 1e-14), "\n")

prec <- 1e-12

cat("calibrate's are the same:",
    all.equal(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
              calibrate_cpp(theta, l, mM, n, w, w_jplus, w_plusk, prec),
              tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(calibrate1(theta, n, w_jplus),
#                                calibrate1_cpp(theta, n, w_jplus),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate2(theta, n, w_plusk),
#                                calibrate2_cpp(theta, n, w_plusk),
#                                times = 1e2)
# microbenchmark::microbenchmark(calibrate(theta, n, w, w_jplus, w_plusk, PP, prec),
#                                calibrate_cpp(theta, l, mM, n, w, w_jplus, w_plusk, prec),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST local_search_functions                                            ####

res_R <- local.search1(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ul)
res_cpp <- local_search1_cpp(theta, l, m, n, lL, mM, w, w_ul)
cat("Psi1's are the same:",
    all.equal(res_R$Psi, res_cpp$Psi, tolerance = 1e-12), "\n")
cat("delta1's are the same:",
    all.equal(res_R$delta, res_cpp$delta, tolerance = 1e-12), "\n")

res_R <- local.search2(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ol)
res_cpp <- local_search2_cpp(theta, l, m, n, lL, mM, w, w_ol)
cat("Psi2's are the same:",
    all.equal(res_R$Psi, res_cpp$Psi, tolerance = 1e-12), "\n")
cat("delta2's are the same:",
    all.equal(res_R$delta, res_cpp$delta, tolerance = 1e-12), "\n")

Psi <- res_cpp$Psi
delta <- res_cpp$delta

res_R <- simple.step(theta, Psi, delta, l, m, n, w, PP)
res_cpp <- simple_step_cpp(theta, Psi, delta, l, mM, n, w)
cat("simple.step's are the same:",
    all.equal(res_R[PP], res_cpp[PP], tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(local.search1(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ul),
#                                local_search1_cpp(theta, l, m, n, lL, mM, w, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(local.search2(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ol),
#                                local_search2_cpp(theta, l, m, n, lL, mM, w, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(simple.step(theta, Psi, delta, l, m, n, w, PP),
#                                simple_step_cpp(theta, Psi, delta, l, mM, n, w),
#                                times = 1e2)


####____________________________________________________________________________
#### TEST main                                                              ####
delta0 <- 1e-9

res_R <- dist.reg(X, Y, W,
                  show.design = FALSE, indices = FALSE, suggest.delta0 = FALSE,
                  delta0 = delta0, x0 = NULL, ST = TRUE)
res_cpp <- dist_reg_cpp(X, Y, W, delta0, x0 = numeric(), ST = TRUE)

cat("h_TP2's are the same:",
    all.equal(res_R$h_TP2, res_cpp$h_TP2, tolerance = 1e-10), "\n")
cat("CDF_LR's are the same:",
    all.equal(res_R$CDF_LR, res_cpp$CDF_LR, tolerance = 1e-10), "\n")
cat("CDF_ST's are the same:",
    all.equal(res_R$CDF_ST, res_cpp$CDF_ST, tolerance = 1e-10), "\n")
cat("CDF_EMP's are the same:",
    all.equal(res_R$CDF_EMP, res_cpp$CDF_EMP, tolerance = 1e-10), "\n")

####____________________________________________________________________________
#### TEST main with interpolate                                             ####
delta0 <- 1e-9
x0 <- unique(sort(floor(X)))

res_R <- dist.reg(X, Y, W,
                  show.design = FALSE, indices = FALSE, suggest.delta0 = FALSE,
                  delta0 = delta0, x0 = x0, ST = TRUE)
res_cpp <- dist_reg_cpp(X, Y, W, delta0, x0 = x0, ST = TRUE)

cat("CDF0_LR's are the same:",
    all.equal(res_R$CDF_LR, res_cpp$CDF_LR, tolerance = 1e-10), "\n")
cat("CDF0_ST's are the same:",
    all.equal(res_R$CDF_ST, res_cpp$CDF_ST, tolerance = 1e-10), "\n")
cat("CDF0_EMP's are the same:",
    all.equal(res_R$CDF_EMP, res_cpp$CDF_EMP, tolerance = 1e-10), "\n")


####____________________________________________________________________________
#### TEST interpolate                                                       ####
l <- 1e1
m <- 1e1

x <- 1:l
x0 <- c(runif(1), x + sort(runif(l)))
# x0 <- (x + sort(runif(l)))[-l]

CDF <- t(apply(matrix(rnorm(l * m), nrow = l, ncol = m), 1, sort))

res_R <- interpolate(x0, x, CDF)
res_CPP <- interpolate_cpp(x0, x, CDF)

cat("Interpolations are the same:",
    all.equal(res_R, res_CPP, tolerance = 1e-12), "\n")

# microbenchmark::microbenchmark(interpolate(x0, x, CDF),
#                                interpolate_cpp(x0, x, CDF),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST SS and CRPS                                                       ####
a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)

r.cond.dist <- function(x)        rgamma(1,     shape = a(x), scale = b(x))
d.cond.dist <- function(x, alpha) dgamma(alpha, shape = a(x), scale = b(x))
p.cond.dist <- function(x, alpha) pgamma(alpha, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))

xx <- seq(1, 4, length.out = 2e2)
yy <- seq(q.cond.dist(1, 0.05), q.cond.dist(4, 0.95), length.out = 2e2)

l0 <- 1e3                 # Number of potential covariates
x0 <- 1 + 3 * (1:l0) / l0 # Potential covariates

n <- 1e3                  # Sample size
X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
Y <- round(Y, 2)          # Creates some ties
W <- rep(1, n)            # Sample weights

delta0 <- 1e-1            # Threshold for estimation precision

res <- dist_reg_cpp(X, Y, W, delta0, x0, TRUE)

SS_CRPS <- SS_CRPS_gamma(x0, l0, res, a, b)

plot(x0, SS_CRPS$SS[,1], type = "l", ylim = range(SS_CRPS$SS), ylab = "Simple score")
lines(x0, SS_CRPS$SS[,2], col = 2)
lines(x0, SS_CRPS$SS[,3], col = 3)

plot(x0, SS_CRPS$CRPS[,1], type = "l", ylim = range(SS_CRPS$CRPS), ylab = "CRPS")
lines(x0, SS_CRPS$CRPS[,2], col = 2)
lines(x0, SS_CRPS$CRPS[,3], col = 3)

# microbenchmark::microbenchmark(
#   dist_reg_cpp(X, Y, W, 1e1, x0, TRUE),
#   dist_reg_cpp(X, Y, W, 1e0, x0, TRUE),
#   dist_reg_cpp(X, Y, W, 1e-1, x0, TRUE),
#   times = 10
# )

# res1 <- dist_reg_cpp(X, Y, W, 10, x0, TRUE)
# res2 <- dist_reg_cpp(X, Y, W, 1, x0, TRUE)
# res3 <- dist_reg_cpp(X, Y, W, 0.1, x0, TRUE)
#
# SS_CRPS1 <- SS_CRPS_gamma(x0, l0, res1, a, b)
# SS_CRPS2 <- SS_CRPS_gamma(x0, l0, res2, a, b)
# SS_CRPS3 <- SS_CRPS_gamma(x0, l0, res3, a, b)
#
# all.equal(SS_CRPS1$SS, SS_CRPS2$SS)
# all.equal(SS_CRPS1$SS, SS_CRPS3$SS)
# all.equal(SS_CRPS2$SS, SS_CRPS3$SS)
#
# all.equal(SS_CRPS1$CRPS, SS_CRPS2$CRPS)
# all.equal(SS_CRPS1$CRPS, SS_CRPS3$CRPS)
# all.equal(SS_CRPS2$CRPS, SS_CRPS3$CRPS)

