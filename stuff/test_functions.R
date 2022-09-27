library(LRDistReg)
my_equal <- function(z1, z2) all.equal(z1, z2, tolerance = 1e-13)

####____________________________________________________________________________
#### SETUP                                                                  ####

n <- 1e4; l0 <- 1e2; m0 <- 1e2
# set.seed(1)
X <- sample(rnorm(l0), n, replace = TRUE)
Y <- sample(rnorm(m0), n, replace = TRUE)
W <- rep(1, n)

####____________________________________________________________________________
#### TEST prepare_data                                                      ####

res_R <- prepare_data_R(X, Y, W)
res_C <- prepare_data_C(X, Y, W)

cat("All parameters are the same:",
    all(c(
      res_R$l       == res_C$l,
      res_R$lL      == res_C$lL + 1,
      res_R$m       == res_C$m,
      res_R$mM      == res_C$mM + 1,
      res_R$n       == res_C$n,
      res_R$PP      == res_C$PP,
      res_R$w       == res_C$w,
      res_R$w_jplus == res_C$w_jplus,
      res_R$w_plusk == res_C$w_plusk,
      res_R$w_ul    == res_C$w_ul,
      res_R$w_ol    == res_C$w_ol,
      res_R$W       == res_C$W,
      res_R$x       == res_C$x,
      res_R$X       == res_C$X,
      res_R$y       == res_C$y,
      res_R$Y       == res_C$Y
    )), "\n")

# microbenchmark::microbenchmark(prepare_data_R(X, Y, W),
#                                prepare_data_C(X, Y, W),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST reparametrize                                                     ####

l <- res_C$l
m <- res_C$m
lL <- res_C$lL
mM <- res_C$mM

lambda <- matrix(rnorm(l * m), nrow = l, ncol = m)

cat("lambda1_to_theta's are the same:",
    my_equal(lambda1_to_theta_R(lambda, l, m, mM + 1),
             lambda1_to_theta_C(lambda, l, m, mM)), "\n")
cat("lambda2_to_theta's are the same:",
    my_equal(lambda2_to_theta_R(lambda, l, m, lL + 1),
             lambda2_to_theta_C(lambda, l, m, lL)), "\n")

# microbenchmark::microbenchmark(lambda1_to_theta_R(lambda, l, m, mM + 1),
#                                lambda1_to_theta_C(lambda, l, m, mM),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST likelihood_functions                                              ####

theta <- matrix(rnorm(l * m), nrow = l, ncol = m)

w_ul <- res_C$w_ul
w_ol <- res_C$w_ol
w <- res_C$w
PP <- res_R$PP

tmp_R <- vg1_R(theta, l, m, n, mM + 1, w_ul)
tmp_C <- vg1_C(theta, l, m, n, mM,     w_ul)
cat("v_tilde1's are the same:", my_equal(tmp_R$v, tmp_C$v), "\n")
cat("g_tilde1's are the same:", my_equal(tmp_R$g, tmp_C$g), "\n")

tmp_R <- vg2_R(theta, l, m, n, lL + 1, w_ol)
tmp_C <- vg2_C(theta, l, m, n, lL,     w_ol)
cat("v_tilde2's are the same:", my_equal(tmp_R$v, tmp_C$v), "\n")
cat("g_tilde2's are the same:", my_equal(tmp_R$g, tmp_C$g), "\n")

cat("f_theta's are the same:", my_equal(ftheta_R(theta, n, w, PP),
                                        ftheta_C(theta, l, n, mM, w)), "\n")

# microbenchmark::microbenchmark(vg1_R(theta, l, m, n, mM + 1, w_ul),
#                                vg1_C(theta, l, m, n, mM, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(vg2_R(theta, l, m, n, lL + 1, w_ol),
#                                vg2_C(theta, l, m, n, lL, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(ftheta_R(theta, n, w, PP),
#                                ftheta_C(theta, l, n, mM, w),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST calibrate                                                         ####

w_jplus <- c(res_C$w_jplus)
w_plusk <- c(res_C$w_plusk)

prec <- 1e-12

cat("calibrate's are the same:",
    my_equal(calibrate_R(theta,        n, w, w_jplus, w_plusk, PP, prec),
             calibrate_C(theta, l, mM, n, w, w_jplus, w_plusk, prec)), "\n")

# microbenchmark::microbenchmark(calibrate_R(theta, n, w, w_jplus, w_plusk, PP, prec),
#                                calibrate_C(theta, l, mM, n, w, w_jplus, w_plusk, prec),
#                                times = 1e2)

####____________________________________________________________________________
#### TEST local_search_functions                                            ####

res_R <- local_search1_R(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ul)
res_C <- local_search1_C(theta, l, m, n, lL, mM, w, w_ul)
cat("Psi1's are the same:",   my_equal(res_R$Psi,   res_C$Psi),   "\n")
cat("delta1's are the same:", my_equal(res_R$delta, res_C$delta), "\n")

res_R <- local_search2_R(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ol)
res_C <- local_search2_C(theta, l, m, n, lL, mM, w, w_ol)
cat("Psi2's are the same:",   my_equal(res_R$Psi,   res_C$Psi),   "\n")
cat("delta2's are the same:", my_equal(res_R$delta, res_C$delta), "\n")

Psi <- res_C$Psi
delta <- res_C$delta

res_R <- simple_step_R(theta, Psi, delta, l, m, n, w, PP)
res_C <- simple_step_C(theta, Psi, delta, l, mM, n, w)
cat("simple_step's are the same:", my_equal(res_R[PP], res_C[PP]), "\n")

# microbenchmark::microbenchmark(local_search1_R(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ul),
#                                local_search1_C(theta, l, m, n, lL, mM, w, w_ul),
#                                times = 1e2)
# microbenchmark::microbenchmark(local_search2_R(theta, l, m, n, mM + 1, lL + 1, PP, w, w_ol),
#                                local_search2_C(theta, l, m, n, lL, mM, w, w_ol),
#                                times = 1e2)
# microbenchmark::microbenchmark(simple_step_R(theta, Psi, delta, l, m, n, w, PP),
#                                simple_step_C(theta, Psi, delta, l, mM, n, w),
#                                times = 1e2)


####____________________________________________________________________________
#### TEST main                                                              ####
delta0 <- 1e-9

res_R <- dist_reg_R(X, Y, W, delta0, NULL, TRUE)
res_C <- dist_reg_C(X, Y, W, delta0, numeric(), TRUE)

cat("theta's are the same:",   my_equal(res_R$theta,   res_C$theta),   "\n")
cat("h_TP2's are the same:",   my_equal(res_R$h_TP2,   res_C$h_TP2),   "\n")
cat("q_LR's are the same:",    my_equal(res_R$q_LR,    res_C$q_LR),    "\n")
cat("CDF_LR's are the same:",  my_equal(res_R$CDF_LR,  res_C$CDF_LR),  "\n")
cat("CDF_ST's are the same:",  my_equal(res_R$CDF_ST,  res_C$CDF_ST),  "\n")
cat("CDF_EMP's are the same:", my_equal(res_R$CDF_EMP, res_C$CDF_EMP), "\n")

####____________________________________________________________________________
#### TEST main with interpolate                                             ####
delta0 <- 1e-9
x0 <- unique(sort(floor(X)))

res_R <- dist_reg_R(X, Y, W, delta0, x0, TRUE)
res_C <- dist_reg_C(X, Y, W, delta0, x0, TRUE)

cat("CDF0_LR's are the same:",  my_equal(res_R$CDF_LR,  res_C$CDF_LR),  "\n")
cat("CDF0_ST's are the same:",  my_equal(res_R$CDF_ST,  res_C$CDF_ST),  "\n")
cat("CDF0_EMP's are the same:", my_equal(res_R$CDF_EMP, res_C$CDF_EMP), "\n")


####____________________________________________________________________________
#### TEST interpolate                                                       ####
l <- 1e1
m <- 1e1

x <- 1:l
x0 <- c(runif(1), x + sort(runif(l)))
# x0 <- (x + sort(runif(l)))[-l]

CDF <- t(apply(matrix(rnorm(l * m), nrow = l, ncol = m), 1, sort))

res_R <- interpolate_R(x0, x, CDF)
res_C <- interpolate_C(x0, x, CDF)

cat("Interpolations are the same:", my_equal(res_R, res_C), "\n")

# microbenchmark::microbenchmark(interpolate_R(x0, x, CDF),
#                                interpolate_C(x0, x, CDF),
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

delta0 <- 1e1             # Threshold for estimation precision

res <- dist_reg_C(X, Y, W, delta0, x0, TRUE)

SS_CRPS <- SS_CRPS_gamma(x0, l0, res, a, b)

plot(x0, SS_CRPS$SS[,1], type = "l", ylim = range(SS_CRPS$SS), ylab = "Simple score")
lines(x0, SS_CRPS$SS[,2], col = 2)
lines(x0, SS_CRPS$SS[,3], col = 3)

plot(x0, SS_CRPS$CRPS[,1], type = "l", ylim = range(SS_CRPS$CRPS), ylab = "CRPS")
lines(x0, SS_CRPS$CRPS[,2], col = 2)
lines(x0, SS_CRPS$CRPS[,3], col = 3)
