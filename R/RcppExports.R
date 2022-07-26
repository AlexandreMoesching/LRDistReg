# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Approximate calibration of log-parameter, C++ version
#'
#' @description Computes an optimal additive adjustment of theta, i.e. a matrix
#' theta.new with values
#'     theta.new_jk = theta_jk + x_j + y_k
#' such that f(theta.new) is minimal.
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param mM (m_j,M_j) index pairs
#' @param n Sample size
#' @param w Sample weights
#' @param w_jplus Row-sums of sample weights
#' @param w_plusk Column sums of sample weights
#' @param prec Precision for the approximate calibration
#'
#' @return Approximately calibrated log-parameter
#'
#' @export
calibrate_C <- function(theta, l, mM, n, w, w_jplus, w_plusk, prec) {
    .Call(`_LRDistReg_calibrate_C`, theta, l, mM, n, w, w_jplus, w_plusk, prec)
}

#' Linear interpolation of the cdf's, C++ version
#'
#' @param x0 Set of covariates on which to extend CDF
#' @param x Set of covariates of CDF
#' @param CDF Step-function matrix, its j-th row contains the cdf for X = x_j
#'
#' @return Linear interpolation of the cdf's on the new set x0
#'
#' @export
interpolate_C <- function(x0, x, CDF) {
    .Call(`_LRDistReg_interpolate_C`, x0, x, CDF)
}

#' v-tilde and gamma-tilde functions (row), C++ version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param w_ul Cumulative row-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#'
#' @export
vg1_C <- function(theta, l, m, n, mM, w_ul) {
    .Call(`_LRDistReg_vg1_C`, theta, l, m, n, mM, w_ul)
}

#' v-tilde and gamma-tilde functions (col), C++ version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param lL (l_k,L_k) index pairs
#' @param w_ol Cumulative column-sums of sample weights
#'
#' @return v-tilde and gamma-tilde functions
#'
#' @export
vg2_C <- function(theta, l, m, n, lL, w_ol) {
    .Call(`_LRDistReg_vg2_C`, theta, l, m, n, lL, w_ol)
}

#' Negative log-likelihood in terms of log-parameter, C++ version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param n Sample size
#' @param mM (m_j,M_j) index pairs
#' @param w Sample weights
#'
#' @return Negative log-likelihood in terms of log-parameter
#'
#' @export
ftheta_C <- function(theta, l, n, mM, w) {
    .Call(`_LRDistReg_ftheta_C`, theta, l, n, mM, w)
}

#' Function to take a simple step, C++ version
#'
#' @param theta Log-parameter
#' @param Psi Proposal
#' @param delta Delta
#' @param l Number of unique covariates
#' @param mM (m_j,M_j) index pairs
#' @param n Sample size
#' @param w Sample weights
#'
#' @return Updated theta parameter
#'
#' @export
simple_step_C <- function(theta, Psi, delta, l, mM, n, w) {
    .Call(`_LRDistReg_simple_step_C`, theta, Psi, delta, l, mM, n, w)
}

#' Local search (row), C++ version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param lL (l_k,L_k) index pairs
#' @param mM (m_j,M_j) index pairs
#' @param w Sample weights
#' @param w_ul Cumulative row-sums of sample weights
#'
#' @return New proposal Psi and step-size delta
#'
#' @export
local_search1_C <- function(theta, l, m, n, lL, mM, w, w_ul) {
    .Call(`_LRDistReg_local_search1_C`, theta, l, m, n, lL, mM, w, w_ul)
}

#' Local search (col), C++ version
#'
#' @param theta Log-parameter
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param n Sample size
#' @param lL (l_k,L_k) index pairs
#' @param mM (m_j,M_j) index pairs
#' @param w Sample weights
#' @param w_ol Cumulative column-sums of sample weights
#'
#' @return New proposal Psi and step-size delta
#'
#' @export
local_search2_C <- function(theta, l, m, n, lL, mM, w, w_ol) {
    .Call(`_LRDistReg_local_search2_C`, theta, l, m, n, lL, mM, w, w_ol)
}

#' Isotonic distributional regression (LR, ST, EMP), C++ version
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified sample weights
#' @param delta0 Threshhold
#' @param x0 Set of covariates on which to estimate the distributions
#' @param ST Boolean indicating whether or not the classical isotonic
#' distributional regression will also be computed
#'
#' @return Isotonic distributional regression(s) and estimation parameters
#'
#' @export
dist_reg_C <- function(X, Y, W, delta0, x0, ST = FALSE) {
    .Call(`_LRDistReg_dist_reg_C`, X, Y, W, delta0, x0, ST)
}

#' Prepare the data, C++ version
#'
#' @param X Covariates
#' @param Y Responses
#' @param W User-specified weights
#'
#' @return A list of pre-computed parameters necessary for the estimation
#'
#' @export
prepare_data_C <- function(X, Y, W) {
    .Call(`_LRDistReg_prepare_data_C`, X, Y, W)
}

#' Transforms lambda (row) into theta, C++ version
#'
#' @param lambda Row-wise differences
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param mM (m_j,M_j) index pairs
#'
#' @return Transformed parameter
#'
#' @export
lambda1_to_theta_C <- function(lambda, l, m, mM) {
    .Call(`_LRDistReg_lambda1_to_theta_C`, lambda, l, m, mM)
}

#' Transforms lambda (col) into theta, C++ version
#'
#' @param lambda Column-wise differences
#' @param l Number of unique covariates
#' @param m Number of unique responses
#' @param lL (l_k,L_k) index pairs
#'
#' @return Transformed parameter
#'
#' @export
lambda2_to_theta_C <- function(lambda, l, m, lL) {
    .Call(`_LRDistReg_lambda2_to_theta_C`, lambda, l, m, lL)
}

