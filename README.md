
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LRDistReg

<!-- badges: start -->
<!-- badges: end -->

The goal of LRDistReg is to …

## Installation

You can install the development version of LRDistReg from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("AlexandreMoesching/LRDistReg")
```

## PART 1: True model (parametric)

The present demo uses a gamma family of distributions
![(Q_x)\_{x\in \mathfrak{X}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Q_x%29_%7Bx%5Cin%20%5Cmathfrak%7BX%7D "(Q_x)_{x\in \mathfrak{X}").
More precisely

![Q_x := \mathrm{Gamma}\bigl(a(x), b(x)\bigr),\quad \text{for all}\\ x \in \mathfrak{X} := \[1,4\],](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Q_x%20%3A%3D%20%5Cmathrm%7BGamma%7D%5Cbigl%28a%28x%29%2C%20b%28x%29%5Cbigr%29%2C%5Cquad%20%5Ctext%7Bfor%20all%7D%5C%20x%20%5Cin%20%5Cmathfrak%7BX%7D%20%3A%3D%20%5B1%2C4%5D%2C "Q_x := \mathrm{Gamma}\bigl(a(x), b(x)\bigr),\quad \text{for all}\ x \in \mathfrak{X} := [1,4],")

with some shape
![a: \mathfrak{X} \to (0,\infty)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%3A%20%5Cmathfrak%7BX%7D%20%5Cto%20%280%2C%5Cinfty%29 "a: \mathfrak{X} \to (0,\infty)")
and scale
![b: \mathfrak{X} \to (0,\infty)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%3A%20%5Cmathfrak%7BX%7D%20%5Cto%20%280%2C%5Cinfty%29 "b: \mathfrak{X} \to (0,\infty)").

``` r
rm(list = ls())
set.seed(111)
a <- function(x) 2 + (x+1)^2
b <- function(x) 1 - exp(-10*x)
r.cond.dist <- function(x) rgamma(1, shape = a(x), scale = b(x))
p.cond.dist <- function(x, y) pgamma(y, shape = a(x), scale = b(x))
d.cond.dist <- function(x, y) dgamma(y, shape = a(x), scale = b(x))
q.cond.dist <- function(x, alpha) qgamma(alpha, shape = a(x), scale = b(x))
```

Visual output of the true family of distributions. For each
![(x,y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x%2Cy%29 "(x,y)")
in a certain rectangle, the value of
![(\mathrm{d}P_x/\mathrm{d}y)(y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Cmathrm%7Bd%7DP_x%2F%5Cmathrm%7Bd%7Dy%29%28y%29 "(\mathrm{d}P_x/\mathrm{d}y)(y)")
is given by the color scale.

``` r
xx <- seq(1, 4, length.out = 2e2)
yy <- seq(q.cond.dist(1, 0.05), q.cond.dist(4, 0.95), length.out = 2e2)
lattice::levelplot(outer(xx, yy, FUN = "d.cond.dist"), 
                   nlevels = 1e2, aspect = "fill",
                   col.regions = hcl.colors(1e2),
                   row.values = xx, column.values = yy,
                   xlab = expression(italic(x)), ylab = expression(italic(y)),
                   xlim = range(xx), ylim = range(yy))
```

<img src="man/figures/README-visual_output-1.png" width="100%" />

``` r

contour(xx, yy, outer(xx, yy, FUN = "d.cond.dist"),
        nlevels = 20,
        xlab = expression(italic(x)), ylab = expression(italic(y)),
        xlim = range(xx), ylim = range(yy))
```

<img src="man/figures/README-visual_output-2.png" width="100%" />

## PART 2: Small data example, first (nonparametric) fit

Let us start with a small sample:
![n = 30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3D%2030 "n = 30")
observations with covariates in a set
![\mathfrak{X}\_o := 1 + 3\*\\{1,2,...,\ell_o\\}/\ell_o](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathfrak%7BX%7D_o%20%3A%3D%201%20%2B%203%2A%5C%7B1%2C2%2C...%2C%5Cell_o%5C%7D%2F%5Cell_o "\mathfrak{X}_o := 1 + 3*\{1,2,...,\ell_o\}/\ell_o"),
for
![\ell_o = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell_o%20%3D%2010 "\ell_o = 10").

``` r
n <- 30
ell0 <- 10
x0 <- 1+(1:ell0)/ell0*3
```

Generate observation pairs
![(X_1,Y_1),(X_2,Y_2),...,(X_N,Y_N)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28X_1%2CY_1%29%2C%28X_2%2CY_2%29%2C...%2C%28X_N%2CY_N%29 "(X_1,Y_1),(X_2,Y_2),...,(X_N,Y_N)").

``` r
X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
```

We now estimate the family of distribution. We use the function
`dist.reg()` with default options, except one plotting option.

``` r
library(LRDistReg)
res <- dist.reg(X, Y, show.design = TRUE)
```

<img src="man/figures/README-fit-1.png" width="100%" />

The design of the experiment should be shown. The color of a pair
![(x,y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28x%2Cy%29 "(x,y)")
is equal to the number of observations at that location, plus 1. In
consequence, black points contain no observations, red points contain
one observation pair, green points contain two observation pairs, etc.

If `indices = TRUE`, then the values
![1,2,...,ell](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2C2%2C...%2Cell "1,2,...,ell")
are used for the plot instead of the unique elements
![x_1 \< x_2 \< ... \< x\_\ell](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%20%3C%20x_2%20%3C%20...%20%3C%20x_%5Cell "x_1 < x_2 < ... < x_\ell")
of
![\\{X_1,X_2,...,X_n\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7BX_1%2CX_2%2C...%2CX_n%5C%7D "\{X_1,X_2,...,X_n\}"),
and the values
![1,2,...,m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;1%2C2%2C...%2Cm "1,2,...,m")
are used instead of the unique elements
![y_1 \< y_2 \< ... \< y_m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;y_1%20%3C%20y_2%20%3C%20...%20%3C%20y_m "y_1 < y_2 < ... < y_m")
of
![\\{Y_1,Y_2,...,Y_n\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7BY_1%2CY_2%2C...%2CY_n%5C%7D "\{Y_1,Y_2,...,Y_n\}").
This improves readability of the design plot.

``` r
res <- dist.reg(X, Y, show.design = TRUE, indices = TRUE)
```

<img src="man/figures/README-fit_with_indices-1.png" width="100%" />

The family of distributions is estimated at each points of this grid.
The estimated conditional distribution functions are given by
`res$CDF.LR`, an
![\ell](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell "\ell")-by-![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
matrix.

``` r
res$CDF.LR
```

One can also retrieve the jumps only.

``` r
res$q.LR
rowSums(res$q.LR)

res$h.TP2
rowSums(res$h.TP2); res$D.Setup$w_j.plus/res$D.Setup$n
colSums(res$h.TP2); res$D.Setup$w_plus.k/res$D.Setup$n
```

To see how the stopping criteria is reached along with the corresponding
values of the negative log-likelihood, the echo option can be activated.

``` r
res <- dist.reg(X, Y, echo = TRUE)
#> iteration 0, delta = Inf, f(theta) = 164.9943
#> iteration 1, delta = 0.7791124, f(theta) = 162.3362
#> iteration 2, delta = 0.3810281, f(theta) = 161.9416
#> iteration 3, delta = 0.09050999, f(theta) = 161.8116
#> iteration 4, delta = 0.03770548, f(theta) = 161.7677
#> iteration 5, delta = 0.009150678, f(theta) = 161.7546
res$delta0
#> [1] 0.01
```

## PART 3: Medium data example, specification of precision parameter

Options:

``` r
n <- 1e2
ell0 <- 10
x0 <- 1+(1:ell0)/ell0*3

X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
Y <- round(Y, 1) # Should create some ties
```

This time we obtain the fit using a self-specified value of the
precision parameter epsilon. It is normally automatically specified via
some kind of rule of thumb.

``` r
res <- dist.reg(X, Y, show.design = TRUE, delta0 = 1e-4)
```

<img src="man/figures/README-fit_medium-1.png" width="100%" />

``` r
res$tot.time # Less than a tenth of a second on a 7th generation i7 CPU
#> Time difference of 0.05826712 secs
```

## PART 4: A larger data example, comparison between Likelihood-Ratio ordering, usual STochastic ordering and the EMPirical

``` r
n <- 1e3
ell0 <- 1e2
x0 <- 1+(1:ell0)/ell0*3

X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])
```

We let the program advise us a “reasonable” value for epsilon.

``` r
res <- dist.reg(X, Y, suggest.delta0 = TRUE, IDR = TRUE, echo = TRUE)
#> iteration 0, delta = Inf, f(theta) = 12082.71
#> iteration 1, delta = 24581.73, f(theta) = 11888.85
#> iteration 2, delta = 1702.676, f(theta) = 11811.1
#> iteration 3, delta = 5127.211, f(theta) = 11775.15
#> iteration 4, delta = 372.9345, f(theta) = 11754.68
#> iteration 5, delta = 614.9382, f(theta) = 11745.79
#> iteration 6, delta = 60.72897, f(theta) = 11740.94
#> iteration 7, delta = 81.10216, f(theta) = 11738.33
#> iteration 8, delta = 20.47268, f(theta) = 11731.7
#> iteration 9, delta = 31.94151, f(theta) = 11728.15
#> iteration 10, delta = 16.38659, f(theta) = 11726.98
#> iteration 11, delta = 11.20272, f(theta) = 11724.46
#> iteration 12, delta = 5.709401, f(theta) = 11721.37
res$delta0
#> [1] 10
res$tot.time
#> Time difference of 4.874189 secs
```

Retrieve all CDF’s and some useful parameters.

``` r
x <- res$D.Setup$x
y <- res$D.Setup$y
ell <- res$D.Setup$ell
m <- res$D.Setup$m

CDF.LR <- res$CDF.LR
CDF.ST <- res$CDF.ST
CDF.EMP <- res$CDF.EMP
CDF.TRUE <- outer(x, y, p.cond.dist)
```

Evaluate estimation quality.

``` r
DIFF.LR <- sum(abs(CDF.LR - CDF.TRUE))/(ell*m)
DIFF.ST <- sum(abs(CDF.ST - CDF.TRUE))/(ell*m)
DIFF.EMP <- sum(abs(CDF.EMP - CDF.TRUE))/(ell*m)
c(DIFF.LR, DIFF.ST, DIFF.EMP)
#> [1] 0.01936299 0.02352673 0.05607067
```

Plot the true CDF
![F_x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;F_x "F_x")
as well as its estimators for two values of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x"),
one middle covariate
(![x = 2.5](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%202.5 "x = 2.5"))
and one boundary covariate
(![x = 4](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x%20%3D%204 "x = 4")).

``` r
par(mfrow = c(2,1), mar = c(4.1,4.1,0.2,0.2))
xx <- c(2.5, 4)
for (xj in xx) {
  j <- which(xj == x)
  plot(y, CDF.TRUE[j,], type = "l",
       xlab = expression(italic(y)),
       ylab = expression(italic(f[x](y))),
       lwd = 2)
  lines(y, CDF.EMP[j,], lty = 2)
  lines(y, CDF.ST[j,], col = "blue", lwd = 2)
  lines(y, CDF.LR[j,], col = "red", lwd = 2)
  legend("topleft", 
         legend = c("True", "EMP", "ST", "LR"),
         col = c(1,1,4,2), lty = c(1,2,1,1), lwd = c(2,1,2,2))
}
```

<img src="man/figures/README-plot_fit-1.png" width="100%" />

The LR-estimator is in general smoother than the ST-estimator. In this
specific example, the LR-estimator is also closer to the true
distribution than the ST-estimator. This was confirmed already earlier
when looking at average absolute differences. Try deactivating the
set.seed to see other outputs.

Plot differences between the true CDF and LR/ST-estimators for all
values of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x").
N.b.: it yields 4 pages of plots.

``` r
n.plot <- min(ceiling(sqrt(ell)), 5)
par(mfrow = c(n.plot, n.plot), mar = c(2,2,1,1))
unif.ylim <- range(c(CDF.LR - CDF.TRUE, CDF.ST - CDF.TRUE))
for (j in 1:ell) {
  plot(y, CDF.ST[j,] - CDF.TRUE[j,], type = "l", col = "blue", 
       lwd = 1, ylim = unif.ylim)
  abline(h = 0)
  lines(y, CDF.LR[j,] - CDF.TRUE[j,], col = "red", lwd = 1)
}
```

## PART 5: Interpolation feature

In case the estimators are desired for other values of
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")
than those provided by the set of covariates, it is still possible to
obtain the estimated CDFs at these other
![x](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x "x")’s
using linear interpolation

``` r
n <- 50
ell0 <- 50
x0 <- 1+(1:ell0)/ell0*3

X <- sort(sample(x0, size = n, replace = TRUE))
Y <- rep(0, n)
for (i in 1:n) Y[i] <- r.cond.dist(X[i])

res <- dist.reg(X, Y, suggest.delta0 = TRUE, IDR = TRUE)
dim(res$CDF.LR); dim(res$CDF.ST); dim(res$CDF.EMP)
#> [1] 28 50
#> [1] 28 50
#> [1] 28 50
res <- dist.reg(X, Y, suggest.delta0 = TRUE, IDR = TRUE, x0 = x0)
dim(res$CDF.LR); dim(res$CDF.ST); dim(res$CDF.EMP)
#> [1] 50 50
#> [1] 50 50
#> [1] 50 50
```
