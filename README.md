
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
![(Q_x)\_{x\\in \\mathfrak{X}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28Q_x%29_%7Bx%5Cin%20%5Cmathfrak%7BX%7D "(Q_x)_{x\in \mathfrak{X}").
More precisely

![
Q_x := \\mathrm{Gamma}\\bigl(a(x), b(x)\\bigr),\\quad \\text{for all}\\ x \\in \\mathfrak{X} := \[1,4\],
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AQ_x%20%3A%3D%20%5Cmathrm%7BGamma%7D%5Cbigl%28a%28x%29%2C%20b%28x%29%5Cbigr%29%2C%5Cquad%20%5Ctext%7Bfor%20all%7D%5C%20x%20%5Cin%20%5Cmathfrak%7BX%7D%20%3A%3D%20%5B1%2C4%5D%2C%0A "
Q_x := \mathrm{Gamma}\bigl(a(x), b(x)\bigr),\quad \text{for all}\ x \in \mathfrak{X} := [1,4],
")

with some shape
![a: \\mathfrak{X} \\to (0,\\infty)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a%3A%20%5Cmathfrak%7BX%7D%20%5Cto%20%280%2C%5Cinfty%29 "a: \mathfrak{X} \to (0,\infty)")
and scale
![b: \\mathfrak{X} \\to (0,\\infty)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b%3A%20%5Cmathfrak%7BX%7D%20%5Cto%20%280%2C%5Cinfty%29 "b: \mathfrak{X} \to (0,\infty)").

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
![(\\mathrm{d}P_x/\\mathrm{d}y)(y)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Cmathrm%7Bd%7DP_x%2F%5Cmathrm%7Bd%7Dy%29%28y%29 "(\mathrm{d}P_x/\mathrm{d}y)(y)")
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

<img src="man/figures/README-Visual output-1.png" width="100%" />

``` r
contour(xx, yy, outer(xx, yy, FUN = "d.cond.dist"),
        nlevels = 20,
        xlab = expression(italic(x)), ylab = expression(italic(y)),
        xlim = range(xx), ylim = range(yy))
```

<img src="man/figures/README-Visual output-2.png" width="100%" />

## PART 2: Small data example, first (nonparametric) fit

Let us start with a small sample:
![n = 30](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n%20%3D%2030 "n = 30")
observations with covariates in a set
![\\mathfrak{X}\_o := 1 + 3\*\\{1,2,...,\\ell_o\\}/\\ell_o](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathfrak%7BX%7D_o%20%3A%3D%201%20%2B%203%2A%5C%7B1%2C2%2C...%2C%5Cell_o%5C%7D%2F%5Cell_o "\mathfrak{X}_o := 1 + 3*\{1,2,...,\ell_o\}/\ell_o"),
for
![\\ell_o = 10](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell_o%20%3D%2010 "\ell_o = 10").

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
![x_1 \< x_2 \< ... \< x\_\\ell](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;x_1%20%3C%20x_2%20%3C%20...%20%3C%20x_%5Cell "x_1 < x_2 < ... < x_\ell")
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

<img src="man/figures/README-fit with indices-1.png" width="100%" />

The family of distributions is estimated at each points of this grid.
The estimated conditional distribution functions are given by
`res$CDF.LR`, an
![\\ell](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cell "\ell")-by-![m](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;m "m")
matrix.

    #>             [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
    #>  [1,] 0.13322490 0.2663839 0.3994703 0.5324763 0.6653919 0.7770358 0.8885786
    #>  [2,] 0.06787462 0.1357156 0.2035197 0.2712828 0.3389998 0.4126500 0.4862333
    #>  [3,] 0.06708856 0.1341439 0.2011628 0.2681411 0.3350738 0.4091614 0.4831819
    #>  [4,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [5,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [6,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [7,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [8,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [9,] 0.00000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>            [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
    #>  [1,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [2,] 0.5597367 0.6705905 0.7437475 0.8167914 0.8896848 0.9223035 0.9482503
    #>  [3,] 0.5571218 0.6686340 0.7422255 0.8157033 0.8890296 0.9218420 0.9479430
    #>  [4,] 0.0000000 0.0000000 0.1689034 0.3375456 0.5058402 0.6519566 0.7681861
    #>  [5,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.2050386 0.3681382
    #>  [6,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.1022953
    #>  [7,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [8,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>  [9,] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
    #>           [,15]     [,16]     [,17]      [,18]     [,19]     [,20]     [,21]
    #>  [1,] 1.0000000 1.0000000 1.0000000 1.00000000 1.0000000 1.0000000 1.0000000
    #>  [2,] 0.9741539 1.0000000 1.0000000 1.00000000 1.0000000 1.0000000 1.0000000
    #>  [3,] 0.9740004 1.0000000 1.0000000 1.00000000 1.0000000 1.0000000 1.0000000
    #>  [4,] 0.8842219 1.0000000 1.0000000 1.00000000 1.0000000 1.0000000 1.0000000
    #>  [5,] 0.5309660 0.6934323 1.0000000 1.00000000 1.0000000 1.0000000 1.0000000
    #>  [6,] 0.2044201 0.3063182 0.4985959 0.74932294 1.0000000 1.0000000 1.0000000
    #>  [7,] 0.0000000 0.0000000 0.0000000 0.08106500 0.1621138 0.3254477 0.4887398
    #>  [8,] 0.0000000 0.0000000 0.0000000 0.03555594 0.0711048 0.1427447 0.2143663
    #>  [9,] 0.0000000 0.0000000 0.0000000 0.00000000 0.0000000 0.0000000 0.0000000
    #>           [,22]     [,23]     [,24]     [,25]     [,26]     [,27]     [,28]
    #>  [1,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [2,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [3,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [4,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [5,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [6,] 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000 1.0000000
    #>  [7,] 0.5744934 0.6602162 0.7459007 0.8018720 0.8578433 0.9138146 0.9425431
    #>  [8,] 0.3461405 0.4778675 0.6095354 0.6955444 0.7815533 0.8675622 0.9117082
    #>  [9,] 0.0000000 0.0000000 0.0000000 0.1141816 0.2283633 0.3425449 0.5616966
    #>           [,29] [,30]
    #>  [1,] 1.0000000     1
    #>  [2,] 1.0000000     1
    #>  [3,] 1.0000000     1
    #>  [4,] 1.0000000     1
    #>  [5,] 1.0000000     1
    #>  [6,] 1.0000000     1
    #>  [7,] 0.9712715     1
    #>  [8,] 0.9558541     1
    #>  [9,] 0.7808483     1

One can also retrieve the jumps only.

    #>             [,1]       [,2]       [,3]       [,4]       [,5]       [,6]
    #>  [1,] 0.13322490 0.13315897 0.13308648 0.13300600 0.13291551 0.11164397
    #>  [2,] 0.06787462 0.06784103 0.06780410 0.06776309 0.06771699 0.07365015
    #>  [3,] 0.06708856 0.06705535 0.06701885 0.06697832 0.06693276 0.07408759
    #>  [4,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [5,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [6,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [7,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [8,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [9,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>             [,7]       [,8]      [,9]      [,10]      [,11]      [,12]
    #>  [1,] 0.11154275 0.11142142 0.0000000 0.00000000 0.00000000 0.00000000
    #>  [2,] 0.07358337 0.07350334 0.1108538 0.07315705 0.07304388 0.07289335
    #>  [3,] 0.07402042 0.07393991 0.1115122 0.07359157 0.07347772 0.07332630
    #>  [4,] 0.00000000 0.00000000 0.0000000 0.16890343 0.16864215 0.16829460
    #>  [5,] 0.00000000 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000
    #>  [6,] 0.00000000 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000
    #>  [7,] 0.00000000 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000
    #>  [8,] 0.00000000 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000
    #>  [9,] 0.00000000 0.00000000 0.0000000 0.00000000 0.00000000 0.00000000
    #>            [,13]      [,14]      [,15]      [,16]     [,17]      [,18]
    #>  [1,] 0.00000000 0.00000000 0.00000000 0.00000000 0.0000000 0.00000000
    #>  [2,] 0.03261873 0.02594683 0.02590359 0.02584608 0.0000000 0.00000000
    #>  [3,] 0.03281247 0.02610094 0.02605744 0.02599959 0.0000000 0.00000000
    #>  [4,] 0.14611644 0.11622947 0.11603576 0.11577814 0.0000000 0.00000000
    #>  [5,] 0.20503864 0.16309960 0.16282777 0.16246627 0.3065677 0.00000000
    #>  [6,] 0.00000000 0.10229529 0.10212480 0.10189806 0.1922778 0.25072699
    #>  [7,] 0.00000000 0.00000000 0.00000000 0.00000000 0.0000000 0.08106500
    #>  [8,] 0.00000000 0.00000000 0.00000000 0.00000000 0.0000000 0.03555594
    #>  [9,] 0.00000000 0.00000000 0.00000000 0.00000000 0.0000000 0.00000000
    #>            [,19]      [,20]      [,21]      [,22]      [,23]      [,24]
    #>  [1,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [2,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [3,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [4,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [5,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [6,] 0.25067706 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [7,] 0.08104885 0.16333389 0.16329208 0.08575357 0.08572284 0.08568444
    #>  [8,] 0.03554886 0.07163992 0.07162158 0.13177419 0.13172696 0.13166796
    #>  [9,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>            [,25]      [,26]      [,27]      [,28]      [,29]      [,30]
    #>  [1,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [2,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [3,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [4,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [5,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [6,] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [7,] 0.05597131 0.05597131 0.05597131 0.02872846 0.02872846 0.02872846
    #>  [8,] 0.08600894 0.08600894 0.08600894 0.04414592 0.04414592 0.04414592
    #>  [9,] 0.11418164 0.11418164 0.11418164 0.21915170 0.21915170 0.21915170
    #> [1] 1 1 1 1 1 1 1 1 1
    #>              [,1]        [,2]        [,3]        [,4]        [,5]        [,6]
    #>  [1,] 0.013248061 0.013241505 0.013234297 0.013226293 0.013217295 0.011102025
    #>  [2,] 0.006763408 0.006760061 0.006756381 0.006752295 0.006747701 0.007338915
    #>  [3,] 0.013306609 0.013300023 0.013292784 0.013284744 0.013275707 0.014694825
    #>  [4,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [5,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [6,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [7,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [8,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [9,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>             [,7]        [,8]       [,9]       [,10]       [,11]       [,12]
    #>  [1,] 0.01109196 0.011079894 0.00000000 0.000000000 0.000000000 0.000000000
    #>  [2,] 0.00733226 0.007324285 0.01104610 0.007289779 0.007278503 0.007263503
    #>  [3,] 0.01468150 0.014665532 0.02211777 0.014596441 0.014573861 0.014543827
    #>  [4,] 0.00000000 0.000000000 0.00000000 0.011291109 0.011273642 0.011250409
    #>  [5,] 0.00000000 0.000000000 0.00000000 0.000000000 0.000000000 0.000000000
    #>  [6,] 0.00000000 0.000000000 0.00000000 0.000000000 0.000000000 0.000000000
    #>  [7,] 0.00000000 0.000000000 0.00000000 0.000000000 0.000000000 0.000000000
    #>  [8,] 0.00000000 0.000000000 0.00000000 0.000000000 0.000000000 0.000000000
    #>  [9,] 0.00000000 0.000000000 0.00000000 0.000000000 0.000000000 0.000000000
    #>             [,13]       [,14]       [,15]       [,16]      [,17]       [,18]
    #>  [1,] 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000 0.000000000
    #>  [2,] 0.003250314 0.002585488 0.002581179 0.002575448 0.00000000 0.000000000
    #>  [3,] 0.006508155 0.005176963 0.005168335 0.005156861 0.00000000 0.000000000
    #>  [4,] 0.009767810 0.007769881 0.007756932 0.007739710 0.00000000 0.000000000
    #>  [5,] 0.013602968 0.010820588 0.010802554 0.010778570 0.02033876 0.000000000
    #>  [6,] 0.000000000 0.006758784 0.006747519 0.006732539 0.01270405 0.016565862
    #>  [7,] 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000 0.010569306
    #>  [8,] 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000 0.005976417
    #>  [9,] 0.000000000 0.000000000 0.000000000 0.000000000 0.00000000 0.000000000
    #>             [,19]      [,20]      [,21]      [,22]      [,23]      [,24]
    #>  [1,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [2,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [3,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [4,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [5,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [6,] 0.016562563 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>  [7,] 0.010567201 0.02129558 0.02129013 0.01118061 0.01117660 0.01117159
    #>  [8,] 0.005975227 0.01204159 0.01203851 0.02214925 0.02214132 0.02213140
    #>  [9,] 0.000000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    #>             [,25]       [,26]       [,27]       [,28]       [,29]       [,30]
    #>  [1,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [2,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [3,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [4,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [5,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [6,] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
    #>  [7,] 0.007297575 0.007297575 0.007297575 0.003745636 0.003745636 0.003745636
    #>  [8,] 0.014456807 0.014456807 0.014456807 0.007420263 0.007420263 0.007420263
    #>  [9,] 0.011553030 0.011553030 0.011553030 0.022174023 0.022174023 0.022174023
    #> [1] 0.09944133 0.09964561 0.19834394 0.06684949 0.06634343 0.06607131 0.13038064
    #> [8] 0.16808492 0.10118116
    #> [1] 0.10000000 0.10000000 0.20000000 0.06666667 0.06666667 0.06666667 0.13333333
    #> [8] 0.16666667 0.10000000
    #>  [1] 0.03331808 0.03330159 0.03328346 0.03326333 0.03324070 0.03313576
    #>  [7] 0.03310572 0.03306971 0.03316387 0.03317733 0.03312601 0.03305774
    #> [13] 0.03312925 0.03311170 0.03305652 0.03298313 0.03304280 0.03311159
    #> [19] 0.03310499 0.03333717 0.03332863 0.03332986 0.03331792 0.03330299
    #> [25] 0.03330741 0.03330741 0.03330741 0.03333992 0.03333992 0.03333992
    #>  [1] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
    #>  [7] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
    #> [13] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
    #> [19] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333
    #> [25] 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333 0.03333333

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

<img src="man/figures/README-fit medium-1.png" width="100%" />

``` r
res$tot.time # Less than a tenth of a second on a 7th generation i7 CPU
#> Time difference of 0.08035493 secs
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
#> Time difference of 5.162873 secs
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

<img src="man/figures/README-plot fit-1.png" width="100%" />

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

<img src="man/figures/README-plot differences between truth and estimators-1.png" width="100%" /><img src="man/figures/README-plot differences between truth and estimators-2.png" width="100%" /><img src="man/figures/README-plot differences between truth and estimators-3.png" width="100%" /><img src="man/figures/README-plot differences between truth and estimators-4.png" width="100%" />

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
