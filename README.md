
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rcea <img src="man/figures/package-sticker.png" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->

[![R CMD
Check](https://github.com/david-beauchesne/rcea/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/david-beauchesne/rcea/actions/workflows/R-CMD-check.yaml)
[![Website](https://github.com/david-beauchesne/rcea/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/david-beauchesne/rcea/actions/workflows/pkgdown.yaml)
[![Test
coverage](https://github.com/david-beauchesne/rcea/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/david-beauchesne/rcea/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/david-beauchesne/rcea/branch/master/graph/badge.svg)](https://codecov.io/gh/david-beauchesne/rcea)
[![CRAN
status](https://www.r-pkg.org/badges/version/rcea)](https://CRAN.R-project.org/package=rcea)
[![License: GPL (&gt;=
2)](https://img.shields.io/badge/License-GPL%20%28%3E%3D%202%29-blue.svg)](https://choosealicense.com/licenses/gpl-2.0/)
[![LifeCycle](https://img.shields.io/badge/lifecycle-experimental-orange)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Project Status:
Concept](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![Dependencies](https://img.shields.io/badge/dependencies-0/0-brightgreen?style=flat)](#)
<!-- badges: end -->

The goal of the R package `rcea` is to **{{ PLEASE ADD A FEW LINES }}**

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("david-beauchesne/rcea")
```

Then you can attach the package `rcea`:

``` r
library("rcea")
```

## Overview

``` r
library(rcea)
library(stars)
#> Loading required package: abind
#> Loading required package: sf
#> Linking to GEOS 3.11.2, GDAL 3.6.3, PROJ 9.2.0; sf_use_s2() is TRUE

# Data
drivers <- rcea:::drivers 
vc <- rcea:::vc
sensitivity <- rcea:::sensitivity
metaweb <- rcea:::metaweb
trophic_sensitivity <- rcea::trophic_sensitivity
pal <- viridis::viridis

# Plots 
plot(merge(drivers), col = pal)
```

<img src="man/figures/unnamed-chunk-4-1.png" width="100%" />

``` r
plot(merge(vc), col = pal(1))
#> Warning in plot.stars(merge(vc), col = pal(1)): breaks="quantile" leads to a
#> single class; maybe try breaks="equal" instead?
```

<img src="man/figures/unnamed-chunk-4-2.png" width="100%" />

``` r
# Cumulative footprint
foot_dr <- cea_extract(drivers, cumul_fun = "footprint")
foot_vc <- cea_extract(vc, cumul_fun = "footprint")
plot(foot_dr, col = pal)
```

<img src="man/figures/unnamed-chunk-4-3.png" width="100%" />

``` r
plot(foot_vc, breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-4.png" width="100%" />

``` r
# Cumulative exposure 
expo <- exposure(drivers, vc, "stars")

# Extract specific attributes and evaluate cumulative exposure
dr_sel <- c("driver1","driver5")
vc_sel <- c("vc4","vc7","vc10","vc12")
dat <- cea_extract(expo, dr_sel = dr_sel, vc_sel = vc_sel) 
plot(dat["vc4"], col = pal) # Exposure of vc4 to driver1 and driver5
```

<img src="man/figures/unnamed-chunk-4-5.png" width="100%" />

``` r
# Cumulative effects assessment (Halpern et al. 2008)
halpern <- cea(drivers, vc, sensitivity, "stars")

# Cumulative effects of all drivers on all vc
dat <- cea_extract(halpern, cumul_fun = "drivers")
plot(merge(dat), breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-6.png" width="100%" />

``` r
# Cumulative effects of all drivers on each vc
dat <- cea_extract(halpern, cumul_fun = "vc") 
plot(merge(dat), breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-7.png" width="100%" />

``` r
# Full cumulative effects
dat <- cea_extract(dat, cumul_fun = "full") 
plot(dat, breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-8.png" width="100%" />

``` r
# Network-scale cumulative effects assessment (Beauchesne et al. 2021)
beauchesne <- ncea(drivers, vc, sensitivity, metaweb, trophic_sensitivity)

# Net cumulative effects
dat <- cea_extract(beauchesne$net, cumul_fun = "full") 
plot(dat, breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-9.png" width="100%" />

``` r
# Direct cumulative effects
dat <- cea_extract(beauchesne$direct, cumul_fun = "full") 
plot(dat, breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-10.png" width="100%" />

``` r
# Indirect cumulative effects
dat <- cea_extract(beauchesne$indirect, cumul_fun = "full") 
plot(dat, breaks = "equal", col = pal)
```

<img src="man/figures/unnamed-chunk-4-11.png" width="100%" />

## Citation

Please cite this package as:

> Beauchesne David (2023) rcea: An R package to perform cumulative
> effects assessments. R package version 0.0.0.9000.

## Code of Conduct

Please note that the `rcea` project is released with a [Contributor Code
of
Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
