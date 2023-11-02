
<!-- README.md is generated from README.Rmd. Please edit that file -->

`{r, include = FALSE}d knitr::opts_chunk$set(   collapse = TRUE,   comment = "#>",   fig.path = "man/figures/README-",   out.width = "100%" )`

# mxfda

<!-- badges: start -->
<!-- badges: end -->

The goal of mxfda is to …

## Installation

You can install the development version of mxfda from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("julia-wrobel/mxfda")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(mxfda)
```

    ## Loading required package: spatstat.explore

    ## Loading required package: spatstat.data

    ## Loading required package: spatstat.geom

    ## spatstat.geom 3.2-7

    ## Loading required package: spatstat.random

    ## spatstat.random 3.2-1

    ## Loading required package: nlme

    ## spatstat.explore 3.2-5

``` r
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
