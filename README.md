<!-- README.md is generated from README.Rmd. Please edit that file -->
rrfields
========

[![Travis-CI Build Status](https://travis-ci.org/seananderson/rrfields.svg?branch=master)](https://travis-ci.org/seananderson/rrfields)

The rrfields R package implements spatiotemporal models that allow for extreme spatial deviations through time. It uses a predictive process approach with random fields implemented through a multivariate-t distribution instead of a multivariate normal.

You can install the development version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("seananderson/rrfields")
```

An example model
----------------

Simulate data:

``` r
library(rrfields)
set.seed(1)
s <- sim_rrfield(df = 3, n_draws = 12, n_knots = 10, gp_scale = 0.5, 
  gp_sigma = 0.2, sd_obs = 0.1)
head(s$dat)
#>   time pt           y      lon      lat
#> 1    1  1 -0.05743592 2.655087 6.547239
#> 2    1  2 -0.12775526 3.721239 3.531973
#> 3    1  3  0.58242865 5.728534 2.702601
#> 4    1  4 -0.25740965 9.082078 9.926841
#> 5    1  5 -0.10636971 2.016819 6.334933
#> 6    1  6 -0.02949620 8.983897 2.132081
```

``` r
print(s$plot)
```

![](README-figs/plot-sim-1.png)

Fit the model:

``` r
options(mc.cores = parallel::detectCores()) # for parallel processing
m <- rrfield(y ~ 1, data = s$dat, time = "time",
  lat = "lat", lon = "lon", nknots = 10, iter = 400, chains = 4)
```

``` r
pars <- c("df[1]", "gp_sigma", "sigma[1]", "gp_scale")
print(m$model, pars = pars)
#> Inference for Stan model: rrfield.
#> 4 chains, each with iter=400; warmup=200; thin=1; 
#> post-warmup draws per chain=200, total post-warmup draws=800.
#> 
#>          mean se_mean   sd 2.5%  25%  50%  75% 97.5% n_eff Rhat
#> df[1]    4.15    0.07 1.88 2.16 2.86 3.62 4.88  9.17   800    1
#> gp_sigma 0.24    0.00 0.04 0.17 0.21 0.24 0.27  0.32   800    1
#> sigma[1] 0.11    0.00 0.00 0.10 0.10 0.11 0.11  0.11   800    1
#> gp_scale 0.50    0.00 0.02 0.47 0.49 0.50 0.51  0.54   800    1
#> 
#> Samples were drawn using NUTS(diag_e) at Sat Nov 19 19:01:26 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Plot:

``` r
library(bayesplot)
posterior <- rstan::extract(m$model, inc_warmup = FALSE, permuted = FALSE)
mcmc_trace(posterior,  pars = pars)
```

![](README-figs/plot-1.png)

``` r

mm <- as.matrix(m$model)
mcmc_areas(mm, pars = pars[ 1])
```

![](README-figs/plot-2.png)

``` r
mcmc_areas(mm, pars = pars[-1])
```

![](README-figs/plot-3.png)

References
==========

Predictive-process models:

Latimer, A. M., S. Banerjee, H. Sang Jr, E. S. Mosher, and J. A. Silander Jr. 2009. Hierarchical models facilitate spatial analysis of large data sets: a case study on invasive plant species in the northeastern United States. Ecology Letters 12:144–154.

Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014. Spatial semiparametric models improve estimates of species abundance and distribution. Canadian Journal of Fisheries and Aquatic Sciences 71:1655–1666.

...
