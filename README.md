<!-- README.md is generated from README.Rmd. Please edit that file -->
rrfields
========

[![Travis-CI Build Status](https://travis-ci.org/seananderson/rrfields.svg?branch=master)](https://travis-ci.org/seananderson/rrfields) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/seananderson/rrfields?branch=master&svg=true)](https://ci.appveyor.com/project/seananderson/rrfields) <!-- [![codecov](https://codecov.io/github/seananderson/rrfields/branch/master/graphs/badge.svg)](https://codecov.io/github/seananderson/rrfields) -->

The rrfields R package implements spatiotemporal models that allow for extreme spatial deviations through time. It uses a predictive process approach with random fields implemented through a multivariate-t distribution instead of a multivariate normal. The models are fit with [Stan](http://mc-stan.org/).

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
s <- sim_rrfield(df = 2, n_draws = 15, n_knots = 12, gp_scale = 2.5, 
  gp_sigma = 0.2, sd_obs = 0.1)
head(s$dat)
#>   time pt           y      lon      lat station_id
#> 1    1  1 0.562277290 2.655087 6.547239          1
#> 2    1  2 0.226782328 3.721239 3.531973          2
#> 3    1  3 0.024898812 5.728534 2.702601          3
#> 4    1  4 0.004810619 9.082078 9.926841          4
#> 5    1  5 0.440714012 2.016819 6.334933          5
#> 6    1  6 0.172293349 8.983897 2.132081          6
```

``` r
print(s$plot)
```

![](README-figs/plot-sim-1.png)

Fit the model:

``` r
options(mc.cores = parallel::detectCores()) # for parallel processing
m <- rrfield(y ~ 1, data = s$dat, time = "time",
  lat = "lat", lon = "lon", station = "station_id", nknots = 12, iter = 600)
```

``` r
print(m)
#> Inference for Stan model: rrfield.
#> 4 chains, each with iter=600; warmup=300; thin=1; 
#> post-warmup draws per chain=300, total post-warmup draws=1200.
#> 
#>             mean se_mean    sd    2.5%     25%     50%     75%   97.5%
#> df[1]       2.85    0.02  0.72    2.03    2.31    2.66    3.18    4.69
#> gp_sigma    0.22    0.00  0.03    0.16    0.20    0.22    0.24    0.29
#> gp_scale    2.45    0.00  0.07    2.33    2.41    2.45    2.50    2.60
#> B[1]       -0.03    0.01  0.02   -0.07   -0.04   -0.03   -0.01    0.02
#> sigma[1]    0.10    0.00  0.00    0.10    0.10    0.10    0.10    0.11
#> lp__     2773.45    0.55 10.20 2751.60 2766.96 2774.10 2780.34 2792.45
#>          n_eff Rhat
#> df[1]     1200 1.00
#> gp_sigma  1200 1.00
#> gp_scale  1200 1.00
#> B[1]        20 1.21
#> sigma[1]  1200 1.00
#> lp__       347 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Thu Dec  8 14:31:14 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Plot:

``` r
plot(m)
```

![](README-figs/plot-predictions-1.png)

``` r
library(bayesplot)
posterior <- rstan::extract(m$model, inc_warmup = FALSE, permuted = FALSE)
pars <- c("df[1]", "gp_sigma", "sigma[1]", "gp_scale")
mcmc_trace(posterior,  pars = pars)
```

![](README-figs/plot-1.png)

``` r

mm <- as.matrix(m$model)
mcmc_areas(mm, pars = pars)
```

![](README-figs/plot-2.png)

References
==========

Predictive-process models:

Latimer, A. M., S. Banerjee, H. Sang Jr, E. S. Mosher, and J. A. Silander Jr. 2009. Hierarchical models facilitate spatial analysis of large data sets: a case study on invasive plant species in the northeastern United States. Ecology Letters 12:144–154.

Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014. Spatial semiparametric models improve estimates of species abundance and distribution. Canadian Journal of Fisheries and Aquatic Sciences 71:1655–1666.

...
