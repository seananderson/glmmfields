<!-- README.md is generated from README.Rmd. Please edit that file -->
rrfields
========

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
set.seed(999)
s <- sim_rrfield(df = 3, n_draws = 15)
```

``` r
library(ggplot2)
ggplot(s$dat, aes(x = lon, y = lat, colour = y)) +
  facet_wrap(~time, nrow = 3) +
  geom_point(size = 2) +
  scale_color_gradient2()
```

![](README-figs/plot-sim-1.png)

Fit the model:

``` r
options(mc.cores = parallel::detectCores())
m <- rrfield(y ~ 1, data = s$dat, time = "time",
  lat = "lat", lon = "lon", nknots = 15, iter = 400, chains = 4)
```

``` r
print(m, pars = c("df[1]", "gp_sigma", "sigma", "gp_scale", "lp__"))
#> Inference for Stan model: rrfield.
#> 4 chains, each with iter=400; warmup=200; thin=1; 
#> post-warmup draws per chain=200, total post-warmup draws=800.
#> 
#>             mean se_mean    sd    2.5%     25%     50%     75%   97.5%
#> df[1]       5.73    0.12  2.47    2.33    3.93    5.29    6.94   11.77
#> gp_sigma    0.18    0.00  0.02    0.14    0.17    0.18    0.19    0.22
#> sigma       0.10    0.00  0.00    0.10    0.10    0.10    0.10    0.11
#> gp_scale    0.49    0.00  0.02    0.44    0.47    0.49    0.50    0.53
#> lp__     2824.67    0.68 11.97 2801.19 2816.97 2824.89 2833.28 2847.63
#>          n_eff Rhat
#> df[1]      418 1.01
#> gp_sigma   800 1.00
#> sigma      800 1.00
#> gp_scale   800 1.00
#> lp__       313 1.00
#> 
#> Samples were drawn using NUTS(diag_e) at Sat Nov 19 11:50:32 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Plot:

``` r
library(bayesplot)
posterior <- rstan::extract(m, inc_warmup = FALSE, permuted = FALSE)
mcmc_trace(posterior,  pars = c("df[1]", "gp_sigma", "sigma", "gp_scale"))
```

![](README-figs/plot-1.png)

``` r

mm <- as.matrix(m)
mcmc_areas(mm, pars = c("df[1]"))
```

![](README-figs/plot-2.png)

``` r
mcmc_areas(mm, pars = c("gp_sigma", "sigma", "gp_scale"))
```

![](README-figs/plot-3.png)

References
==========

Predictive-process models:

Latimer, A. M., S. Banerjee, H. Sang Jr, E. S. Mosher, and J. A. Silander Jr. 2009. Hierarchical models facilitate spatial analysis of large data sets: a case study on invasive plant species in the northeastern United States. Ecology Letters 12:144–154.

Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014. Spatial semiparametric models improve estimates of species abundance and distribution. Canadian Journal of Fisheries and Aquatic Sciences 71:1655–1666.

...
