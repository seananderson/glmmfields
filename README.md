<!-- README.md is generated from README.Rmd. Please edit that file -->

# glmmfields <img src="inst/logo.png" align="right" />

[![Travis-CI Build
Status](https://travis-ci.org/seananderson/glmmfields.svg?branch=master)](https://travis-ci.org/seananderson/glmmfields)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/seananderson/glmmfields?branch=master&svg=true)](https://ci.appveyor.com/project/seananderson/glmmfields)
[![CRAN
status](https://www.r-pkg.org/badges/version/glmmfields)](https://cran.r-project.org/package=glmmfields)
<!-- [![codecov](https://codecov.io/github/seananderson/glmmfields/branch/master/graphs/badge.svg)](https://codecov.io/github/seananderson/glmmfields) -->

The glmmfields R package implements Bayesian spatiotemporal models that
allow for extreme spatial deviations through time. It uses a predictive
process approach with random fields implemented through a multivariate-t
distribution instead of a multivariate normal. The models are fit with
[Stan](http://mc-stan.org/).

We published a paper describing the model and package in *Ecology*:

Anderson, S. C., Ward, E. J. 2018. Black swans in space: modelling
spatiotemporal processes with extremes. In press at Ecology.
<https://doi.org/10.1002/ecy.2403>

You can install the [CRAN
version](https://cran.r-project.org/package=glmmfields) of the package
with:

``` r
install.packages("glmmfields")
```

If you have a C++ compiler installed, you can install the development
version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("seananderson/glmmfields", build_vignettes = TRUE)
```

glmmfields can also fit spatial GLMs with Stan. See the vignette:

``` r
vignette("spatial-glms", package = "glmmfields")
```

## An example spatiotemporal model

``` r
library(glmmfields)
#> Loading required package: Rcpp
library(ggplot2)
```

Simulate data:

``` r
set.seed(42)
s <- sim_glmmfields(
  df = 2.8, n_draws = 12, n_knots = 12, gp_theta = 2.5,
  gp_sigma = 0.2, sd_obs = 0.1
)
head(s$dat)
#>   time pt           y      lon      lat station_id
#> 1    1  1  0.02818963 9.148060 6.262453          1
#> 2    1  2 -0.21924739 9.370754 2.171577          2
#> 3    1  3 -0.34719485 2.861395 2.165673          3
#> 4    1  4 -0.15785483 8.304476 3.889450          4
#> 5    1  5 -0.04703617 6.417455 9.424557          5
#> 6    1  6 -0.23904924 5.190959 9.626080          6
```

``` r
print(s$plot)
```

![](README-figs/plot-sim-1.png)

Fit the model:

``` r
options(mc.cores = parallel::detectCores()) # for parallel processing
m <- glmmfields(y ~ 0,
  data = s$dat, time = "time",
  lat = "lat", lon = "lon",
  nknots = 12, estimate_df = TRUE, iter = 800, seed = 1
)
```

``` r
print(m)
#> Inference for Stan model: glmmfields.
#> 4 chains, each with iter=800; warmup=400; thin=1; 
#> post-warmup draws per chain=400, total post-warmup draws=1600.
#> 
#>             mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#> df[1]       3.80    0.04 1.48    2.10    2.74    3.46    4.43    7.58  1440 1.00
#> gp_sigma    0.30    0.00 0.04    0.22    0.27    0.30    0.33    0.39   569 1.00
#> gp_theta    2.59    0.00 0.07    2.47    2.54    2.59    2.63    2.72  1600 1.00
#> sigma[1]    0.10    0.00 0.00    0.09    0.10    0.10    0.10    0.10  1600 1.00
#> lp__     2291.58    0.41 9.35 2272.47 2285.72 2291.79 2298.18 2308.96   516 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Thu May  3 15:36:11 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Plot:

``` r
plot(m, type = "prediction") + scale_color_gradient2()
```

![](README-figs/plot-predictions-1.png)<!-- -->

``` r
plot(m, type = "spatial-residual")
```

![](README-figs/plot-predictions-2.png)<!-- -->

Predictions:

``` r
# link scale:
p <- predict(m)
head(p)
#> # A tibble: 6 x 3
#>   estimate conf_low conf_high
#>      <dbl>    <dbl>     <dbl>
#> 1  -0.0287  -0.0871    0.0300
#> 2  -0.290   -0.360    -0.218 
#> 3  -0.396   -0.447    -0.344 
#> 4  -0.197   -0.266    -0.129 
#> 5  -0.0356  -0.108     0.0338
#> 6  -0.215   -0.285    -0.142

# posterior predictive intervals on new observations (include observation error):
p <- predictive_interval(m)
head(p)
#> # A tibble: 6 x 3
#>   estimate conf_low conf_high
#>      <dbl>    <dbl>     <dbl>
#> 1  -0.0287   -0.229   0.180  
#> 2  -0.290    -0.490  -0.0931 
#> 3  -0.396    -0.603  -0.196  
#> 4  -0.197    -0.391   0.00359
#> 5  -0.0356   -0.238   0.166  
#> 6  -0.215    -0.429  -0.00893
```

Use the `tidy` method to extract parameter estimates as a data frame:

``` r
x <- tidy(m, conf.int = TRUE)
head(x)
#>                       term    estimate   std.error    conf.low   conf.high
#> 1                    df[1]  3.80038696 1.484906573  2.09881985  7.57545322
#> 2                 gp_sigma  0.30206707 0.042116979  0.22379764  0.38976618
#> 3                 gp_theta  2.58866483 0.065031189  2.46802568  2.72087920
#> 4                 sigma[1]  0.09797513 0.002108619  0.09381234  0.10216097
#> 5 spatialEffectsKnots[1,1] -0.10977109 0.035134765 -0.17942243 -0.04241248
#> 6 spatialEffectsKnots[2,1] -0.23122316 0.034852007 -0.29929125 -0.16103601
```

Make predictions on a fine-scale spatial grid:

``` r
pred_grid <- expand.grid(
  lat = seq(min(s$dat$lat), max(s$dat$lat), length.out = 25),
  lon = seq(min(s$dat$lon), max(s$dat$lon), length.out = 25),
  time = unique(s$dat$time)
)

pred_grid$prediction <- predict(m,
  newdata = pred_grid, type = "response", iter = 100, estimate_method = "median"
)$estimate

ggplot(pred_grid, aes(lon, lat, fill = prediction)) +
  facet_wrap(~time) +
  geom_raster() +
  scale_fill_gradient2()
```

![](README-figs/grid-predictions-1.png)<!-- -->

# References

Anderson, S. C., Ward, E. J. 2018. Black swans in space: modelling
spatiotemporal processes with extremes. In press at *Ecology*.
<https://doi.org/10.1002/ecy.2403>

Latimer, A. M., S. Banerjee, H. Sang Jr, E. S. Mosher, and J. A.
Silander Jr. 2009. Hierarchical models facilitate spatial analysis of
large data sets: a case study on invasive plant species in the
northeastern United States. Ecology Letters 12:144–154.

Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014.
Spatial semiparametric models improve estimates of species abundance and
distribution. Canadian Journal of Fisheries and Aquatic Sciences
71:1655–1666.
