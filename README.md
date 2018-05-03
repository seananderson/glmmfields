<!-- README.md is generated from README.Rmd. Please edit that file -->
glmmfields
==========

[![Travis-CI Build
Status](https://travis-ci.org/seananderson/glmmfields.svg?branch=master)](https://travis-ci.org/seananderson/glmmfields)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/seananderson/glmmfields?branch=master&svg=true)](https://ci.appveyor.com/project/seananderson/glmmfields)
<!-- [![codecov](https://codecov.io/github/seananderson/glmmfields/branch/master/graphs/badge.svg)](https://codecov.io/github/seananderson/glmmfields) -->

The glmmfields R package implements Bayesian spatiotemporal models that
allow for extreme spatial deviations through time. It uses a predictive
process approach with random fields implemented through a multivariate-t
distribution instead of a multivariate normal. The models are fit with
[Stan](http://mc-stan.org/).

You can install the development version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("seananderson/glmmfields")
```

glmmfields can also fit spatial GLMs with Stan. See the vignette:

``` r
vignette("spatial-glms", package = "glmmfields")
```

An example spatiotemporal model
-------------------------------

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
  nknots = 12, estimate_df = TRUE, iter = 800
)
```

``` r
print(m)
#> Inference for Stan model: glmmfields.
#> 4 chains, each with iter=800; warmup=400; thin=1; 
#> post-warmup draws per chain=400, total post-warmup draws=1600.
#> 
#>             mean se_mean   sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
#> df[1]       3.71    0.04 1.42    2.05    2.67    3.38    4.38    7.40  1060 1.00
#> gp_sigma    0.30    0.00 0.04    0.22    0.27    0.30    0.33    0.39   397 1.01
#> gp_theta    2.59    0.00 0.07    2.46    2.54    2.58    2.63    2.72  1289 1.00
#> sigma[1]    0.10    0.00 0.00    0.09    0.10    0.10    0.10    0.10  1600 1.00
#> lp__     2290.00    0.41 9.70 2269.22 2283.94 2290.55 2296.85 2307.38   571 1.01
#> 
#> Samples were drawn using NUTS(diag_e) at Thu May  3 15:32:47 2018.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
```

Plot:

``` r
plot(m, type = "prediction") + scale_color_gradient2()
```

![](README-figs/plot-predictions-1.png)

``` r
plot(m, type = "spatial-residual")
```

![](README-figs/plot-predictions-2.png)

Predictions:

``` r
# link scale:
p <- predict(m)
head(p)
#> # A tibble: 6 x 3
#>   estimate conf_low conf_high
#>      <dbl>    <dbl>     <dbl>
#> 1  -0.0295  -0.0863    0.0318
#> 2  -0.291   -0.356    -0.226 
#> 3  -0.398   -0.446    -0.345 
#> 4  -0.195   -0.268    -0.120 
#> 5  -0.0363  -0.114     0.0407
#> 6  -0.216   -0.293    -0.140

# posterior predictive intervals on new observations (include observation error):
p <- predictive_interval(m)
head(p)
#> # A tibble: 6 x 3
#>   estimate conf_low conf_high
#>      <dbl>    <dbl>     <dbl>
#> 1  -0.0295   -0.221   0.178  
#> 2  -0.291    -0.498  -0.0925 
#> 3  -0.398    -0.599  -0.182  
#> 4  -0.195    -0.406   0.00577
#> 5  -0.0363   -0.249   0.162  
#> 6  -0.216    -0.425  -0.00767
```

Use the `tidy` method to extract parameter estimates as a data frame:

``` r
x <- tidy(m, conf.int = TRUE)
head(x)
#>                       term   estimate   std.error    conf.low   conf.high
#> 1                    df[1]  3.7093487 1.418096126  2.05397716  7.39855605
#> 2                 gp_sigma  0.3000014 0.044343226  0.22202151  0.39188751
#> 3                 gp_theta  2.5857411 0.067728885  2.45893749  2.71856485
#> 4                 sigma[1]  0.0981333 0.002130289  0.09417861  0.10239067
#> 5 spatialEffectsKnots[1,1] -0.1095841 0.035728942 -0.17816094 -0.03775841
#> 6 spatialEffectsKnots[2,1] -0.2302546 0.037457393 -0.30500031 -0.15813113
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

![](README-figs/grid-predictions-1.png)

References
==========

Anderson, S. A., Ward, E. J. In press. Black swans in space: modelling
spatiotemporal processes with extremes.
[Code](https://github.com/seananderson/spatial-extremes)

Latimer, A. M., S. Banerjee, H. Sang Jr, E. S. Mosher, and J. A.
Silander Jr. 2009. Hierarchical models facilitate spatial analysis of
large data sets: a case study on invasive plant species in the
northeastern United States. Ecology Letters 12:144–154.

Shelton, A. O., J. T. Thorson, E. J. Ward, and B. E. Feist. 2014.
Spatial semiparametric models improve estimates of species abundance and
distribution. Canadian Journal of Fisheries and Aquatic Sciences
71:1655–1666.
