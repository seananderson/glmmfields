<!-- README.md is generated from README.Rmd. Please edit that file -->

# glmmfields <img src="inst/logo.png" align="right" />

[![Travis-CI Build
Status](https://travis-ci.org/seananderson/glmmfields.svg?branch=master)](https://travis-ci.org/seananderson/glmmfields)
[![CRAN
status](https://www.r-pkg.org/badges/version/glmmfields)](https://cran.r-project.org/package=glmmfields)
<!-- [![codecov](https://codecov.io/github/seananderson/glmmfields/branch/master/graphs/badge.svg)](https://codecov.io/github/seananderson/glmmfields) -->

The glmmfields R package implements Bayesian spatiotemporal models that
allow for extreme spatial deviations through time. It uses a predictive
process approach with random fields implemented through a multivariate-t
distribution instead of a multivariate normal. The models are fit with
[Stan](http://mc-stan.org/).

We published a paper describing the model and package in *Ecology*:

Anderson, S. C., Ward, E. J. 2019. Black swans in space: modelling
spatiotemporal processes with extremes. 100(1):e02403.
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

![](README-figs/plot-sim-1.png)<!-- -->

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
#> df[1]       3.79    0.04 1.39    2.15    2.74    3.47    4.44    7.26  1169    1
#> gp_sigma    0.30    0.00 0.04    0.22    0.27    0.30    0.33    0.39   464    1
#> gp_theta    2.59    0.00 0.07    2.46    2.54    2.58    2.63    2.72  1189    1
#> sigma[1]    0.10    0.00 0.00    0.09    0.10    0.10    0.10    0.10  1805    1
#> lp__     2291.54    0.38 9.45 2271.16 2285.75 2291.95 2297.94 2308.92   623    1
#> 
#> Samples were drawn using NUTS(diag_e) at Wed Jan 16 15:22:32 2019.
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
#> 1  -0.0294  -0.0894    0.0343
#> 2  -0.290   -0.362    -0.221 
#> 3  -0.397   -0.450    -0.346 
#> 4  -0.195   -0.268    -0.124 
#> 5  -0.0353  -0.109     0.0389
#> 6  -0.215   -0.291    -0.135

# posterior predictive intervals on new observations (include observation error):
p <- predictive_interval(m)
head(p)
#> # A tibble: 6 x 3
#>   estimate conf_low conf_high
#>      <dbl>    <dbl>     <dbl>
#> 1  -0.0294   -0.243   0.181  
#> 2  -0.290    -0.508  -0.0944 
#> 3  -0.397    -0.607  -0.209  
#> 4  -0.195    -0.395   0.00251
#> 5  -0.0353   -0.245   0.179  
#> 6  -0.215    -0.429  -0.0121
```

Use the `tidy` method to extract parameter estimates as a data frame:

``` r
x <- tidy(m, conf.int = TRUE)
head(x)
#> # A tibble: 6 x 5
#>   term                     estimate std.error conf.low conf.high
#>   <chr>                       <dbl>     <dbl>    <dbl>     <dbl>
#> 1 df[1]                      3.79     1.39      2.15      7.26  
#> 2 gp_sigma                   0.300    0.0428    0.218     0.389 
#> 3 gp_theta                   2.59     0.0656    2.46      2.72  
#> 4 sigma[1]                   0.0979   0.00216   0.0935    0.102 
#> 5 spatialEffectsKnots[1,1]  -0.110    0.0366   -0.179    -0.0373
#> 6 spatialEffectsKnots[2,1]  -0.231    0.0341   -0.296    -0.164
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
