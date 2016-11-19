<!-- README.md is generated from README.Rmd. Please edit that file -->
rrfields
========

The rrfields R package implements predictive process spatiotemporal models that allow for extreme spatial deviations through time. It uses random fields implemented with a multivariate-t distribution instead of a multivariate normal.

You can install the development version of the package with:

    ```R
    # install.packages("devtools")
    devtools::install_github("seananderson/rrfields")
    ```

An example model
----------------

Simulate data:

``` r
library(rrfields)
#> Loading required package: Rcpp
#> Warning: package 'Rcpp' was built under R version 3.3.2
set.seed(123)
s <- sim_rrfield(df = 4, n_draws = 15)
```

``` r
library(ggplot2)
#> Warning: package 'ggplot2' was built under R version 3.3.2
ggplot(s$dat, aes(x = lon, y = lat, z = y, colour = y)) +
  facet_wrap(~time, nrow = 3) +
  geom_point(size = 3) +
  scale_color_gradient2()
```

![](README-figs/plot-sim-1.png)

Fit the model:

``` r
options(mc.cores = parallel::detectCores())
m <- rrfield(y ~ 1, data = s$dat, time = "time",
  lat = "lat", lon = "lon", nknots = 15, iter = 400, chains = 4)
```

Plot:

``` r
library(bayesplot)
#> Warning: package 'bayesplot' was built under R version 3.3.2
#> This is bayesplot version 1.0.0
mm <- as.matrix(m)
mcmc_areas(mm, pars = c("df[1]"))
```

![](README-figs/plot-1.png)

``` r
mcmc_areas(mm, pars = c("gp_sigma", "sigma", "gp_scale"))
```

![](README-figs/plot-2.png)

``` r
posterior <- rstan::extract(m, inc_warmup = FALSE, permuted = FALSE)
mcmc_trace(posterior,  pars = c("df[1]", "gp_sigma", "sigma", "gp_scale"))
```

![](README-figs/plot-3.png)
