<!-- README.md is generated from README.Rmd. Please edit that file -->
rrfields
========

[![Travis-CI Build Status](https://travis-ci.org/seananderson/rrfields.svg?branch=master)](https://travis-ci.org/seananderson/rrfields) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/seananderson/rrfields?branch=master&svg=true)](https://ci.appveyor.com/project/seananderson/rrfields) [![codecov](https://codecov.io/github/seananderson/rrfields/branch/master/graphs/badge.svg)](https://codecov.io/github/seananderson/rrfields)

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
print(m, pars = pars)
#> $model
#> Inference for Stan model: rrfield.
#> 4 chains, each with iter=400; warmup=200; thin=1; 
#> post-warmup draws per chain=200, total post-warmup draws=800.
#> 
#>                               mean se_mean   sd    2.5%     25%     50%
#> df[1]                         4.15    0.07 1.88    2.16    2.86    3.62
#> gp_sigma                      0.24    0.00 0.04    0.17    0.21    0.24
#> gp_scale                      0.50    0.00 0.02    0.47    0.49    0.50
#> B[1]                         -0.01    0.00 0.01   -0.02   -0.01   -0.01
#> spatialEffectsKnots[1,1]     -0.16    0.00 0.05   -0.26   -0.20   -0.16
#> spatialEffectsKnots[1,2]     -0.10    0.00 0.04   -0.18   -0.13   -0.09
#> spatialEffectsKnots[1,3]      0.64    0.00 0.06    0.52    0.60    0.64
#> spatialEffectsKnots[1,4]     -0.25    0.00 0.05   -0.34   -0.28   -0.25
#> spatialEffectsKnots[1,5]     -0.04    0.00 0.05   -0.14   -0.07   -0.04
#> spatialEffectsKnots[1,6]     -0.05    0.00 0.05   -0.14   -0.08   -0.05
#> spatialEffectsKnots[1,7]     -0.15    0.00 0.05   -0.25   -0.19   -0.15
#> spatialEffectsKnots[1,8]     -0.08    0.00 0.05   -0.19   -0.12   -0.08
#> spatialEffectsKnots[1,9]      0.11    0.00 0.05    0.01    0.07    0.11
#> spatialEffectsKnots[1,10]    -0.15    0.00 0.06   -0.27   -0.19   -0.15
#> spatialEffectsKnots[2,1]     -0.05    0.00 0.05   -0.15   -0.09   -0.06
#> spatialEffectsKnots[2,2]      0.26    0.00 0.04    0.18    0.23    0.26
#> spatialEffectsKnots[2,3]      0.01    0.00 0.05   -0.10   -0.03    0.01
#> spatialEffectsKnots[2,4]     -0.14    0.00 0.05   -0.23   -0.17   -0.14
#> spatialEffectsKnots[2,5]     -0.02    0.00 0.05   -0.11   -0.05   -0.02
#> spatialEffectsKnots[2,6]      0.14    0.00 0.04    0.05    0.11    0.14
#> spatialEffectsKnots[2,7]      0.01    0.00 0.05   -0.09   -0.02    0.01
#> spatialEffectsKnots[2,8]     -0.05    0.00 0.05   -0.13   -0.08   -0.05
#> spatialEffectsKnots[2,9]     -0.17    0.00 0.05   -0.26   -0.20   -0.17
#> spatialEffectsKnots[2,10]    -0.01    0.00 0.06   -0.13   -0.05   -0.02
#> spatialEffectsKnots[3,1]     -0.04    0.00 0.05   -0.13   -0.07   -0.04
#> spatialEffectsKnots[3,2]     -0.07    0.00 0.04   -0.15   -0.10   -0.07
#> spatialEffectsKnots[3,3]      0.15    0.00 0.06    0.05    0.12    0.15
#> spatialEffectsKnots[3,4]     -0.18    0.00 0.05   -0.27   -0.21   -0.18
#> spatialEffectsKnots[3,5]      0.04    0.00 0.05   -0.06    0.00    0.04
#> spatialEffectsKnots[3,6]     -0.15    0.00 0.04   -0.24   -0.17   -0.15
#> spatialEffectsKnots[3,7]      0.00    0.00 0.05   -0.09   -0.04    0.00
#> spatialEffectsKnots[3,8]     -0.03    0.00 0.05   -0.13   -0.06   -0.03
#> spatialEffectsKnots[3,9]     -0.03    0.00 0.05   -0.14   -0.06   -0.03
#> spatialEffectsKnots[3,10]     0.03    0.00 0.06   -0.08   -0.01    0.04
#> spatialEffectsKnots[4,1]     -0.56    0.00 0.05   -0.66   -0.60   -0.56
#> spatialEffectsKnots[4,2]      0.36    0.00 0.04    0.28    0.33    0.36
#> spatialEffectsKnots[4,3]     -0.47    0.00 0.06   -0.59   -0.51   -0.47
#> spatialEffectsKnots[4,4]     -0.08    0.00 0.05   -0.17   -0.12   -0.08
#> spatialEffectsKnots[4,5]     -0.29    0.00 0.05   -0.39   -0.32   -0.29
#> spatialEffectsKnots[4,6]     -0.23    0.00 0.04   -0.32   -0.26   -0.23
#> spatialEffectsKnots[4,7]      0.65    0.00 0.05    0.55    0.61    0.65
#> spatialEffectsKnots[4,8]     -0.03    0.00 0.05   -0.13   -0.07   -0.03
#> spatialEffectsKnots[4,9]     -0.35    0.00 0.05   -0.44   -0.38   -0.35
#> spatialEffectsKnots[4,10]    -0.44    0.00 0.06   -0.57   -0.49   -0.44
#> spatialEffectsKnots[5,1]      0.09    0.00 0.05   -0.01    0.05    0.09
#> spatialEffectsKnots[5,2]     -0.02    0.00 0.04   -0.10   -0.05   -0.02
#> spatialEffectsKnots[5,3]     -0.04    0.00 0.05   -0.14   -0.08   -0.04
#> spatialEffectsKnots[5,4]     -0.07    0.00 0.05   -0.16   -0.10   -0.07
#> spatialEffectsKnots[5,5]     -0.27    0.00 0.05   -0.37   -0.30   -0.27
#> spatialEffectsKnots[5,6]     -0.22    0.00 0.05   -0.32   -0.25   -0.22
#> spatialEffectsKnots[5,7]      0.08    0.00 0.05   -0.02    0.04    0.08
#> spatialEffectsKnots[5,8]     -0.04    0.00 0.05   -0.14   -0.08   -0.05
#> spatialEffectsKnots[5,9]     -0.23    0.00 0.05   -0.32   -0.27   -0.23
#> spatialEffectsKnots[5,10]     0.27    0.00 0.06    0.14    0.23    0.27
#> spatialEffectsKnots[6,1]      0.12    0.00 0.05    0.02    0.08    0.12
#> spatialEffectsKnots[6,2]     -0.02    0.00 0.04   -0.10   -0.05   -0.02
#> spatialEffectsKnots[6,3]      0.25    0.00 0.06    0.14    0.21    0.25
#> spatialEffectsKnots[6,4]      0.20    0.00 0.05    0.11    0.17    0.20
#> spatialEffectsKnots[6,5]     -0.18    0.00 0.06   -0.28   -0.22   -0.18
#> spatialEffectsKnots[6,6]      0.54    0.00 0.04    0.45    0.51    0.54
#> spatialEffectsKnots[6,7]      0.07    0.00 0.06   -0.05    0.02    0.07
#> spatialEffectsKnots[6,8]     -0.36    0.00 0.05   -0.46   -0.39   -0.36
#> spatialEffectsKnots[6,9]     -0.02    0.00 0.05   -0.12   -0.05   -0.02
#> spatialEffectsKnots[6,10]     0.10    0.00 0.06   -0.02    0.06    0.10
#> spatialEffectsKnots[7,1]      1.13    0.00 0.05    1.03    1.10    1.14
#> spatialEffectsKnots[7,2]      0.04    0.00 0.04   -0.05    0.01    0.04
#> spatialEffectsKnots[7,3]      0.23    0.00 0.06    0.13    0.19    0.23
#> spatialEffectsKnots[7,4]     -0.11    0.00 0.05   -0.20   -0.14   -0.11
#> spatialEffectsKnots[7,5]     -0.15    0.00 0.05   -0.26   -0.19   -0.15
#> spatialEffectsKnots[7,6]     -0.05    0.00 0.05   -0.14   -0.08   -0.05
#> spatialEffectsKnots[7,7]      0.44    0.00 0.06    0.33    0.41    0.44
#> spatialEffectsKnots[7,8]      1.02    0.00 0.05    0.91    0.98    1.02
#> spatialEffectsKnots[7,9]      0.48    0.00 0.05    0.39    0.45    0.48
#> spatialEffectsKnots[7,10]     0.66    0.00 0.06    0.54    0.62    0.66
#> spatialEffectsKnots[8,1]     -0.22    0.00 0.05   -0.32   -0.25   -0.22
#> spatialEffectsKnots[8,2]      0.23    0.00 0.04    0.15    0.20    0.23
#> spatialEffectsKnots[8,3]      0.16    0.00 0.06    0.04    0.12    0.16
#> spatialEffectsKnots[8,4]     -0.27    0.00 0.05   -0.36   -0.30   -0.27
#> spatialEffectsKnots[8,5]      0.15    0.00 0.05    0.05    0.12    0.15
#> spatialEffectsKnots[8,6]      0.01    0.00 0.05   -0.07   -0.02    0.01
#> spatialEffectsKnots[8,7]      0.27    0.00 0.05    0.17    0.24    0.27
#> spatialEffectsKnots[8,8]     -0.12    0.00 0.05   -0.23   -0.16   -0.12
#> spatialEffectsKnots[8,9]     -0.15    0.00 0.05   -0.24   -0.18   -0.15
#> spatialEffectsKnots[8,10]    -0.16    0.00 0.06   -0.28   -0.20   -0.16
#> spatialEffectsKnots[9,1]     -0.04    0.00 0.05   -0.14   -0.08   -0.04
#> spatialEffectsKnots[9,2]      0.17    0.00 0.04    0.10    0.15    0.17
#> spatialEffectsKnots[9,3]     -0.25    0.00 0.06   -0.36   -0.29   -0.25
#> spatialEffectsKnots[9,4]      0.40    0.00 0.05    0.31    0.37    0.40
#> spatialEffectsKnots[9,5]     -0.51    0.00 0.05   -0.61   -0.55   -0.51
#> spatialEffectsKnots[9,6]     -0.40    0.00 0.05   -0.49   -0.43   -0.40
#> spatialEffectsKnots[9,7]      0.58    0.00 0.05    0.47    0.54    0.57
#> spatialEffectsKnots[9,8]     -0.39    0.00 0.05   -0.50   -0.42   -0.39
#> spatialEffectsKnots[9,9]      0.18    0.00 0.05    0.08    0.14    0.18
#> spatialEffectsKnots[9,10]    -0.18    0.00 0.06   -0.30   -0.22   -0.18
#> spatialEffectsKnots[10,1]     0.28    0.00 0.05    0.18    0.24    0.27
#> spatialEffectsKnots[10,2]     0.91    0.00 0.04    0.83    0.88    0.91
#> spatialEffectsKnots[10,3]     0.79    0.00 0.06    0.68    0.75    0.79
#> spatialEffectsKnots[10,4]    -0.11    0.00 0.05   -0.21   -0.15   -0.11
#> spatialEffectsKnots[10,5]    -1.24    0.00 0.05   -1.34   -1.27   -1.24
#> spatialEffectsKnots[10,6]     1.38    0.00 0.05    1.29    1.34    1.38
#> spatialEffectsKnots[10,7]     0.44    0.00 0.05    0.33    0.40    0.43
#> spatialEffectsKnots[10,8]     0.29    0.00 0.05    0.18    0.25    0.29
#> spatialEffectsKnots[10,9]     0.04    0.00 0.05   -0.05    0.01    0.04
#> spatialEffectsKnots[10,10]    0.35    0.00 0.07    0.21    0.31    0.35
#> spatialEffectsKnots[11,1]    -0.12    0.00 0.05   -0.21   -0.15   -0.12
#> spatialEffectsKnots[11,2]     0.16    0.00 0.04    0.08    0.13    0.16
#> spatialEffectsKnots[11,3]    -0.14    0.00 0.06   -0.25   -0.17   -0.14
#> spatialEffectsKnots[11,4]    -0.46    0.00 0.05   -0.55   -0.49   -0.45
#> spatialEffectsKnots[11,5]     0.33    0.00 0.05    0.22    0.29    0.33
#> spatialEffectsKnots[11,6]     0.49    0.00 0.05    0.40    0.46    0.49
#> spatialEffectsKnots[11,7]    -0.15    0.00 0.06   -0.26   -0.19   -0.15
#> spatialEffectsKnots[11,8]    -0.39    0.00 0.05   -0.50   -0.43   -0.39
#> spatialEffectsKnots[11,9]     0.21    0.00 0.05    0.12    0.17    0.21
#> spatialEffectsKnots[11,10]    0.05    0.00 0.06   -0.08    0.01    0.05
#> spatialEffectsKnots[12,1]    -0.49    0.00 0.05   -0.58   -0.52   -0.49
#> spatialEffectsKnots[12,2]     0.02    0.00 0.04   -0.06   -0.01    0.02
#> spatialEffectsKnots[12,3]    -0.12    0.00 0.06   -0.23   -0.15   -0.12
#> spatialEffectsKnots[12,4]    -0.04    0.00 0.05   -0.14   -0.08   -0.04
#> spatialEffectsKnots[12,5]    -0.23    0.00 0.05   -0.34   -0.27   -0.23
#> spatialEffectsKnots[12,6]     0.39    0.00 0.04    0.31    0.36    0.39
#> spatialEffectsKnots[12,7]    -0.19    0.00 0.06   -0.29   -0.23   -0.18
#> spatialEffectsKnots[12,8]    -0.37    0.00 0.05   -0.48   -0.41   -0.37
#> spatialEffectsKnots[12,9]     0.02    0.00 0.05   -0.08   -0.01    0.02
#> spatialEffectsKnots[12,10]    0.09    0.00 0.06   -0.03    0.05    0.10
#> sigma[1]                      0.11    0.00 0.00    0.10    0.10    0.11
#> lp__                       2133.08    0.45 8.36 2115.35 2127.47 2133.20
#>                                75%   97.5% n_eff Rhat
#> df[1]                         4.88    9.17   800 1.00
#> gp_sigma                      0.27    0.32   800 1.00
#> gp_scale                      0.51    0.54   800 1.00
#> B[1]                          0.00    0.01   200 1.01
#> spatialEffectsKnots[1,1]     -0.13   -0.07   800 1.00
#> spatialEffectsKnots[1,2]     -0.07   -0.01   800 1.00
#> spatialEffectsKnots[1,3]      0.68    0.76   800 1.00
#> spatialEffectsKnots[1,4]     -0.22   -0.16   800 1.00
#> spatialEffectsKnots[1,5]      0.00    0.07   800 1.00
#> spatialEffectsKnots[1,6]     -0.02    0.04   800 1.00
#> spatialEffectsKnots[1,7]     -0.12   -0.06   800 1.00
#> spatialEffectsKnots[1,8]     -0.05    0.02   800 1.00
#> spatialEffectsKnots[1,9]      0.14    0.22   800 1.00
#> spatialEffectsKnots[1,10]    -0.11   -0.03   800 1.00
#> spatialEffectsKnots[2,1]     -0.02    0.04   800 1.00
#> spatialEffectsKnots[2,2]      0.29    0.34   800 1.00
#> spatialEffectsKnots[2,3]      0.04    0.11   800 1.00
#> spatialEffectsKnots[2,4]     -0.11   -0.05   800 1.00
#> spatialEffectsKnots[2,5]      0.02    0.07   800 1.00
#> spatialEffectsKnots[2,6]      0.17    0.22   800 1.00
#> spatialEffectsKnots[2,7]      0.05    0.12   800 1.00
#> spatialEffectsKnots[2,8]     -0.02    0.05   800 1.00
#> spatialEffectsKnots[2,9]     -0.13   -0.08   800 1.00
#> spatialEffectsKnots[2,10]     0.03    0.10   800 1.00
#> spatialEffectsKnots[3,1]     -0.01    0.05   800 1.00
#> spatialEffectsKnots[3,2]     -0.04    0.01   800 1.00
#> spatialEffectsKnots[3,3]      0.19    0.26   800 1.00
#> spatialEffectsKnots[3,4]     -0.15   -0.09   800 1.00
#> spatialEffectsKnots[3,5]      0.07    0.12   800 1.00
#> spatialEffectsKnots[3,6]     -0.12   -0.07   800 1.00
#> spatialEffectsKnots[3,7]      0.03    0.09   800 1.00
#> spatialEffectsKnots[3,8]      0.01    0.07   800 1.00
#> spatialEffectsKnots[3,9]      0.00    0.07   800 1.00
#> spatialEffectsKnots[3,10]     0.08    0.15   800 1.00
#> spatialEffectsKnots[4,1]     -0.53   -0.47   800 1.00
#> spatialEffectsKnots[4,2]      0.39    0.44   800 1.00
#> spatialEffectsKnots[4,3]     -0.43   -0.35   800 1.00
#> spatialEffectsKnots[4,4]     -0.05    0.01   800 1.00
#> spatialEffectsKnots[4,5]     -0.26   -0.19   800 1.00
#> spatialEffectsKnots[4,6]     -0.20   -0.14   800 1.00
#> spatialEffectsKnots[4,7]      0.68    0.75   800 1.00
#> spatialEffectsKnots[4,8]      0.00    0.08   800 1.00
#> spatialEffectsKnots[4,9]     -0.32   -0.25   800 1.00
#> spatialEffectsKnots[4,10]    -0.40   -0.32   800 1.00
#> spatialEffectsKnots[5,1]      0.12    0.18   800 1.00
#> spatialEffectsKnots[5,2]      0.01    0.06   800 1.00
#> spatialEffectsKnots[5,3]      0.00    0.06   800 1.00
#> spatialEffectsKnots[5,4]     -0.05    0.02   800 1.00
#> spatialEffectsKnots[5,5]     -0.23   -0.17   800 1.00
#> spatialEffectsKnots[5,6]     -0.19   -0.14   800 1.00
#> spatialEffectsKnots[5,7]      0.11    0.18   800 1.00
#> spatialEffectsKnots[5,8]     -0.01    0.06   800 1.00
#> spatialEffectsKnots[5,9]     -0.20   -0.14   800 1.00
#> spatialEffectsKnots[5,10]     0.31    0.39   800 1.00
#> spatialEffectsKnots[6,1]      0.15    0.21   800 1.00
#> spatialEffectsKnots[6,2]      0.00    0.05   800 1.00
#> spatialEffectsKnots[6,3]      0.29    0.37   800 1.00
#> spatialEffectsKnots[6,4]      0.24    0.31   800 1.00
#> spatialEffectsKnots[6,5]     -0.14   -0.05   800 1.00
#> spatialEffectsKnots[6,6]      0.57    0.63   800 1.00
#> spatialEffectsKnots[6,7]      0.10    0.18   800 1.00
#> spatialEffectsKnots[6,8]     -0.33   -0.27   800 1.00
#> spatialEffectsKnots[6,9]      0.01    0.07   800 1.00
#> spatialEffectsKnots[6,10]     0.14    0.22   800 1.00
#> spatialEffectsKnots[7,1]      1.17    1.23   800 1.00
#> spatialEffectsKnots[7,2]      0.07    0.12   800 1.00
#> spatialEffectsKnots[7,3]      0.27    0.35   800 1.00
#> spatialEffectsKnots[7,4]     -0.07   -0.01   800 1.00
#> spatialEffectsKnots[7,5]     -0.12   -0.05   800 1.00
#> spatialEffectsKnots[7,6]     -0.02    0.04   800 1.00
#> spatialEffectsKnots[7,7]      0.48    0.55   618 1.00
#> spatialEffectsKnots[7,8]      1.06    1.12   800 1.00
#> spatialEffectsKnots[7,9]      0.51    0.57   800 1.00
#> spatialEffectsKnots[7,10]     0.70    0.79   800 1.00
#> spatialEffectsKnots[8,1]     -0.19   -0.13   800 1.00
#> spatialEffectsKnots[8,2]      0.25    0.30   800 1.00
#> spatialEffectsKnots[8,3]      0.19    0.27   800 1.00
#> spatialEffectsKnots[8,4]     -0.24   -0.17   800 1.00
#> spatialEffectsKnots[8,5]      0.19    0.26   800 1.00
#> spatialEffectsKnots[8,6]      0.05    0.10   800 1.00
#> spatialEffectsKnots[8,7]      0.31    0.38   800 1.00
#> spatialEffectsKnots[8,8]     -0.09   -0.02   800 1.00
#> spatialEffectsKnots[8,9]     -0.11   -0.06   800 1.00
#> spatialEffectsKnots[8,10]    -0.12   -0.04   800 1.00
#> spatialEffectsKnots[9,1]     -0.01    0.06   800 1.00
#> spatialEffectsKnots[9,2]      0.20    0.25   800 1.00
#> spatialEffectsKnots[9,3]     -0.21   -0.14   800 1.00
#> spatialEffectsKnots[9,4]      0.44    0.50   800 1.00
#> spatialEffectsKnots[9,5]     -0.47   -0.41   800 1.00
#> spatialEffectsKnots[9,6]     -0.37   -0.31   800 1.00
#> spatialEffectsKnots[9,7]      0.61    0.68   800 1.00
#> spatialEffectsKnots[9,8]     -0.36   -0.28   800 1.00
#> spatialEffectsKnots[9,9]      0.21    0.27   800 1.00
#> spatialEffectsKnots[9,10]    -0.14   -0.05   800 1.00
#> spatialEffectsKnots[10,1]     0.31    0.38   800 1.00
#> spatialEffectsKnots[10,2]     0.94    0.99   800 1.00
#> spatialEffectsKnots[10,3]     0.83    0.91   800 1.00
#> spatialEffectsKnots[10,4]    -0.08   -0.02   800 1.00
#> spatialEffectsKnots[10,5]    -1.20   -1.14   800 1.00
#> spatialEffectsKnots[10,6]     1.41    1.47   800 1.00
#> spatialEffectsKnots[10,7]     0.47    0.55   800 1.00
#> spatialEffectsKnots[10,8]     0.32    0.39   800 1.00
#> spatialEffectsKnots[10,9]     0.08    0.15   800 1.00
#> spatialEffectsKnots[10,10]    0.40    0.48   800 1.00
#> spatialEffectsKnots[11,1]    -0.08   -0.03   800 1.00
#> spatialEffectsKnots[11,2]     0.18    0.24   800 1.00
#> spatialEffectsKnots[11,3]    -0.10   -0.02   800 1.00
#> spatialEffectsKnots[11,4]    -0.42   -0.36   800 1.00
#> spatialEffectsKnots[11,5]     0.36    0.43   800 1.00
#> spatialEffectsKnots[11,6]     0.52    0.57   800 1.00
#> spatialEffectsKnots[11,7]    -0.11   -0.04   800 1.00
#> spatialEffectsKnots[11,8]    -0.36   -0.29   800 1.00
#> spatialEffectsKnots[11,9]     0.24    0.30   800 1.00
#> spatialEffectsKnots[11,10]    0.10    0.17   800 1.00
#> spatialEffectsKnots[12,1]    -0.46   -0.40   800 1.01
#> spatialEffectsKnots[12,2]     0.05    0.10   800 1.00
#> spatialEffectsKnots[12,3]    -0.08    0.00   800 1.00
#> spatialEffectsKnots[12,4]    -0.01    0.05   800 1.00
#> spatialEffectsKnots[12,5]    -0.20   -0.14   800 1.00
#> spatialEffectsKnots[12,6]     0.42    0.48   800 1.00
#> spatialEffectsKnots[12,7]    -0.15   -0.07   800 1.00
#> spatialEffectsKnots[12,8]    -0.34   -0.27   800 1.00
#> spatialEffectsKnots[12,9]     0.05    0.11   800 1.00
#> spatialEffectsKnots[12,10]    0.14    0.21   800 1.00
#> sigma[1]                      0.11    0.11   800 1.00
#> lp__                       2138.53 2148.67   347 1.02
#> 
#> Samples were drawn using NUTS(diag_e) at Sat Nov 19 19:01:26 2016.
#> For each parameter, n_eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor on split chains (at 
#> convergence, Rhat=1).
#> 
#> $knots
#>             lon       lat
#>  [1,] 1.8621760 6.4279549
#>  [2,] 3.8003518 4.4628435
#>  [3,] 6.6200508 2.7775593
#>  [4,] 8.7626921 9.2730209
#>  [5,] 8.6969085 2.2865814
#>  [6,] 7.7732070 6.0530345
#>  [7,] 6.9273156 8.6454495
#>  [8,] 3.4668349 8.5613166
#>  [9,] 5.1863426 0.7527575
#> [10,] 0.9946616 1.8086636
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
