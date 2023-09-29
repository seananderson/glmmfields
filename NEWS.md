# glmmfields 0.1.8

* Update for rstan 2.26 (#16)

# glmmfields 0.1.7

* Rebuild with newest rstantools to avoid CRAN NOTE about C++ version

# glmmfields 0.1.6

* Add back RcppParallel as an import

# glmmfields 0.1.5

* Fully delegate installation to rstantools (#15)

* Other minor fixes to pass R CMD check on R devel

* Add offset functionality

# glmmfields 0.1.4

* Make compatible with updated broom.mixed

# glmmfields 0.1.3

* Make compatible with C++14.

# glmmfields 0.1.2

* Make compatible with R 3.6.0 staged installation and latest rstantools.

# glmmfields 0.1.1

* Changed how Stan finds directories / files

# glmmfields 0.1.0.9000

* Add support for random walk year effects with covariates. There are a few
  specific cases where covariates or the random year effect are not estimable.
  Examples are:

    1. estimating an intercept and AR year effects (intercept confounded with
       1st year effect)
    2. estimating AR year effects and estimating 'phi' — the AR parameter in
       spatial field (also confounded)

* Import S3 methods from rstantools instead of rstanarm (#5)

* Adjust calculation of year index values to better allow for missing time slices

# glmmfields 0.1.0

* Initial submission to CRAN.
