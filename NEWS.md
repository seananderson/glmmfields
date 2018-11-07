# glmmfields 0.1.1

* Changed how Stan files are found, per Stan developers

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
