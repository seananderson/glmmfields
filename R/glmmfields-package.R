#' The 'glmmfields' package.
#'
#' @description Implements Bayesian spatial and spatiotemporal models that
#'   optionally allow for extreme spatial deviations through time. 'glmmfields'
#'   uses a predictive process approach with random fields implemented through a
#'   multivariate-t distribution instead of the usual multivariate normal.
#'   Sampling is conducted with 'Stan'.
#'
#' @docType package
#' @name glmmfields-package
#' @useDynLib glmmfields, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @import rstantools
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2018). RStan: the R interface to Stan. R package
#' version 2.18.2. http://mc-stan.org
#'
NULL
