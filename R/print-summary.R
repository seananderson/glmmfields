#' @export
#' @import methods
print.rrfield <- function(x, pars = "spatialEffectsKnots", include = FALSE, ...) {
  print(x$model, pars = pars, include = include, ...)
}

#' Tidy model output
#'
#' @param x Output from \code{\link{rrfield}}
#' @param ... Other arguments
#' @export
tidy <- function(x, ...){
  UseMethod("tidy")
}

#' Tidy model output
#'
#' @param x Output from \code{\link{rrfield}}
#' @param ... Other arguments
#' @export
tidy.rrfield <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}

