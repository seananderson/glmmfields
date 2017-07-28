#' @export
#' @import methods
print.glmmfields <- function(x, pars = c("spatialEffectsKnots", "log_lik"),
  include = FALSE, ...) {
  print(x$model, pars = pars, include = include, ...)
}

#' Tidy model output
#'
#' @param x Output from \code{\link{glmmfields}}
#' @param ... Other arguments
#' @export
tidy <- function(x, ...){
  UseMethod("tidy")
}

#' Tidy model output
#'
#' @param x Output from \code{\link{glmmfields}}
#' @param ... Other arguments
#' @export
tidy.glmmfields <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}

