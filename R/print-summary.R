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
#' @rdname tidy
tidy <- function(x, ...) {
  UseMethod("tidy")
}

#' @export
#' @rdname tidy
tidy.glmmfields <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}
