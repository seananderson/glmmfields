#' @export
#' @import methods
print.rrfield <- function(x, ...) {
  print(x$model, ...)
}

#' Tidy model output
#'
#' @param x Output from \code{\link{rrfield}}
#' @export
tidy <- function(x, ...){
  UseMethod("tidy")
}

#' Tidy model output
#'
#' @param x Output from \code{\link{rrfield}}
#' @export
tidy.rrfield <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}

