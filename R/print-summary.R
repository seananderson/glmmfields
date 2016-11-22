#' @export
#' @import methods
print.rrfield <- function(x, ...) {
  print(x$model, ...)
}

tidy <- function(x, ...){
  UseMethod("tidy")
}

#' @export
tidy.rrfield <- function(x, ...) {
  broom::tidyMCMC(x$model, ...)
}

