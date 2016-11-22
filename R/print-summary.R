#' @export
print.rrfield <- function(object, ...) {
  print(object$model, ...)
}

#' @export
tidy <- function(object, ...){
  UseMethod("tidy")
}

#' @export
tidy.rrfield <- function(object, ...) {
  broom::tidyMCMC(object$model, ...)
}

