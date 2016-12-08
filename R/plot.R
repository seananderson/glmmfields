#' Plot predictions from an rrfield model
#'
#' @param object An object returned by \code{\link{rrfield}}
#' @param ... Other arguments passed to \code{\link{predict.rrfield}}
#'
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point scale_color_gradient2

#' @export

plot.rrfield <- function(object, ...)  {
  d <- data.frame(object$data, predict(object, ...))
  ggplot(d, aes_string(object$lon, object$lat, colour = "estimate")) +
    geom_point(size = 2) +
    scale_color_gradient2() +
    facet_wrap(object$time)
}
