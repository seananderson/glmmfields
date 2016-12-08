#' Plot predictions from an rrfield model
#'
#' @param x An object returned by \code{\link{rrfield}}
#' @param ... Other arguments passed to \code{\link{predict.rrfield}}
#'
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point scale_color_gradient2

#' @export

plot.rrfield <- function(x, ...)  {
  d <- data.frame(x$data, predict(x, ...))
  ggplot(d, aes_string(x$lon, x$lat, colour = "estimate")) +
    geom_point(size = 2) +
    scale_color_gradient2() +
    facet_wrap(x$time)
}
