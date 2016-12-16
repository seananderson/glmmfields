#' Plot predictions from an rrfield model
#'
#' @param x An object returned by \code{\link{rrfield}}
#' @param type Type of plot
#' @param ... Other arguments passed to \code{\link{predict.rrfield}}
#'
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point
#'   scale_color_gradient2 geom_smooth geom_hline facet_wrap
#'   geom_contour

#' @export

plot.rrfield <- function(x, type = c("prediction", "contour", "spatial-residual", "residual-vs-fitted"), ...)  {

  p <- predict(x, ...)
  d <- data.frame(x$data, p)
  d$residual <- x$y - p$estimate

  g <- NULL

  if (type[[1]] == "prediction") {
    g <- ggplot(d, aes_string(x$lon, x$lat, colour = "estimate")) +
      geom_point(size = 2) +
      facet_wrap(x$time)
  }

  if (type[[1]] == "spatial-residual") {
    g <- ggplot(d, aes_string(x$lon, x$lat, colour = "residual")) +
      geom_point(size = 2) +
      scale_color_gradient2() +
      facet_wrap(x$time)
  }

  if (type[[1]] == "residual-vs-fitted") {
    g <- ggplot(d, aes_string("estimate", "residual")) +
      geom_point() +
      facet_wrap(x$time) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_smooth(method = "loess", se = FALSE, colour = "red")
  }

  if (type[[1]] == "contour") {
    g <- ggplot(d, aes_string(x$lon, x$lat, z = "estimate")) +
      geom_contour() +
      facet_wrap(x$time)
  }

  g
}
