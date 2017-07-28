#' Plot predictions from an glmmfields model
#'
#' @param x An object returned by \code{\link{glmmfields}}
#' @param type Type of plot
#' @param link Logical: should the plots be made on the link scale
#'  or on the natural scale?
#' @param ... Other arguments passed to \code{\link{predict.glmmfields}}
#'
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point
#'   scale_color_gradient2 geom_smooth geom_hline facet_wrap

#' @export

plot.glmmfields <- function(x,
  type = c("prediction", "spatial-residual", "residual-vs-fitted"),
  link = TRUE, ...)  {

  p <- predict(x, type = ifelse(link, "link", "response"), ...)
  d <- data.frame(x$data, p)
  y <- x$y
  if (link) y <- do.call(x$family$link, list(y))
  d$residual <- y - p$estimate

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

  g
}
