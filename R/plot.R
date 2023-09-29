#' Plot predictions from an glmmfields model
#'
#' @param x An object returned by \code{\link{glmmfields}}
#' @param type Type of plot
#' @param link Logical: should the plots be made on the link scale
#'  or on the natural scale?
#' @param ... Other arguments passed to \code{\link{predict.glmmfields}}
#'
#' @importFrom ggplot2 aes ggplot facet_wrap geom_point
#'   scale_color_gradient2 geom_smooth geom_hline facet_wrap
#' @export
#' @examples
#' \donttest{
#' # Spatiotemporal example:
#' set.seed(1)
#' s <- sim_glmmfields(n_draws = 12, n_knots = 12, gp_theta = 1.5,
#' gp_sigma = 0.2, sd_obs = 0.1)
#' # options(mc.cores = parallel::detectCores()) # for parallel processing
#' m <- glmmfields(y ~ 0, time = "time",
#'  lat = "lat", lon = "lon", data = s$dat,
#'  nknots = 12, iter = 600, chains = 1)
#' x <- plot(m, type = "prediction")
#' x
#' x + ggplot2::scale_color_gradient2()
#' plot(m, type = "spatial-residual")
#' plot(m, type = "residual-vs-fitted")
#' }

plot.glmmfields <- function(x,
                            type = c("prediction", "spatial-residual", "residual-vs-fitted"),
                            link = TRUE, ...) {
  type <- match.arg(type)

  p <- predict(x, type = ifelse(link, "link", "response"), ...)
  d <- data.frame(x$data, p)
  y <- x$y
  if (link) y <- do.call(x$family$link, list(y))
  d$residual <- y - p$estimate

  g <- NULL

  if (type == "prediction") {
    g <- ggplot(d, aes(d[[x$lon]], d[[x$lat]], colour = "estimate")) +
      geom_point(size = 2) +
      facet_wrap(x$time)
  }

  if (type == "spatial-residual") {
    g <- ggplot(d, aes(d[[x$lon]], d[[x$lat]], colour = "residual")) +
      geom_point(size = 2) +
      scale_color_gradient2() +
      facet_wrap(x$time)
  }

  if (type == "residual-vs-fitted") {
    g <- ggplot(d, aes(d[["estimate"]], d[["residual"]])) +
      geom_point() +
      facet_wrap(x$time) +
      geom_hline(yintercept = 0, lty = 2) +
      geom_smooth(method = "loess", se = FALSE, colour = "red")
  }

  g
}
