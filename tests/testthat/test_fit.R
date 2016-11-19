skip_on_cran()

test_that("mvt-norm ar1 model fits", {
  library(rrfields)
  set.seed(99)
  n_draws <- 10
  s <- sim_mvt_rf(df = 7, n_draws = n_draws, n_knots = 12)
  out <- reshape2::melt(s$proj)
  names(out) <- c("time_slice", "pt", "y")
  out <- dplyr::arrange(out, time_slice, pt)
  out$lon <- rep(s$g$lon, n_draws)
  out$lat <- rep(s$g$lat, n_draws)
  # ggplot(out, aes(x = lon, y = lat, z = y, colour = y)) +
  #   facet_wrap(~time_slice, ncol = n_draws) +
  #   geom_point(size = 3) +
  #   viridis::scale_color_viridis() +
  #   theme_light()

  m <- rrfield(y ~ 1, data = out, time = "time_slice",
    lat = "lat", lon = "lon", nKnots = 12,
    iter = 500, chains = 1, control = list(adapt_delta = 0.9))
  expect_is(m, "stanfit")
})
