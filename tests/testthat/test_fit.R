skip_on_cran()

test_that("mvt-norm ar1 model fits", {
  library(rrfields)
  set.seed(1234)
  n_draws <- 5
  s <- sim_mvt_rf(df = 5, n_draws = n_draws, n_knots = 8)
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

  d <- format_data(out, y = "y", time = "time_slice", nKnots = 8)
  m <- rrfield(data = d, iter = 300, chains = 1, control = list(adapt_delta = 0.8))
  expect_is(m, "stanfit")
})
