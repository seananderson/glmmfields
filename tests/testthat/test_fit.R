skip_on_cran()

test_that("mvt-norm model fits", {
  set.seed(123)
  n_draws <- 5
  s <- sim_rrfield(df = 3, n_draws = n_draws)
  library(ggplot2)
  g <- ggplot(s$dat, aes(x = lon, y = lat, z = y, colour = y)) +
    facet_wrap(~time, ncol = n_draws) +
    geom_point(size = 3) +
    viridis::scale_color_viridis() +
    theme_light()
  print(g)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = 15, estimate_df = FALSE,
    iter = 300, chains = 1)
  expect_is(m, "stanfit")
})
