library(rrfields)
# library(rstan)
if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 400
CHAINS <- 2
SEED <- 123

skip_on_cran()
test_that("mvt-norm model fits", {
  set.seed(SEED)
  n_draws <- 5
  s <- sim_rrfield(df = 3, n_draws = n_draws)
  library(ggplot2)
  g <- ggplot(s$dat, aes(x = lon, y = lat, z = y, colour = y)) +
    facet_wrap(~time, ncol = n_draws) +
    geom_point(size = 3) +
    viridis::scale_color_viridis() +
    theme_light()
  # print(g)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = 15,
    iter = ITER, chains = CHAINS)

  # expect_equal(fixef(fit), fixef(ans), tol = FIXEF_tol)
})
