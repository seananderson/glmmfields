if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.2 # %
TOL_df <- .25 # %

gp_sigma <- 0.2
sigma <- 0.1
df <- 4
gp_scale <- 1.2
n_draws <- 15
nknots <- 10

test_that("predict.rrfield works", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  p <- predict(m)
  p_newdata <- predict(m, newdata = s$dat)
  p_newdata2 <- predict(m, newdata = m$data)

  plot(s$dat$y, p$estimate)
  plot(s$dat$y, p_newdata$estimate)
  plot(s$dat$y, p_newdata2$estimate)

  expect_identical(p, p_newdata)

  expect_gte(cor(s$dat$y, p$estimate), 0.75)
  expect_gte(cor(s$dat$y, p_newdata$estimate), 0.75)
})
