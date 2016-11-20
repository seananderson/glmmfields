library(rrfields)
library(ggplot2)
library(testthat)
if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 400
CHAINS <- 2
SEED <- 123

gp_sigma <- 0.2
sigma <- 0.1
df <- 3
gp_scale <- 1.2
n_draws <- 5
nknots <- 15
TOL <- 0.2 # %
TOL_df <- .25 # %

skip_on_cran()
test_that("mvt-norm model fits", {
  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, obs_error = "normal", seed = SEED)

  b <- broom::tidyMCMC(m$model, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
})
