if (interactive()) options(mc.cores = parallel::detectCores())

# expect_parameter <- function(x, stan_term, true_value, tolerance) {
#   expect_equal(x[x[,term] == stan_term, "estimate"], true_value, tol = tolerance)
# }

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.2 # %
TOL_df <- .25 # %

# ------------------------------------------------------
# a basic fit

gp_sigma <- 0.2
sigma <- 0.1
df <- 4
gp_scale <- 1.2
n_draws <- 15
nknots <- 10

# ------------------------------------------------------
# with repeat stations

test_that("mvt-norm model fits with repeat stations (plus other main functions)", {
  skip_on_cran()
  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)
  # s$plot

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  expect_output(print(m), "Inference for Stan model")

  p <- predict(m)
  pp <- predict(m, interval = "prediction")
  plot(s$dat$y, p$estimate)
  segments(s$dat$y, pp$conf_low, s$dat$y, pp$conf_high, lwd = 0.5, col = "#00000020")
  segments(s$dat$y, p$conf_low, s$dat$y, p$conf_high, lwd = 1, col = "#00000060")
  abline(a = 0, b = 1)

  expect_equal(mean((p$estimate - s$dat$y)^2), 0, tol = 0.01)

  plot(m)
  plot(m, type = "spatial-residual")
  plot(m, type = "residual-vs-fitted")

  coverage <- mean(s$dat$y > pp$conf_low & s$dat$y < pp$conf_high)
  expect_equal(coverage, 0.95, tol = 0.025)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# without repeat stations

test_that("mvt-norm model fits without station argument", {
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
  b <- broom::tidyMCMC(m$model, estimate.method = "median")

  m_nostation <- rrfield(y ~ 0, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)
  b_nostation <- broom::tidyMCMC(m1$model, estimate.method = "median")

  expect_equal(b$estimate, b_nostation$estimate, tol = 0.02) # w or w/o station arg
})

# ------------------------------------------------------
# a Gaussian observation model exponential covariance function

test_that("mvt-norm model fits with an exponential covariance function", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_scale <- 1.2
  n_draws <- 4
  nknots <- 9

  set.seed(SEED)
  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    covariance = "exponential")
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots, station = "station_id",
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    covariance = "exponential")

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# with repeat stations but missing in some years

test_that("mvt-norm model fits with repeat stations but missing in some years", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 20
  gp_scale <- 1.2
  n_draws <- 15
  nknots <- 10

  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)

  s$dat <- s$dat[-1, ] # remove some

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  b <- broom::tidyMCMC(m$model, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})
