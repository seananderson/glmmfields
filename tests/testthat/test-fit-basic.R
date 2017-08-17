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

gp_eta <- 0.2
sigma <- 0.1
df <- 4
gp_scale <- 1.2
n_draws <- 15
nknots <- 7
n_data_points <- 50

# ------------------------------------------------------
# with repeat stations

test_that("mvt-norm model fits with repeat stations (plus other main functions)", {
  skip_on_cran()
  set.seed(SEED)

  s <- sim_glmmfields(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_eta = gp_eta, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points)
  # s$plot

  m <- glmmfields(y ~ 0, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  expect_output(print(m), "Inference for Stan model")

  p <- predict(m)
  pp <- predict(m, type = "response", interval = "prediction")
  # plot(s$dat$y, p$estimate)
  # segments(s$dat$y, pp$conf_low, s$dat$y, pp$conf_high, lwd = 0.5, col = "#00000020")
  # segments(s$dat$y, p$conf_low, s$dat$y, p$conf_high, lwd = 1, col = "#00000060")
  # abline(a = 0, b = 1)

  expect_equal(mean((p$estimate - s$dat$y)^2), 0, tol = 0.01)

  plot(m)
  plot(m, type = "spatial-residual")
  plot(m, type = "residual-vs-fitted")

  coverage <- mean(s$dat$y > pp$conf_low & s$dat$y < pp$conf_high)
  expect_equal(coverage, 0.95, tol = 0.025)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_eta", "estimate"], gp_eta, tol = gp_eta * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# a Gaussian observation model exponential covariance function

test_that("mvt-norm model fits with an exponential covariance function", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_eta <- 0.2
  sigma <- 0.1
  df <- 10
  gp_scale <- 1.2
  n_draws <- 4
  nknots <- 9

  set.seed(SEED)
  s <- sim_glmmfields(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_eta = gp_eta, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points,
    covariance = "exponential")
  # print(s$plot)

  m <- glmmfields(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    covariance = "exponential")

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_eta", "estimate"], gp_eta, tol = gp_eta * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

test_that("predictions work with one time slice", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  s <- sim_glmmfields(df = df, n_draws = 1, gp_scale = gp_scale,
    gp_eta = gp_eta, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points)

  m <- glmmfields(y ~ 0, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  p <- predict(m)

})

# -------------------------------------------------
# make sure large degrees of freedom values
# return values very close to the true MVN distribution

test_that("true MVN model closely resembles MVT model with a large fixed df", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_eta <- 0.2
  sigma <- 0.1
  gp_scale <- 1.2
  n_draws <- 4
  nknots <- 9

  set.seed(SEED)
  s <- sim_glmmfields(n_draws = n_draws, gp_scale = gp_scale,
    gp_eta = gp_eta, sd_obs = sigma, n_knots = nknots,
    df = 800, n_data_points = n_data_points)

  m_mvt <- glmmfields(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = 800)

  m_mvn <- glmmfields(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = 1e9) # internally switched to true MVN

  b_mvt <- tidy(m_mvt, estimate.method = "median")
  b_mvn <- tidy(m_mvn, estimate.method = "median")
  expect_equal(b_mvn$estimate, b_mvt$estimate, tol = 0.02)
})
