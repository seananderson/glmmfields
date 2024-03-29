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
gp_theta <- 1.2
n_draws <- 15
nknots <- 7
n_data_points <- 50

# ------------------------------------------------------
# with repeat stations

test_that("mvt-norm model fits with repeat stations (plus other main functions)", {
  skip_on_cran()
  set.seed(SEED)

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points
  )
  # s$plot

  suppressWarnings({
  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df
  )
  })

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
  expect_equal(as.numeric(b[b$term == "sigma[1]", "estimate", drop = TRUE]), sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
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
  gp_theta <- 1.2
  n_draws <- 4
  nknots <- 9

  set.seed(SEED)
  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points,
    covariance = "exponential"
  )
  # print(s$plot)

  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    covariance = "exponential"
  )
  })

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "sigma[1]", "estimate", drop = TRUE]), sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
})

# ------------------------------------------------------
# a Gaussian observation model matern covariance function

test_that("mvn-norm model fits with an matern covariance function", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_theta <- 1.2
  n_draws <- 4
  nknots <- 9
  matern_kappa <- 1.5

  set.seed(SEED)
  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points,
    covariance = "matern", matern_kappa = matern_kappa
  )
  # print(s$plot)


  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    covariance = "matern", matern_kappa = matern_kappa
  )
  })

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "sigma[1]", "estimate", drop = TRUE]), sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
})

test_that("predictions work with one time slice", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  s <- sim_glmmfields(
    df = df, n_draws = 1, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points
  )

  suppressWarnings({
  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df
  )
  })

  p <- predict(m)
})

# -------------------------------------------------
# make sure large degrees of freedom values
# return values very close to the true MVN distribution

test_that("true MVN model closely resembles MVT model with a large fixed df", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  gp_theta <- 1.2
  n_draws <- 4
  nknots <- 9

  set.seed(SEED)
  s <- sim_glmmfields(
    n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    df = 800, n_data_points = n_data_points
  )

  suppressWarnings({
  m_mvt <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = 800
  )
  })

  suppressWarnings({
  m_mvn <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = 1e9
  ) # internally switched to true MVN
  })

  b_mvt <- tidy(m_mvt, estimate.method = "median")
  b_mvn <- tidy(m_mvn, estimate.method = "median")
  expect_equal(as.numeric(b_mvn$estimate), as.numeric(b_mvt$estimate), tol = 0.02)
})

# -------------------------------------------------------------------

test_that("A basic model fits with missing time elements", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_theta <- 1.2
  n_draws <- 10
  nknots <- 9

  set.seed(SEED * 2)
  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points,
    covariance = "squared-exponential"
  )
  # print(s$plot)

  s$dat <- s$dat[s$dat$time != 4, , drop = FALSE]

  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    covariance = "squared-exponential"
  )
  })

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "sigma[1]", "estimate", drop = TRUE]), sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL * 1.5)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
})


test_that("A offset in formula works", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_theta <- 1.2
  n_draws <- 2
  nknots <- 9

  set.seed(SEED * 2)
  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points,
    covariance = "squared-exponential"
  )
  # print(s$plot)
  s$dat$offset <- rnorm(nrow(s$dat), mean = 0, sd = 0.1)

  suppressWarnings({
  m <- glmmfields(y ~ 1,
                  data = s$dat, time = "time",
                  lat = "lat", lon = "lon", nknots = nknots,
                  iter = ITER, chains = CHAINS, seed = SEED,
                  estimate_df = FALSE, fixed_df_value = df,
                  covariance = "squared-exponential"
  )
  })
  suppressWarnings({
  m_offset <- glmmfields(y ~ 1, offset = s$dat$offset,
                  data = s$dat, time = "time",
                  lat = "lat", lon = "lon", nknots = nknots,
                  iter = ITER, chains = CHAINS, seed = SEED,
                  estimate_df = FALSE, fixed_df_value = df,
                  covariance = "squared-exponential"
  )
  })
  b <- tidy(m, estimate.method = "median")
  b_offset <- tidy(m_offset, estimate.method = "median")

  p1 <- predict(m_offset)
  p2 <- predict(m_offset, newdata = s$dat, offset = s$dat$offset)
  expect_identical(p1, p2)
  # expect_error(p3 <- predict(m_offset, newdata = s$dat), regexp = "offset")
})
