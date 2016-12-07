if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 999
TOL <- 0.2 # %
TOL_df <- .25 # %

# ------------------------------------------------------
# a basic fit

gp_sigma <- 0.2
sigma <- 0.1
df <- 4
gp_scale <- 1.2
n_draws <- 10
nknots <- 10

# ------------------------------------------------------
# with repeat stations

test_that("mvt-norm model fits with repeat stations", {
  skip_on_cran()
  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)

  m <- rrfield(y ~ 1, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  b <- broom::tidyMCMC(m$model, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# without repeat stations

test_that("mvt-norm model fits without station argument", {
  skip_on_cran()
  skip_on_travis()

  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)
  # print(s$plot)

  m1 <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  print(m)
  p <- predict(m, newdata = s$dat, mcmc_draws = 200)

  b1 <- broom::tidyMCMC(m1$model, estimate.method = "median")
  expect_equal(b1[b1$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b1[b1$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b1[b1$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)

  expect_equal(b$estimate, b1$estimate, tol = 0.01) # w or w/o station arg
})

# ------------------------------------------------------
# a negative binomial model

test_that("mvt-nb2 model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  sigma <- 8
  df <- 5
  b0 <- 7
  n_draws <- 15
  gp_scale <- 1.6

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "nb2", B = b0)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, obs_error = "nb2",
    estimate_df = FALSE, fixed_df_value = df,
    control = list(adapt_delta = 0.9), seed = SEED)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "nb2_phi[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[b$term == "B[1]", "estimate"], b0, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# a Gamma observation model

test_that("mvt-gamma model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  sigma <- 0.3
  df <- 10
  b0 <- 2
  n_draws <- 15
  gp_scale <- 1.6

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "gamma", B = b0)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, obs_error = "gamma",
    estimate_df = FALSE, fixed_df_value = df, seed = SEED)
  print(m)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "CV[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[b$term == "B[1]", "estimate"], b0, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# a Gaussian observation model with factor-level predictors for years

test_that("mvt-norm estimates betas", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_scale <- 1.2
  n_draws <- 15
  nknots <- 9
  set.seed(SEED)
  B <- rnorm(n_draws, 0, 1)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, B = B,
    X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  m <- rrfield(y ~ as.factor(time) - 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots, station = "station_id",
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    prior_beta = student_t(50, 0, 2))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[grep("B\\[*", b$term), "estimate"], B, tol = 0.05)
})

# ------------------------------------------------------
# a Gaussian observation model exponential covariance function

test_that("mvt-norm model fits with an exponential covariance function", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 4
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
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
})

# ------------------------------------------------------
# a Gaussian observation model with random walk year effects

test_that("mvt-norm estimates random walk year effects", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_scale <- 1.2
  n_draws <- 15
  nknots <- 11
  year_sigma <- 0.9
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0.5
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, B = B,
    X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  m <- rrfield(y ~ 1, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df, year_re = TRUE,
    # TODO fake tight prior until I disable the intercept:
    prior_intercept = student_t(999, 0, 0.3))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[grep("yearEffects\\[*", b$term), "estimate"], B, tol = B * TOL)
  expect_equal(b[grep("year_sigma", b$term), "estimate"], year_sigma, tol = B * TOL)
})
