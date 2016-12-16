if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.25 # %
TOL_df <- .25 # %

# ------------------------------------------------------
# a Gaussian observation model with random walk year effects

test_that("mvt-norm estimates random walk year effects", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED*2)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 10
  gp_scale <- 1.8
  n_draws <- 15
  nknots <- 10
  year_sigma <- 0.5
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, B = B,
    X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df, year_re = TRUE,
    prior_intercept = student_t(999, 0, 5), control = list(adapt_delta = 0.9),
    prior_rw_sigma = half_t(1e6, 0, 1))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[grep("yearEffects\\[*", b$term), "estimate"], B, tol = 0.1)
  expect_equal(b[grep("year_sigma", b$term), "estimate"], year_sigma, tol = 0.1)
})

# ---------------
# AR process

test_that("mvt-norm estimates ar process", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED*2)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 8
  gp_scale <- 1.8
  n_draws <- 20
  nknots <- 10
  ar <- 0.5

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar)
  print(s$plot)
  library(ggplot2); ggplot(s$dat, aes(time, y)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = TRUE,
    control = list(adapt_delta = 0.9),
    estimate_ar = TRUE)
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
})

# -------------------
# AR + year random effects

test_that("mvt-norm estimates ar process *with* year random walk effects", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED*2)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 8
  gp_scale <- 1.8
  n_draws <- 20
  nknots <- 10
  ar <- 0.5
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  year_sigma <- 0.3
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
    B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  print(s$plot)
  library(ggplot2); ggplot(s$dat, aes(time, y)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER * 2, chains = CHAINS, seed = SEED,
    estimate_df = TRUE, year_re = TRUE,
    control = list(adapt_delta = 0.9),
    estimate_ar = TRUE)
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
  expect_equal(b[b$term == "year_sigma[1]", "estimate"], year_sigma, tol = year_sigma * TOL_df)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
})

# --------------------------
# global int + AR RF

test_that("mvt-norm estimates global int + AR RF", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED*2)

  gp_sigma <- 0.2
  sigma <- 0.2
  df <- 6
  gp_scale <- 1.8
  n_draws <- 20
  nknots <- 10
  ar <- 0.75
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 6
  year_sigma <- 0.0001
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
    B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  print(s$plot)
  library(ggplot2); ggplot(s$dat, aes(time, y)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- rrfield(y ~ 1, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = TRUE, year_re = FALSE,
    control = list(adapt_delta = 0.95),
    estimate_ar = TRUE, prior_intercept = student_t(99, 0, 30))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
})
