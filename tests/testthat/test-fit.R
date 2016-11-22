if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 400
CHAINS <- 2
SEED <- 123

gp_sigma <- 0.2
sigma <- 0.1
df <- 4
gp_scale <- 1.2
n_draws <- 4
nknots <- 9
TOL <- 0.2 # %
TOL_df <- .25 # %

test_that("mvt-norm model fits", {
  skip_on_cran()

  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df)

  print(m)
  p <- predict(m, newdata = s$dat, mcmc_draws = 200)

  # med <- apply(p, 1, median)
  # l <- apply(p, 1, quantile, probs = 0.025)
  # u <- apply(p, 1, quantile, probs = 0.975)
  # plot(s$dat$y, med);abline(a = 0, b = 1)
  # segments(s$dat$y, l, s$dat$y, u)

  b <- broom::tidyMCMC(m$model, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  # expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
})

test_that("mvt-nb2 model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  sigma <- 8
  df <- 5
  b0 <- 7
  n_draws <- 8
  gp_scale <- 1.6

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "nb2", b0 = b0)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER * 2, chains = CHAINS, obs_error = "nb2",
    estimate_df = FALSE, fixed_df_value = df,
    control = list(adapt_delta = 0.9), seed = SEED)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "nb2_phi[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[b$term == "B[1]", "estimate"], b0, tol = gp_scale * TOL)
})

test_that("mvt-gamma model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  sigma <- 0.3
  df <- 10
  b0 <- 2
  n_draws <- 8
  gp_scale <- 1.6

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "gamma", b0 = b0)
  # print(s$plot)

  m <- rrfield(y ~ 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER * 1.5, chains = CHAINS, obs_error = "gamma",
    estimate_df = FALSE, fixed_df_value = df, seed = SEED)
  print(m)

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "CV[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[b$term == "B[1]", "estimate"], b0, tol = gp_scale * TOL)
})


gp_sigma <- 0.2
sigma <- 0.1
df <- 10
gp_scale <- 1.2
n_draws <- 5
nknots <- 9
B <- c(0.5, 2.2, 3.8, 2.6, -0.9)
TOL <- 0.2 # %

test_that("mvt-norm estimates betas", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, B = B,
    X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  m <- rrfield(y ~ as.factor(time) - 1, data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df,
    prior_beta = rstanarm::student_t(3, 0, 10))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[grep("B\\[*", b$term), "estimate"], B, tol = B * TOL)
})
