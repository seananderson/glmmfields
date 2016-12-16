if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.2 # %
TOL_df <- .25 # %

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

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_scale", "estimate"], gp_scale, tol = gp_scale * TOL)
  expect_equal(b[grep("B\\[*", b$term), "estimate"], B, tol = 0.05)
})
