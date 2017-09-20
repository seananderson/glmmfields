if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.2 # %
TOL_df <- .25 # %

gp_sigma <- 0.2
sigma <- 0.1
df <- 1000
gp_theta <- 1.2
n_draws <- 15
nknots <- 8
n_data_points <- 50

test_that("predict.glmmfields works", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  set.seed(SEED)

  s <- sim_glmmfields(df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, n_data_points = n_data_points)

  m <- glmmfields(y ~ 0, data = s$dat, time = "time",
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
  expect_gte(cor(s$dat$y, p_newdata2$estimate), 0.75)

  # with a subset of data
  random_subset <- sample(seq_len(nrow(s$dat)), size = 200)
  p_newdata <- predict(m, newdata = s$dat[random_subset, ])
  plot(s$dat$y[random_subset], p_newdata$estimate)
  expect_gte(cor(s$dat$y[random_subset], p_newdata$estimate), 0.75)

  nd <- s$dat
  nd$y <- NULL
  p <- predict(m, newdata = nd)
})
