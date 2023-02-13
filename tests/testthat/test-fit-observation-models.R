if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.25 # %
TOL_df <- .25 # %

nknots <- 10
gp_sigma <- 0.3

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
  n_draws <- 8
  gp_theta <- 1.6

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "nb2", B = b0
  )
  # print(s$plot)

  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, family = nbinom2(link = "log"),
    estimate_df = FALSE, fixed_df_value = df,
    control = list(adapt_delta = 0.9), seed = SEED
  )
  })

  p <- predict(m)

  b <- tidy(m, estimate.method = "median")
  # expect_equal(b[b$term == "nb2_phi[1]", "estimate", drop = TRUE], sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
  expect_equal(as.numeric(b[b$term == "B[1]", "estimate", drop = TRUE]), b0, tol = gp_theta * TOL)
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
  gp_theta <- 1.6

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "gamma", B = b0
  )
  # print(s$plot)

  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, family = Gamma(link = "log"),
    estimate_df = FALSE, fixed_df_value = df, seed = SEED
  )
  })

  p <- predict(m)

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "CV[1]", "estimate", drop = TRUE]), sigma, tol = sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
  expect_equal(as.numeric(b[b$term == "B[1]", "estimate", drop = TRUE]), b0, tol = gp_theta * TOL)
})

# ------------------------------------------------------
# a binomial observation model

test_that("mvt-binomial model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  nknots <- 12
  gp_sigma <- 1.2
  df <- 10
  b0 <- 0
  n_draws <- 15
  gp_theta <- 2.1

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "binomial", B = b0
  )
  # print(s$plot)
  # out <- reshape2::melt(s$proj)
  # names(out) <- c("time", "pt", "y")
  # out <- dplyr::arrange_(out, "time", "pt")
  # out$lon <- s$dat$lon
  # out$lat <- s$dat$lat
  # ggplot2::ggplot(out, ggplot2::aes(lon, lat, colour = plogis(y))) + ggplot2::geom_point() +
  #   ggplot2::facet_wrap(~time)

  suppressWarnings({
  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, family = binomial(link = "logit"),
    estimate_df = FALSE, fixed_df_value = df, seed = SEED
  )
  })
  # m
  #
  # p <- predict(m)
  # p <- predict(m, interval = "prediction")

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
})

# ------------------------------------------------------
# a poisson observation model

test_that("mvt-poisson model fits", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  nknots <- 12
  gp_sigma <- 0.8
  df <- 10
  b0 <- 3
  n_draws <- 15
  gp_theta <- 2.1

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    obs_error = "poisson", B = b0
  )
  # print(s$plot)
  # hist(s$dat$y)

  suppressWarnings({
  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, family = poisson(link = "log"),
    estimate_df = FALSE, fixed_df_value = df, seed = SEED
  )
  })
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(as.numeric(b[b$term == "gp_sigma", "estimate", drop = TRUE]), gp_sigma, tol = gp_sigma * TOL)
  expect_equal(as.numeric(b[b$term == "gp_theta", "estimate", drop = TRUE]), gp_theta, tol = gp_theta * TOL)
  expect_equal(as.numeric(b[b$term == "B[1]", "estimate", drop = TRUE]), b0, tol = gp_theta * TOL)
})
