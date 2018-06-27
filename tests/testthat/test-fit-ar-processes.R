if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 600
CHAINS <- 2
SEED <- 9999
TOL <- 0.1 # %
TOL_df <- .25 # %
n_data_points <- 50

# ------------------------------------------------------
# a Gaussian observation model with random walk year effects

test_that("mvt-norm estimates random walk year effects", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED * 2)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 12
  nknots <- 5
  year_sigma <- 0.5
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, B = B,
    X = model.matrix(~a - 1, data.frame(a = gl(n_draws, 100)))
  )
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df, year_re = TRUE,
    prior_intercept = student_t(999, 0, 5), control = list(adapt_delta = 0.9),
    prior_rw_sigma = half_t(1e6, 0, 1)
  )
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_theta", "estimate"], gp_theta, tol = gp_theta * TOL)
  expect_equal(b[grep("yearEffects\\[*", b$term), "estimate"], B, tol = 0.1)
  expect_equal(b[grep("year_sigma", b$term), "estimate"], year_sigma, tol = 0.1)
})

# ------------------------------------------------------
# a Gaussian observation model with random walk year effects with covariate

test_that("mvt-norm estimates random walk year effects with covariate", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED * 2)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 12
  nknots <- 5
  year_sigma <- 0.5
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  cov_vec = rnorm(n_draws*100,0,1)
  model_matrix = model.matrix(~a - 1 + cov + cov2,
    data.frame(a = gl(n_draws, 100), cov=cov_vec, cov2=cov_vec^2))

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots,
    B = c(B, 3, -0.1),
    X = model_matrix)
  s$dat$cov = model_matrix[,"cov"]
  s$dat$cov2 = model_matrix[,"cov2"]
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) + geom_point()

  # include formula with the covariate and transformation
  m <- glmmfields(y ~ -1+ cov + cov2,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = FALSE, fixed_df_value = df, year_re = TRUE,
    prior_intercept = student_t(999, 0, 5), control = list(adapt_delta = 0.9),
    prior_rw_sigma = half_t(1e6, 0, 1)
  )
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "gp_theta", "estimate"], gp_theta, tol = gp_theta * TOL)
  expect_equal(b[grep("yearEffects\\[*", b$term), "estimate"], B, tol = 0.1)
  expect_equal(b[grep("year_sigma", b$term), "estimate"], year_sigma, tol = 0.1)
})

# ---------------
# AR process

test_that("mvt-norm estimates ar process", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 20
  nknots <- 7
  phi <- 0.5

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
    n_data_points = 100
  )
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) +
  # geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_ar = TRUE
  )
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "phi[1]", "estimate"], phi, tol = phi * TOL)
})

# -------------------
# AR + year random effects

test_that("mvt-norm estimates ar process *with* year random walk effects", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 20
  nknots <- 7
  phi <- 0.3
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  year_sigma <- 0.3
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
    B = B, X = model.matrix(~a - 1, data.frame(a = gl(n_draws, n_data_points))),
    n_data_points = n_data_points
  )
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) +
  # geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- glmmfields(y ~ 0,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    year_re = TRUE,
    estimate_ar = TRUE
  )
  m

  TOL <- 0.15
  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "year_sigma[1]", "estimate"], year_sigma, tol = year_sigma * 0.2)
  expect_equal(b[b$term == "phi[1]", "estimate"], phi, tol = phi * TOL)
  expect_equal(b[b$term == "phi[1]", "estimate"], phi, tol = phi * TOL)
})

# -------------------
# AR fixed + year random effects + covariate. These checks fail when phi
# is estimated, because the phi parameter in the spatial field is
# confounded with the AR(1) year random effects. Here phi is fixed.

test_that("mvt-norm estimates ar process *with* year random walk effects and covariate", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED)

  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 20
  nknots <- 7
  phi <- 0.7
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  year_sigma <- 0.3
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  cov_vec = rnorm(n_draws*n_data_points,0,1)
  model_matrix = model.matrix(~a - 1 + cov + cov2,
    data.frame(a = gl(n_draws, n_data_points), cov=cov_vec, cov2=cov_vec^2))

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
    B = c(B, 3, -0.1), X = model_matrix, n_data_points = n_data_points
  )
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) +
  # geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  # create dummy covariate
  s$dat$cov = model_matrix[,"cov"]
  s$dat$cov2 = model_matrix[,"cov2"]

  # include formula with the covariate and transformation
  m <- glmmfields(y ~ -1 + cov + cov2,
     data = s$dat, time = "time",
     lat = "lat", lon = "lon", nknots = nknots,
     iter = ITER, chains = CHAINS, seed = SEED,
     year_re = TRUE, fixed_phi_value = phi,
     estimate_ar = FALSE
   )
   m

  TOL <- 0.15
  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "gp_theta", "estimate"], gp_theta, tol = gp_theta * TOL)
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[grep("yearEffects\\[*", b$term), "estimate"], B, tol = 0.1)
  #expect_equal(b[b$term == "year_sigma[1]", "estimate"], year_sigma, tol = 0.1)
  #expect_equal(b[b$term == "phi[1]", "estimate"], phi, tol = phi * TOL)
})

# --------------------------
# global int + AR RF

test_that("mvt-norm estimates global int + AR RF", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED * 2)

  gp_sigma <- 0.2
  sigma <- 0.2
  df <- 1000
  gp_theta <- 1.8
  n_draws <- 20
  nknots <- 10
  phi <- 0.75
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 6
  year_sigma <- 0.0001
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
    B = B, X = model.matrix(~a - 1, data.frame(a = gl(n_draws, n_data_points))),
    n_data_points = n_data_points
  )
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) +
  # geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- glmmfields(y ~ 1,
    data = s$dat, time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    year_re = FALSE,
    control = list(adapt_delta = 0.95),
    estimate_ar = TRUE, prior_intercept = student_t(99, 0, 30)
  )
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "phi[1]", "estimate"], phi, tol = phi * TOL)
})

# --------------------------
# many ints + fixed AR

test_that("mvt-norm estimates many ints + fixed AR", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED * 5)

  gp_sigma <- 0.2
  sigma <- 0.2
  df <- 6
  gp_theta <- 1.8
  n_draws <- 20
  nknots <- 10
  phi <- 1
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 6
  year_sigma <- 0.4
  for (i in 2:length(B)) {
    B[i] <- B[i - 1] + rnorm(1, 0, year_sigma) # random walk
  }

  # plot(B)
  s <- sim_glmmfields(
    df = df, n_draws = n_draws, gp_theta = gp_theta,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, phi = phi,
    B = B, X = model.matrix(~a - 1, data.frame(a = gl(n_draws, n_data_points))),
    n_data_points = n_data_points
  )
  # print(s$plot)

  library(dplyr)
  means <- group_by(s$dat, time) %>% summarise(m = mean(y))
  plot(B, pch = 19)
  points(means$m, col = "red")

  library(ggplot2)
  ggplot(s$dat, aes(time, y)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2)) +
    geom_point(
      data = data.frame(B = B, time = seq_len(n_draws)),
      aes(x = time, y = B), inherit.aes = FALSE, col = "red"
    )

  m <- glmmfields(y ~ 0 + as.factor(time),
    data = s$dat,
    time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    fixed_df_value = 6, estimate_df = FALSE,
    estimate_ar = FALSE, fixed_phi_value = 1, prior_intercept = student_t(99, 0, 30)
  )
  m

  m2 <- glmmfields(y ~ -1,
    data = s$dat,
    time = "time",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    fixed_df_value = 6, estimate_df = FALSE,
    estimate_ar = FALSE, fixed_phi_value = 1, prior_intercept = student_t(99, 0, 30)
  )
  m2

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)

  B_hat <- subset(b, grepl("B", term))
  expect_equal(B, B_hat$estimate, tol = TOL)

  library(dplyr)
  q <- subset(b, grepl("B", term))
  qq <- group_by(s$dat, time) %>% summarise(m = mean(y))
  plot(q$estimate, col = "blue")
  points(seq_len(length(B)), B, pch = 19)
  points(qq$m, col = "red")
})
