if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 400
CHAINS <- 2
SEED <- 9999
TOL <- 0.15 # %
TOL_df <- 0.15 # %

# --------------------------
# fixed effects + fixed AR with and without demeaning

test_that("check demeaning", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(SEED*1234)

  gp_sigma <- 0.6
  sigma <- 0.2
  df <- 2.5
  gp_scale <- 1
  n_draws <- 20
  nknots <- 10
  ar <- 0.95
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 1
  year_sigma <- 0.8
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  plot(B)

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
    B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))))
  print(s$plot)

  # look at the auto regressive spatial field (de-meaned)
  library(dplyr)
  library(ggplot2)
  group_by(s$dat, time) %>%
    mutate(z = y - mean(y)) %>% # A quick hack, this should be offset by a year
    ggplot(aes(lon, lat, color = z)) +
    geom_point() +
    scale_colour_gradient2() +
    facet_wrap(~time)

  # Look at the true auto regressive knot spatial field
  data.frame(reshape2::melt(s$re_knots) %>% arrange(Var1, Var2),
    as.data.frame(s$knots)) %>%
    ggplot(aes(lon, lat, colour = value)) +
    geom_point() +
    scale_colour_gradient2() +
    facet_wrap(~Var1)

  means <- group_by(s$dat, time) %>% summarise(m = mean(y))
  plot(B, pch = 19)
  points(means$m, col = "red")

  ggplot(s$dat, aes(time, y)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.2)) +
    geom_point(data = data.frame(B = B, time = seq_len(n_draws)),
      aes(x = time, y = B), inherit.aes = FALSE, col = "red")

  m <- rrfield(y ~ 1, data = s$dat,
    time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    fixed_df_value = df, estimate_df = T,
    estimate_ar = TRUE,
    year_re = TRUE,
    prior_intercept = student_t(99, 0, 20),
    prior_beta = student_t(99, 0, 20))
  m

  b <- tidy(m, estimate.method = "median")
  expect_equal(b[b$term == "sigma[1]", "estimate"], sigma, tol = sigma * TOL)
  expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
  expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = gp_sigma * TOL)
  B_hat <- subset(b, grepl("yearEffects", term))
  expect_equal(B, B_hat$estimate, tol = TOL)

})
