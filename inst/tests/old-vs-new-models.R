set.seed(42)
if (interactive()) options(mc.cores = parallel::detectCores())

ITER <- 500
CHAINS <- 2
TOL <- 0.2

test_new <- function(i) {
  gp_sigma <- 0.2
  sigma <- 0.1
  df <- 8
  gp_scale <- 1.8
  n_draws <- 18
  nknots <- 10
  ar <- 0.9
  B <- vector(mode = "double", length = n_draws)
  B[1] <- 0
  year_sigma <- 0.2
  for (i in 2:length(B)) {
    B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
  }

  s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
    gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
    B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))),
    obs_error = "gamma")
  # print(s$plot)
  # library(ggplot2); ggplot(s$dat, aes(time, y)) +
    # geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

  m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
    lat = "lat", lon = "lon", nknots = nknots,
    iter = ITER, chains = CHAINS, seed = SEED,
    estimate_df = TRUE, year_re = TRUE, obs_error = "gamma",
    estimate_ar = TRUE, algorithm = "sampling")
  list(model = m, B = B)
}

x <- lapply(seq_len(4), test_new)

library(purrr)
library(dplyr)
par(mfrow = c(2,2))
out <- map_dbl(x, function(y) {
  b <- tidy(y$model, estimate.method = "median")
  est <- filter(b, grepl("yearEffects", term))
  # mean((est$estimate - y$B)^2)
  plot(est$estimate, y$B)
  mean(est$estimate - y$B)
})



# --------------------------------
# old vs new
library(testthat)

set.seed(SEED*42)

ITER <- 600
gp_sigma <- 0.3
sigma <- 1
df <- 15
gp_scale <- 0.07
n_draws <- 18
nknots <- 10
ar <- 1
B <- vector(mode = "double", length = n_draws)
B[1] <- 3
year_sigma <- 0.14
for (i in 2:length(B)) {
  B[i] <- B[i-1] + rnorm(1, 0, year_sigma) # random walk
}

s <- sim_rrfield(df = df, n_draws = n_draws, gp_scale = gp_scale,
  gp_sigma = gp_sigma, sd_obs = sigma, n_knots = nknots, ar = ar,
  B = B, X = model.matrix(~ a - 1, data.frame(a = gl(n_draws, 100))),
  obs_error = "gamma")
print(s$plot)
library(ggplot2); ggplot(s$dat, aes(time, y)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.2))

# --------------
# new:
m <- rrfield(y ~ 0, data = s$dat, time = "time", station = "station_id",
  lat = "lat", lon = "lon", nknots = nknots,
  iter = ITER*1.5, chains = 3, seed = SEED,
  estimate_df = TRUE, year_re = TRUE, obs_error = "gamma",
  estimate_ar = TRUE, algorithm = "sampling", control = list(adapt_delta = 0.95))
m

# --------------
# old:
Y <- s$dat$y
old_dat <- list(
  nKnots = nknots,
  nLocs = nrow(s$dist_knots21_sq),
  nT = length(unique(s$dat$time)),
  N = length(Y), stationID = s$dat$station_id,
  yearID = s$dat$time, y = Y, x = rep(0,length(Y)),
  distKnotsSq = s$dist_knots_sq, distKnots21Sq = s$dist_knots21_sq)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
m_old <- stan("inst/tests/mvtGamma_estSigma_index_cov_yr_ar1.stan", data = old_dat,
  iter = ITER*1.2, chains = 3, control = list(adapt_delta = 0.9))
m_old

b_new <- tidy(m, estimate.method = "median")
b_old <- broom::tidyMCMC(m_old, estimate.method = "median")

check_est <- function(old, new, tol = 0.1) {
  expect_equal(
    b_old[b_old$term == old, "estimate"],
    b_new[b_new$term == new, "estimate"], tol = tol)
}

# check_est("gp_scale", "gp_scale")
check_est("ar", "ar[1]")
check_est("CV", "CV[1]")
check_est("scaledf", "df[1]")
# check_est("gp_sigmaSq", "gp_sigma")
check_est("yearEffects[1]", "yearEffects[1]")
check_est("yearEffects[2]", "yearEffects[2]")
check_est("yearEffects[10]", "yearEffects[10]")

old <- filter(b_old, grepl("spatialEffects", term))
new <- filter(b_new, grepl("spatialEffects", term))
mean(old$estimate - new$estimate)
plot(old$estimate[1:100], new$estimate[1:100])

check_est("spatialEffects[3,2]", "spatialEffectsKnots[3,2]")

check_est("spatialEffects[3,2]", "spatialEffectsKnots[3,2]")

## expect_equal(b[b$term == "gp_sigma", "estimate"], gp_sigma, tol = gp_sigma * TOL)
## expect_equal(b[b$term == "df[1]", "estimate"], df, tol = df * TOL_df)
## expect_equal(b[b$term == "year_sigma[1]", "estimate"], year_sigma, tol = year_sigma * TOL_df)
## expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
## expect_equal(b[b$term == "ar[1]", "estimate"], ar, tol = ar * TOL)
