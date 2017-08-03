if (interactive()) options(mc.cores = parallel::detectCores())

set.seed(42)

s <- sim_glmmfields(df = 1000, n_draws = 2, gp_scale = 1.5,
  gp_sigma = 0.3, sd_obs = 0.1, n_knots = 10, n_data_points = 50)

test_that("Stations in second time slice can be in different order from first-time slice", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  d <- s$dat
  d$ID <- seq_len(nrow(d))
  m <- glmmfields(y ~ 0, data = d, time = "time",
    lat = "lat", lon = "lon", nknots = 10,
    iter = 400, chains = 2, seed = 1)
  d$pred <- predict(m)$estimate

  d2 <- d
  d2[d2$time == 2, ] <- d2[d2$time == 2, ][sample(seq_len(50), size = 50), ] # scramble time 2
  m2 <- glmmfields(y ~ 0, data = d2, time = "time",
    lat = "lat", lon = "lon", nknots = 10,
    iter = 400, chains = 2, seed = 2)
  d2$pred <- predict(m2)$estimate
  d2 <- dplyr::arrange(d2, ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred)
})

test_that("Stations in second time slice introduce new stations", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  d <- s$dat
  d$ID <- seq_len(nrow(d))
  m <- glmmfields(y ~ 0, data = d, time = "time",
    lat = "lat", lon = "lon", nknots = 10,
    iter = 800, chains = 2, seed = 1)

  d2 <- d[-c(2, 10, 25), ]
  m2 <- glmmfields(y ~ 0, data = d2, time = "time",
    lat = "lat", lon = "lon", nknots = 10,
    iter = 800, chains = 2, seed = 1)
  d2$pred <- predict(m2)$estimate
  d$pred <- predict(m)$estimate
  d <- dplyr::filter(d, ID %in% d2$ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred, tolerance = .03)
})

test_that("Ordering of time slices doesn't matter if stations aren't always presents ", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(42)
  s2 <- sim_glmmfields(df = 1000, n_draws = 2, gp_scale = 1.5,
    gp_sigma = 0.3, sd_obs = 0.1, n_knots = 2, n_data_points = 3)
  d <- s2$dat
  d$ID <- seq_len(nrow(d))
  d <- d[-c(2), ]

  m <- glmmfields(y ~ 0, data = d, time = "time",
    lat = "lat", lon = "lon", nknots = 2,
    iter = 800, chains = 2, seed = 1, cores = 1)
  sd <- m$stan_data

  d2 <- rbind(d[d$time == 2,], d[d$time == 1,])
  m2 <- glmmfields(y ~ 0, data = d2, time = "time",
    lat = "lat", lon = "lon", nknots = 2,
    iter = 800, chains = 2, seed = 1, cores = 1)
  sd2 <- m2$stan_data
  d2$pred <- predict(m2)$estimate
  d$pred <- predict(m)$estimate
  d2 <- dplyr::arrange(d2, ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred, tolerance = .01)
})
