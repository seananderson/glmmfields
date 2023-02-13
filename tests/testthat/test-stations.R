set.seed(42)

s <- sim_glmmfields(
  df = 1000, n_draws = 2, gp_theta = 1.5,
  gp_sigma = 0.3, sd_obs = 0.1, n_knots = 8, n_data_points = 30
)

test_that("Stations in second time slice can be in different order from first time slice", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  d <- s$dat
  d$ID <- seq_len(nrow(d))
  suppressWarnings({
    m <- glmmfields(y ~ 0,
      data = d, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 400, chains = 2, seed = 1
    )
  })
  d$pred <- predict(m)$estimate

  d2 <- d
  d2[d2$time == 2, ] <- d2[d2$time == 2, ][sample(seq_len(30), size = 30), ] # scramble time 2
  suppressWarnings({
    m2 <- glmmfields(y ~ 0,
      data = d2, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 400, chains = 2, seed = 1
    )
  })
  d2$pred <- predict(m2)$estimate
  d2 <- dplyr::arrange(d2, ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred, tolerance = 0.000001)
})

test_that("Stations in second time slice introduce new stations", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  d <- s$dat
  d$ID <- seq_len(nrow(d))
  suppressWarnings({
    m <- glmmfields(y ~ 0,
      data = d, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 800, chains = 2, seed = 1
    )
  })

  d2 <- d[-c(2, 10), ]
  suppressWarnings({
    m2 <- glmmfields(y ~ 0,
      data = d2, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 800, chains = 2, seed = 1
    )
  })
  d2$pred <- predict(m2)$estimate
  d$pred <- predict(m)$estimate
  d <- dplyr::filter(d, ID %in% d2$ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred, tolerance = .02)
})

test_that("Ordering of time slices doesn't matter if stations aren't always present", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  d <- s$dat
  d$ID <- seq_len(nrow(d))
  d <- d[-c(2, 10), ]

  suppressWarnings({
    m <- glmmfields(y ~ 0,
      data = d, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 800, chains = 2, seed = 1, cores = 1
    )
  })
  sd <- m$stan_data

  d2 <- rbind(d[d$time == 2, ], d[d$time == 1, ])
  suppressWarnings({
    m2 <- glmmfields(y ~ 0,
      data = d2, time = "time",
      lat = "lat", lon = "lon", nknots = 8,
      iter = 800, chains = 2, seed = 1, cores = 1
    )
  })
  sd2 <- m2$stan_data
  d2$pred <- predict(m2)$estimate
  d$pred <- predict(m)$estimate
  d2 <- dplyr::arrange(d2, ID)

  plot(d2$pred, d$pred)
  expect_equal(d2$pred, d$pred, tolerance = .01)
})
