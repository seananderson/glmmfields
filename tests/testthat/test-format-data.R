test_that("format_data() formats some data", {
  s <- sim_glmmfields(n_data_points = 50, n_knots = 5, n_draws = 2)
  mf <- model.frame(y ~ 1, s$dat)
  X <- model.matrix(y ~ 1, mf)
  y <- model.response(mf, "numeric")
  f <- format_data(s$dat, y = y, X = X, time = "time", nknots = 5)

  expect_equal(nrow(f$knots), 5L)
  expect_equal(dim(f$spatglm_data$distKnots), c(5L, 5L))
})
