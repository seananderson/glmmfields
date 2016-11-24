test_that("priors parse", {
  expect_equal(parse_t_prior(student_t(3, 0, 1)), c(3, 0, 1))
  expect_error(student_t(-99, 0, 1))
  expect_error(student_t(3, 0, -99))
  expect_warning(half_t(3, -99, 1))
})
