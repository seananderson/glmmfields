test_that("families are parsed", {
  expect_equal(check_family(nbinom2(link = "log"))$family, "nbinom2")
  expect_equal(check_family(nbinom2(link = "log"))$link, "log")
  expect_equal(check_family(lognormal(link = "log"))$family, "lognormal")
  expect_equal(check_family(nbinom2(link = "log"))$link, "log")

  expect_error(check_family(gaussian(link = "aaa")))
  expect_error(check_family(nbinom2(link = "aaa")))
  expect_error(check_family(lognormal(link = "aaa")))

  expect_error(check_family(quasibinomial(link = "logit")))

  check_family(binomial(link = logit))
  check_family(nbinom2(link = log))
  check_family(lognormal(link = log))
})
