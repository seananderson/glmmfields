test_that("families are parsed", {
  expect_equal(nbinom2(link = "log")$family, "nbinom2")
  expect_equal(nbinom2(link = "log")$link, "log")
  expect_equal(lognormal(link = "log")$family, "lognormal")
  expect_equal(nbinom2(link = "log")$link, "log")
})
