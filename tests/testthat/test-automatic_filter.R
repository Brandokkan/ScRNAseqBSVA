data("so_mit_rbp")

test_that("returned cell number make sense for worst and best case scenario", {
  filt_so <- automatic_filter(so_mit_rbp)
  expect_true(ncol(filt_so) <= ceiling(0.9 * ncol(so_mit_rbp)) & ncol(filt_so) >= floor(0.85 * ncol(so_mit_rbp))
              )
})

test_that("the input check for the percentiles works", {
  expect_error(automatic_filter(so_mit_rbp, up_per = 0.3, low_per = 0.9),
               "low_per was set as higher")
})

test_that("invalid inputs give errors", {
  expect_error(automatic_filter("hello"))
  expect_error(automatic_filter(so_mit_rbp, up_per = "hu"))
  expect_error(automatic_filter(so_mit_rbp, low_per = "hu"))
  expect_error(automatic_filter(so_mit_rbp, mit_per = "hu"))
})
