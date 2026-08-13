data("so_mit_rbp")

test_that("the log-normalization works", {
  log_so <- sc_normalization(so_mit_rbp)
  no_log_so <- sc_normalization(so_mit_rbp, loga = FALSE)

  expect_true(all(log1p(no_log_so[["RNA"]]$data) == log_so[["RNA"]]$data))
})

test_that("the scale parameter works", {
  so_10k <- sc_normalization(so_mit_rbp, loga = FALSE)
  so_1 <- sc_normalization(so_mit_rbp, loga = FALSE, scale_factor = 1)

  expect_equal(so_1[["RNA"]]$data, so_10k[["RNA"]]$data / 10000)
})

test_that("imput check works", {
  expect_error(sc_normalization("hello"))
  expect_error(sc_normalization(so_mit_rbp, scale_factor = "hi"), "scale_factor must be a positive")
  expect_error(sc_normalization(so_mit_rbp, scale_factor = -10), "scale_factor must be a positive")
})
