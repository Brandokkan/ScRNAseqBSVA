data("so_mit_rbp")

test_that("cell_qc_vis errors on invalid plt_show", {
  expect_error(cell_qc_vis(so_mit_rbp, plt_show = "bad"), "invalid value")
})

test_that("cell_qc_vis returns a ggplot for vln", {
  pdf(NULL); on.exit(dev.off())

  result <- cell_qc_vis(so_mit_rbp, plt_show = "vln")
  expect_s3_class(result, "ggplot")
})

test_that("cell_qc_vis returns a patchwork for scat", {
  pdf(NULL); on.exit(dev.off())

  result <- cell_qc_vis(so_mit_rbp, plt_show = "scat")
  expect_s3_class(result, "patchwork")
})
