csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")

test_that("final data makes sense for example dataset", {

  so_mit_rbp <- convert_to_so(csv_path)
  so_mit_rbp <- calculate_mt_rbp(so_mit_rbp, specie = "mouse")
  so_typo <- calculate_mt_rbp(so_mit_rbp, specie = "mose", overwrite = TRUE)

  expect_true(all(so_mit_rbp@meta.data[["percent.mt"]] > 0))
  expect_false(any(so_typo@meta.data[["percent.mt"]] > 0))

  expect_true(all(so_mit_rbp@meta.data[["percent.rbp"]] > 0))
  expect_false(any(so_typo@meta.data[["percent.rbp"]] > 0))
})

test_that("An error is called when overwrite is set to TRUE wile no data is present", {

  so_no_mit_rbp <- convert_to_so(csv_path)

  expect_error(calculate_mt_rbp(so_no_mit_rbp, overwrite = TRUE), "Error. Overwrite was set to TRUE while one or both metadata are absent")
})

test_that("The method rejects a non-SeuratObject", {
  expect_error(calculate_mt_rbp("asdad"))
})
