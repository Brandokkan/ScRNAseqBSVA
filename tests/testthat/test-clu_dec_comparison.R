suex_dataset <- celldex::MouseRNAseqData()

data("pbmc3k", package = "pbmc3k.SeuratData")
pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
pbmc3k <- sc_normalization(pbmc3k)
pbmc3k_no_clusters <- pbmc3k
pbmc3k <- sc_clustering(pbmc3k)
pbmc3k_seu <- sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE)
pbmc3k_suex <- sc_deconvolute(pbmc3k, suex_dataset, diagnosis = FALSE)
pbmc3k_wrong_name <- pbmc3k_suex
names(pbmc3k_wrong_name@misc)[[1]] <- "wrong"


test_that("The function works properlly", {
  expect_no_error(clu_dec_comparison(pbmc3k_seu))
  expect_no_error(clu_dec_comparison(pbmc3k_suex))
})

test_that("A SeuratObject with no clusters is rejected", {
  expect_error(clu_dec_comparison(pbmc3k_no_clusters),
               "no clusters were found in the SeuratObject")
})

test_that("A SeuratObject with no predicted cell labels or with an invalid name
          for the element containing the predicted label is rejected", {
  expect_error(clu_dec_comparison(pbmc3k_wrong_name),
               "no predicted labels were found")
  expect_error(clu_dec_comparison(pbmc3k),
               "no predicted labels were found")
})
