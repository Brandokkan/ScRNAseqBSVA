skip_if_not_installed("pbmc3k.SeuratData")
skip_if_not_installed("celldex")

e <- new.env()
data("pbmc3k", package = "pbmc3k.SeuratData", envir = e)
pbmc3k <- e$pbmc3k
pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
pbmc3k <- sc_normalization(pbmc3k)
pbmc3k <- sc_clustering(pbmc3k)
suex_dataset <- celldex::MouseRNAseqData()

test_that("prediction_diagnostics works for SingleR (SummarizedExperiment) predictions", {
    suex_pred <- sc_deconvolute(suex_dataset, suex_dataset, diagnosis = FALSE)
    result <- prediction_diagnostics(suex_pred)
    expect_true(is.table(result))
})

test_that("prediction_diagnostics works for TransferData (Seurat-Seurat) predictions", {
    pbmc3k_pred <- sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE)
    result <- prediction_diagnostics(pbmc3k_pred)
    expect_true(is.table(result))
})

test_that("prediction_diagnostics works for SingleR (SummarizedExperiment - Seurat) predictions", {
    suex_pred <- sc_deconvolute(pbmc3k, suex_dataset, diagnosis = FALSE)
    result <- prediction_diagnostics(suex_pred)
    expect_true(is.table(result))
})

test_that("prediction_diagnostics errors when no prediction metadata is present", {
    expect_error(prediction_diagnostics(pbmc3k), "no prediction metadata")
})
