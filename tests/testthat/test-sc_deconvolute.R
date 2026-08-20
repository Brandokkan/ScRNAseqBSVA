data("so_norm_log")
data("so_mit_rbp")
pbmc3k <- SeuratData::LoadData("pbmc3k")
pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
pbmc3k <- sc_normalization(pbmc3k)
pbmc3k <- sc_clustering(pbmc3k)
suex_dataset <- celldex::MouseRNAseqData()
# so_suex <- Seurat::as.Seurat(as(suex_dataset, "SingleCellExperiment"), counts = NULL)
scex <- scRNAseq::ZeiselBrainData()
scex <- scuttle::logNormCounts(scex)
# so_scex <- as.Seurat(scex, data = NULL)

test_that("method for SummarizedExperiment and SingleCellExperiment works and gives
          sound results", {
  suex_labels <- SummarizedExperiment::colData(suex_dataset)$label.main
  suex_labels_p <- sc_deconvolute(suex_dataset, suex_dataset, diagnosis = FALSE)$predictions$labels
  concordance <- mean(suex_labels_p == suex_labels)
  expect_gt(concordance, 0.85)

  scex_labels <- scex$level1class
  scex_labels_p <- sc_deconvolute(scex, scex, diagnosis = FALSE, lab_name = "level1class")$predictions$labels
  concordance <- mean(scex_labels == scex_labels_p)
  expect_gt(concordance, 0.85)
})

test_that("method for SeuratObject works and gives sound results", {
  pbm_labels <- pbmc3k@meta.data$seurat_annotations
  pbm_labels_p <- sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE)@meta.data$predicted.id
  concordance <- mean(pbm_labels == pbm_labels_p, na.rm = TRUE)
  expect_gt(concordance, 0.85)
})
