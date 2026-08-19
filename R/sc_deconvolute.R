#' @describeIn sc_deconvolute method used to predict cell types from a reference
#' SummarizedExperiment dataset
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleR SingleR
#' @export
setMethod("sc_deconvolute", signature(so = "Seurat", ref = "SummarizedExperiment"),
          function(so, ref, diagnosis = TRUE) {
  sce <- as.SingleCellExperiment(so)

  if (!"logcounts" %in% assayNames(sce)) {
    stop("Error: logcounts is not present in the SingleCellExperiment. Make sure
         the SeuratObject was LogNormalized", call. = FALSE )
  }

  predictions <- SingleR(
    test    = sce,
    ref     = ref,
    labels  = ref$label.main,
    de.method = "classic"
  )

  so$"predictions" <- predictions

  if (diagnosis) {
    prediction_diagnostics(so)
  }

  so
})
