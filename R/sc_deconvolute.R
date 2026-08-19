#' @describeIn sc_deconvolute method used to predict cell types of a SeuratObject
#' from a reference SummarizedExperiment dataset
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


#' @describeIn sc_deconvolute method used to predict cell types in a
#' SummarizedExperiment or SingleCellExperiment from a reference
#' SummarizedExperiment or SingleCellExperiment dataset
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleR SingleR
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colnames
#' @export
setMethod("sc_deconvolute", signature(so = "SummarizedExperiment", ref = "SummarizedExperiment"),
          function(so, ref, diagnosis = TRUE, lab_name = "label.main") {
  sce <- so

  if (!"logcounts" %in% assayNames(sce)) {
    stop("Error: logcounts is not present in the SingleCellExperiment. Make sure
          the SeuratObject was LogNormalized", call. = FALSE )
  }

  if (!lab_name %in% colnames(colData(ref))) {
    stop("Error: lab_name value is not in colData. Try looking in colData of ref
         and see what is the name of the vector containing the cell labels",
         call. = FALSE)
  }

  predictions <- SingleR(
    test    = sce,
    ref     = ref,
    labels  = ref[[lab_name]],
    de.method = "classic"
  )

  so$"predictions" <- predictions

  if (diagnosis) {
    prediction_diagnostics(so)
  }

  so
})
