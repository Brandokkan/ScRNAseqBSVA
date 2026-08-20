#' @describeIn sc_deconvolute method used to predict cell types of a SeuratObject
#' from a reference SummarizedExperiment dataset
#' @importFrom Seurat as.SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colnames
#' @importFrom SingleR SingleR
#' @export
setMethod("sc_deconvolute", signature(so = "Seurat", ref = "SummarizedExperiment"),
          function(so, ref, diagnosis = TRUE) {
  sce <- as.SingleCellExperiment(so)

  if (!"logcounts" %in% assayNames(sce)) {
    stop("Error: logcounts is not present in the query. Make sure
         the query was LogNormalized", call. = FALSE )
  }

  if (!lab_name %in% colnames(colData(ref))) {
    stop("Error: lab_name value is not in colData. Try looking in colData of ref
         and see what is the name of the vector containing the cell labels",
         call. = FALSE)
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
#' @importFrom SummarizedExperiment assayNames
#' @importFrom SingleR SingleR
#' @importFrom SummarizedExperiment colData
#' @importFrom SummarizedExperiment colnames
#' @export
setMethod("sc_deconvolute", signature(so = "SummarizedExperiment", ref = "SummarizedExperiment"),
          function(so, ref, diagnosis = TRUE, lab_name = "label.main") {
  sce <- so

  if (!"logcounts" %in% assayNames(sce)) {
    stop("Error: logcounts is not present in the query. Make sure
          the query was LogNormalized", call. = FALSE )
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


#' @describeIn sc_deconvolute method used to predict cell types in a
#' SeuratObject from a reference SeuratObject
#' @importFrom SeuratObject Layers
#' @importFrom Seurat FindTransferAnchors
#' @importFrom Seurat TransferData
#' @importFrom Seurat AddMetaData
#' @export
setMethod("sc_deconvolute", signature(so = "Seurat", ref = "Seurat"),
          function(so, ref, diagnosis = TRUE, lab_name = "seurat_annotations") {

            if (!"data" %in% Layers(so[["RNA"]])) {
              stop("Error: data is not present in the query SeuratObject. Make sure
                   the SeuratObject was LogNormalized", call. = FALSE )
            }

            if (!lab_name %in% colnames(ref@meta.data)) {
              stop("Error: lab_name value is not in meta.data. Try looking in meta.data of ref
                    and see what is the name of the vector containing the cell labels",
                   call. = FALSE)
            }

            anchors <- FindTransferAnchors(reference = ref,
                                           query = so)

            predictions <- TransferData(anchorset = anchors, refdata = ref[[lab_name, drop = TRUE]])
            so <- AddMetaData(object = so, metadata = predictions)

            if (diagnosis) {
              prediction_diagnostics(so)
            }

            so
          })
