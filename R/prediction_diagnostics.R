#' @describeIn prediction_diagnostics method to plot different visualization
#' tools useful to asses the accuracy of the predictions
#' @importFrom SingleR plotScoreHeatmap
#' @importFrom SingleR plotDeltaDistribution
#' @export
setMethod("prediction_diagnostics", singature(so = "Seurat"), function(so) {
  if ("predictions" %in% colnames(so@meta.data)) {
    predictions <- so$predictions
  } else {
    stop("Error: no prediction metadata was present in the SeuratObject.
         Make sure sc_deconvolute was used first on the SeuratObject.",
         call. = FALSE)
  }

  # Heatmap of per-cell scores across reference labels — look for clean diagonal blocks
  plotScoreHeatmap(predictions)

  # Cells where the top label wasn't clearly better than the runner-up get pruned to NA
  summary(is.na(predictions$pruned.labels))

  # Delta distribution helps spot ambiguous calls
  plotDeltaDistribution(predictions, ncol = 3)
})
