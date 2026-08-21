#' @describeIn prediction_diagnostics method to plot different visualization
#' tools useful to asses the accuracy of the predictions in SeuratObjects
#' @importFrom SingleR plotScoreHeatmap
#' @importFrom stats heatmap
#' @export
setMethod("prediction_diagnostics", signature(so = "Seurat"), function(so) {
  if (!is.null(so@misc$predictions)) {
    predictions <- so@misc$predictions

    # Heatmap of per-cell scores across reference labels — look for clean diagonal blocks
    plot_prediction_scores(predictions)

    # Cells where the top label wasn't clearly better than the runner-up get pruned to NA
    summary(is.na(predictions$pruned.labels))
  } else if ("predicted.id" %in% colnames(so@meta.data)) {
    score_cols <- grep("^prediction\\.score\\.", colnames(so@meta.data), value = TRUE)
    score_cols <- setdiff(score_cols, "prediction.score.max")

    if (length(score_cols) == 0) {
      stop("Error: no prediction.score columns were found in the meta.data.
           Make sure sc_deconvolute was used with a SeuratObject reference
           and query.",
           call. = FALSE)
    }

    scores <- as.matrix(so@meta.data[, score_cols, drop = FALSE])
    colnames(scores) <- sub("^prediction\\.score\\.", "", score_cols)

    # Heatmap of per-cell scores across reference labels, cells ordered by
    # their predicted label
    heatmap(t(scores[order(so$predicted.id), , drop = FALSE]), Colv = NA,
            scale = "none", labCol = FALSE, main = "Prediction scores per cell")

    # Cells whose winning label score didn't clear 0.5 are considered low-confidence
    summary(so$prediction.score.max < 0.5)
  } else {
    stop("Error: no prediction metadata was present in the SeuratObject.
         Make sure sc_deconvolute was used first on the SeuratObject.",
         call. = FALSE)
  }

})


#' @describeIn prediction_diagnostics method to plot different visualization
#' tools useful to asses the accuracy of the predictions in SUmmarizedExperiments
#' @importFrom SingleR plotScoreHeatmap
#' @importFrom SummarizedExperiment colData
#' @export
setMethod("prediction_diagnostics", signature(so = "SummarizedExperiment"), function(so) {
  if ("predictions" %in% colnames(colData(so))) {
    predictions <- so$predictions

    # Heatmap of per-cell scores across reference labels — look for clean diagonal blocks
    plot_prediction_scores(predictions)

    # Cells where the top label wasn't clearly better than the runner-up get pruned to NA
    summary(is.na(predictions$pruned.labels))
  } else {
    stop("Error: no prediction metadata was present in the SeuratObject.
         Make sure sc_deconvolute was used first on the SeuratObject.",
         call. = FALSE)
  }

})
