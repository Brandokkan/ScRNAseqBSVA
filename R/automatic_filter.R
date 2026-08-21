#' @describeIn automatic_filter method to filter the cells
#' @importFrom stats quantile
#' @export
setMethod("automatic_filter", signature(so = "Seurat"), function(so, up_per = 0.9,
                                                                 low_per = 0.1,
                                                                 mit_per = 0.95,
                                                                 ...) {
  if (low_per >= up_per) {
    stop("low_per was set as higher or equal compared to up_per",
         call. = FALSE)
  }

  per_feat <- quantile(so@meta.data$nFeature_RNA, probs = c(low_per, up_per))
  per_mit <- quantile(so@meta.data$percent.mt, probs = c(mit_per))

  low_per_val <- per_feat[[1]]
  up_per_val <- per_feat[[2]]
  mit_val <- per_mit[[1]]

  subset(so, subset = nFeature_RNA > low_per_val & nFeature_RNA < up_per_val & percent.mt < mit_val,
         ...)
})
