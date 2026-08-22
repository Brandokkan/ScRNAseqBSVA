#' @describeIn clu_dec_comparison method to compare the clustering and the
#' deconvolution of \code{SeuratObjects}
#' @importFrom pheatmap pheatmap
#' @export
setMethod("clu_dec_comparison", signature(so = "Seurat"), function(so) {

  if (!"seurat_clusters" %in% colnames(so@meta.data)) {
    stop("Error: no clusters were found in the SeuratObject. Make sure
         sc_clustering was used on the object first.",
         call. = FALSE)
  }

  if ("predicted.id" %in% colnames(so@meta.data)) {
    cell_type <- so$predicted.id
  } else if ("predictions" %in% names(so@misc)) {
    cell_type <- so@misc$predictions$labels
  } else {
    stop("Error: no predicted labels were found following standard naming
         conventions of sc_deconvolute",
         call. = FALSE)
  }

  comparison_table <- table(
    cluster    = so$seurat_clusters,
    cell_type
  )

  comparison_table

  proportions_table <- prop.table(comparison_table, margin = 1)
  proportions_table

  pheatmap(
    proportions_table,
    main = "Proportion of predicted cell types per cluster"
  )
})
