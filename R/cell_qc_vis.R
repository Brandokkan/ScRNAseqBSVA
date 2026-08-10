#' @describeIn cell_qc_vis method to visualize (and calculate if necessary)
#'  quality control parameters.
#' @export
setMethod("cell_qc_vis", signature(so = "Seurat"), function(so, specie = "human",
                                                            mit_pat = "^MT-",
                                                            rib_pat = "^RP[LS]",
                                                            ...) {
  if (!"percent.mt" %in% colnames(so@meta.data)) {

  }
})
