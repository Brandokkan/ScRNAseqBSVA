#' @describeIn cell_qc_vis method to visualize (and calculate if necessary)
#'  quality control parameters.
#' @importFrom Seurat VlnPlot
#' @importFrom Seurat FeatureScatter
#' @export
setMethod("cell_qc_vis", signature(so = "Seurat"), function(so, specie = "human",
                                                            mit_pat = "^MT-",
                                                            rib_pat = "^RP[LS]",
                                                            plt_show = c("scat", "vln"),
                                                            ...) {
  so <- calculate_mt_rbp(so, specie, mit_pat, rib_pat)

  if ("vln" %in% plt_show) {
  vln_plt <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"),
             ncol = 4, ...)
  vln_plt
  }

  if ("scat" %in% plt_show) {
  plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.rbp")
  scat_plot <- plot1 + plot2 + plot3
  scat_plot
  }

})
