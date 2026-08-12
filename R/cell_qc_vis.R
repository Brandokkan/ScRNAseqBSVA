#' @describeIn cell_qc_vis method to visualize (and calculate if necessary)
#'  quality control parameters.
#' @importFrom Seurat VlnPlot
#' @importFrom Seurat FeatureScatter
#' @export
setMethod("cell_qc_vis", signature(so = "Seurat"), function(so, plt_show,
                                                            specie = "human",
                                                            mit_pat = "^MT-",
                                                            rib_pat = "^RP[LS]",
                                                            overwrite = FALSE,
                                                            ...) {
  so <- calculate_mt_rbp(so, specie, mit_pat, rib_pat, overwrite)

  if ("vln" == plt_show) {
  vln_plt <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp"),
             ncol = 4, ...)
  print(vln_plt)
  } else if ("scat" == plt_show) {
  plot1 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot3 <- FeatureScatter(so, feature1 = "nCount_RNA", feature2 = "percent.rbp")
  scat_plot <- plot1 + plot2 + plot3
  print(scat_plot)
  } else {
    stop("an invalid value was suplied to plt_show")
  }

})
