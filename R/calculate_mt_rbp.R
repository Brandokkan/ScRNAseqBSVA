#' @describeIn calculate_mt_rbp method to calculate the mitochondrial and ribosomal
#'  percentages.
#' @importFrom Seurat PercentageFeatureSet
#' @export
setMethod("calculate_mt_rbp", signature(so = "Seurat"), function(so, specie = "human",
                                                            mit_pat = "^MT-",
                                                            rib_pat = "^RP[LS]",
                                                            ...) {
  added <- FALSE

  if (!"percent.mt" %in% colnames(so@meta.data)) {
    if (specie == "human") {
      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    } else if (specie == "mouse") {
      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
    } else {
      so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = mit_pat)
    }
    added <- TRUE
    print("Added % of mitochondrial genes to SeuratObject")
  }

  if (!"percent.rbp" %in% colnames(so@meta.data)) {
    if (specie == "human") {
      so[["percent.rbp"]] <- PercentageFeatureSet(so, pattern = "^RP[LS]")
    } else if (specie == "mouse") {
      so[["percent.rbp"]] <- PercentageFeatureSet(so, pattern = "^Rp[ls]")
    } else {
      so[["percent.rbp"]] <- PercentageFeatureSet(so, pattern = rib_pat)
    }
    added <- TRUE
    print("Added % of ribosomal protein genes to SeuratObject")
  }

  if (!added) {
    print("The % of mitochondrial and ribosomal genes are already present.")
  }
  so
})
