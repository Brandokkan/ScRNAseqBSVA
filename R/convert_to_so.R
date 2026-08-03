#' @describeIn convert_to_so Method for a \code{data.frame} of counts
#'   (genes x cells).
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod("convert_to_so", signature(sc_data = "data.frame"), function(sc_data, ...) {
  Seurat::CreateSeuratObject(counts = sc_data, ...)
})

#' @describeIn convert_to_so Method for a file path to a supported
#'   single-cell data file (e.g. a 10x Genomics directory or \code{.h5} file).
#' @export
setMethod("convert_to_so", signature(sc_data = "character"), function(sc_data, ...) {
  data <- read_sc_file(sc_data)
  Seurat::CreateSeuratObject(counts = data, ...)
})

#' @describeIn convert_to_so Method for a sparse counts matrix of class
#'   \code{dgCMatrix}.
#' @importClassesFrom Matrix dgCMatrix
#' @export
setMethod("convert_to_so", signature(sc_data = "dgCMatrix"), function(sc_data, ...) {
  Seurat::CreateSeuratObject(counts = sc_data, ...)
})
