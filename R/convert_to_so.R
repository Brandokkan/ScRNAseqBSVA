#' @describeIn convert_to_so Method for a \code{data.frame} of counts
#'   (genes x cells).
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod("convert_to_so", signature(sc_data = "data.frame"), function(sc_data,
                                                                       min.cells = 3,
                                                                       min.features = 200,
                                                                       ...) {
  CreateSeuratObject(counts = sc_data, min.cells, min.features, ...)
})

#' @describeIn convert_to_so Method for a file path to a supported
#'   single-cell data file (e.g. a .mtx or .csv file).
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod("convert_to_so", signature(sc_data = "character"), function(sc_data, min.cells = 3,
                                                                      min.features = 200,
                                                                      ...) {
  dots <- list(...)
  cso_fixed_names <- c("counts", "min.cells", "min.features", "...")
  cso_formals <- setdiff(names(formals(CreateSeuratObject)), cso_fixed_names)
  cso_args <- dots[names(dots) %in% cso_formals]
  read_args <- dots[!names(dots) %in% cso_formals]

  data <- do.call(read_sc_file, c(list(path = sc_data), read_args))

  cso_fixed_val <- list(counts = data, min.cells = min.cells, min.features = min.features)
  do.call(CreateSeuratObject, c(cso_fixed_val, cso_args))
})

#' @describeIn convert_to_so Method for a sparse counts matrix of class
#'   \code{dgCMatrix}.
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod("convert_to_so", signature(sc_data = "dgCMatrix"), function(sc_data,
                                                                      min.cells = 3,
                                                                      min.features = 200,
                                                                      ...) {
  CreateSeuratObject(counts = sc_data, min.cells, min.features, ...)
})
