#' @describeIn convert_to_so Method for a \code{data.frame} of counts
#'   (genes x cells).
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Matrix Matrix
#' @export
setMethod(
    "convert_to_so", signature(sc_data = "data.frame"),
    function(sc_data, min.cells = 3, min.features = 200, rem_ens_name = TRUE,
             ens_reg = "_[^_]*$", ...) {
        mat_data <- Matrix(as.matrix(sc_data), sparse = TRUE)
        if (rem_ens_name) {
            rownames(mat_data) <- make.unique(
                sub(ens_reg, "", rownames(mat_data))
            )
        }

        CreateSeuratObject(
            counts = mat_data, min.cells = min.cells,
            min.features = min.features, ...
        )
    }
)

#' @describeIn convert_to_so Method for a file path to a supported
#'   single-cell data file (e.g. a .mtx or .csv file).
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod(
    "convert_to_so", signature(sc_data = "character"),
    function(sc_data, min.cells = 3, min.features = 200, rem_ens_name = TRUE,
             ens_reg = "_[^_]*$", ...) {
        dots <- list(...)
        cso_fixed_names <- c("counts", "min.cells", "min.features", "...")
        cso_formals <- setdiff(
            names(formals(CreateSeuratObject)),
            cso_fixed_names
        )
        cso_args <- dots[names(dots) %in% cso_formals]
        read_args <- dots[!names(dots) %in% cso_formals]

        data.mtx <- do.call(read_sc_file, c(list(path = sc_data), read_args))

        if (rem_ens_name) {
            rownames(data.mtx) <- make.unique(
                sub(ens_reg, "", rownames(data.mtx))
            )
        }
        cso_fixed_val <- list(
            counts = data.mtx, min.cells = min.cells,
            min.features = min.features
        )
        do.call(CreateSeuratObject, c(cso_fixed_val, cso_args))
    }
)

#' @describeIn convert_to_so Method for a sparse counts matrix of
#'   class \code{dgCMatrix}.
#' @importClassesFrom Matrix dgCMatrix
#' @importFrom Seurat CreateSeuratObject
#' @export
setMethod(
    "convert_to_so", signature(sc_data = "dgCMatrix"),
    function(sc_data, min.cells = 3, min.features = 200, rem_ens_name = TRUE,
             ens_reg = "_[^_]*$", ...) {
        if (rem_ens_name) {
            rownames(sc_data) <- make.unique(
                sub(ens_reg, "", rownames(sc_data))
            )
        }
        CreateSeuratObject(
            counts = sc_data, min.cells = min.cells,
            min.features = min.features, ...
        )
    }
)
