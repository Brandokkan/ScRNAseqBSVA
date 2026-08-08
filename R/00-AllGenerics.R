#' Creates a Seurat object for single-cell analysis from a file or
#' various variable types.
#'
#' Creates a \code{Seurat} object for downstream single-cell analysis,
#' dispatching on the class of \code{sc_data}.
#'
#' @param sc_data Object or path containing the single-cell data. See the Methods
#'   section below for supported input classes.
#' @param ... Additional arguments passed to
#'   \code{\link[Seurat]{CreateSeuratObject}}, e.g. \code{project},
#'   \code{min.cells}, \code{min.features}.
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}
#'
#' @export
setGeneric("convert_to_so", function(sc_data, ...) {
  standardGeneric("convert_to_so")
})


#' Reads one or multiple single cell files and outputs
#' the corresponding variable containing the information
#' inside the file(s)
#'
#' The output variable can be used for single-cell analysis
#'
#' @param path the path of the file containing the read counts
#' @param ... additional parameters
#'
#' @return a variable containing the information in the file(s)
#'  in the corresponding format recognized formats
#'
#' @export
#' @examples
#' mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
#' gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
#' cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
#' read_sc_file(mat_path, gene_path, cell_path)
#'
#' csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
#' read_sc_file(csv_path, header = FALSE)
setGeneric("read_sc_file", function(path, ...) {
  standardGeneric("read_sc_file")
})
