#' Creates a Seurat object for single-cell analysis from a file or
#' various variable types.
#'
#' Creates a \code{Seurat} object for downstream single-cell analysis,
#' dispatching on the class of \code{sc_data}.
#'
#' @param sc_data Object or path containing the single-cell data. See the Methods
#'   section below for supported input classes.
#' @param min.cells A gene must be expressed in at least this number of cells to
#'  be included. Default is \code{min.cells = 3}.
#' @param rem_ens_name if \code{TRUE}, then an ENSEMBL suffix is assumed to be
#'  present in the gene names and is removed.
#' @param ens_reg what regular expression is used to recognize and substitute
#'  the gene name from the ENSEMBL suffix.
#'  Ignored if \code{rem_ens_name} is \code{FALSE}.
#' @param min.features A cell must express at least this number of genes to be
#'  included. Default is \code{min.features = 200}.
#' @param ... Additional arguments passed to
#'   \code{\link[Seurat]{CreateSeuratObject}}, e.g. \code{project}
#'
#' @return A \code{Seurat} object.
#'
#' @seealso \code{\link[Seurat]{CreateSeuratObject}}
#'
#' @export
setGeneric("convert_to_so", function(sc_data, min.cells = 3,
                                     min.features = 200,
                                     rem_ens_name = TRUE,
                                     ens_reg = "_[^_]*$",
                                     ...) {
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
