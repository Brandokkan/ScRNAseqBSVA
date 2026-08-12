#' Creates a Seurat object for single-cell analysis from a file or
#' various variable types.
#'
#' Creates a \code{Seurat} object for downstream single-cell analysis,
#' dispatching on the class of \code{sc_data}. By default, it filters
#' cells and genes by min.cells and min.features and removes the ensemble
#' suffix.
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
#' @examples
#' mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
#' gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
#' cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
#' mtx_mat <- convert_to_so(mat_path, gene_path = gene_path, bar_path = cell_path)
#' mtx_mat
#'
#' csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
#' csv_mat <- convert_to_so(csv_path)
#' csv_mat
#'
#' tsv_path <- system.file("extdata", "sc_test.tsv", package = "ScRNAseqBSVA")
#' tsv_df <- read.table(tsv_path, header = TRUE, sep = "\t", row.names = 1)
#' tsv_mat <- convert_to_so(tsv_df)
#' tsv_mat
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
#' @param bar_path the path of the .tsv file containing the cells bar-codes
#' @param gene_path the path of the .tsv file containing the names of the genes
#' @param ... additional parameters
#'
#' @return a variable containing the information in the file(s)
#'  in the form of a \code{dgCMatrix}.
#'
#' @export
#' @examples
#' mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
#' gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
#' cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
#' mtx_mat <- read_sc_file(mat_path, gene_path = gene_path, bar_path = cell_path)
#' mtx_mat
#'
#' csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
#' csv_mat <- read_sc_file(csv_path)
#' csv_mat
setGeneric("read_sc_file", function(path, ...) {
  standardGeneric("read_sc_file")
})


#' Adds to a \code{SeuratObject} metadata the percentage of mitochondrial and
#' ribosomal protein genes per cell.
#'
#' This function calculates and adds to the \code{SeuratObject} metadata, if
#' not already present, the percentage of mitochondrial and ribosomal protein genes.
#'
#' @param so the \code{SeuratObject} that contains the data to visualize
#' @param specie the specie from which the data was derived. Allows to know
#'  the correct nomenclature of the gene names. it can be set as: human, mouse
#'  or other. default value is human.
#' @param mit_pat reg-ex pattern for finding the mitochondrial genes.
#'  Only change if \code{specie} is set to "other".
#' @param rib_pat reg-ex pattern for finding the ribosomal protein genes.
#'  Only change if \code{specie} is set to "other".
#' @param overwrite logical value to tell the function if it should overwrite
#' the current values of mitochondrial and ribosomal genes. If set to \code{TRUE}
#' while the values are not yet calculated, it raises an error.
#'
#' @return the \code{SeuratObject} with the new metadata
#'
#' @export
#' @examples
#' csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
#'
#' so_mit_rbp <- convert_to_so(csv_path)
#' so_mit_rbp <- calculate_mt_rbp(so_mit_rbp, specie = "mouse")
#'
setGeneric("calculate_mt_rbp", function(so, specie = "human",
                                        mit_pat = "^MT-", rib_pat = "^RP[LS]",
                                        overwrite = FALSE) {
  standardGeneric("calculate_mt_rbp")
})


#' Visualize, and calculate if needed, common metrics for cell quality control.
#'
#' This function streamlines the pipeline used to visualize cell quality control
#' in Seurat for easier use. It returns violin plots for the number of feature,
#' the number of counts, the percentage of mitochondrial DNA and the percentage
#' of DNA that codes for ribosomal proteins per cell. The last two parameters are
#' calculated if not already present in the \code{SeuratObject}.
#'
#' @param plt_show What kind of plots to show. R markdown is suggested for showing
#'  both the scatter and violoin plots.
#' @inheritParams calculate_mt_rbp
#' @param ... additional parameters
#'
#' @return plots for common quality control parameters.
#'
#' @export
setGeneric("cell_qc_vis", function(so, specie = "human",
                                   mit_pat = "^MT-", rib_pat = "^RP[LS]",
                                   overwrite = FALSE,
                                   plt_show = c("scat", "vln"),
                                   ...) {
  standardGeneric("cell_qc_vis")
})
