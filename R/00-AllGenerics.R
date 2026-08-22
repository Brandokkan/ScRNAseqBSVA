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
#' @param plt_show string that tells the method What kind of plots to show.
#' @inheritParams calculate_mt_rbp
#' @param ... additional parameters
#'
#' @return plots for common quality control parameters.
#'
#' @export
#' @examples
#' data("so_mit_rbp")
#'
#' cell_qc_vis(so_mit_rbp, plt_show = "vln", specie = "mouse")
#'
#' cell_qc_vis(so_mit_rbp, plt_show = "scat", specie = "mouse")
setGeneric("cell_qc_vis", function(so, plt_show, specie = "human",
                                   mit_pat = "^MT-", rib_pat = "^RP[LS]",
                                   overwrite = FALSE,
                                   ...) {
  standardGeneric("cell_qc_vis")
})


#' Automatic filtering of cells based on the percentiles.
#'
#' This function automatically filters cell by using the percentiles of
#' the distributions of \code{nFeatures_RNA} and \code{percent.mt}.
#'
#' @param so the \code{SeuratObject} containing the cells to filter.
#' @param up_per the upper percentile of \code{nFeature_RNA}. Cells that are
#'  placed above this are discarded.
#' @param low_per the lower percentile of \code{nFeature_RNA}. Cells that are
#'  placed below this are discarded.
#' @param mit_per percentile used to filter cells based on mitochondrial genes.
#'  cells that have are placed below this are discarded.
#' @param ... additional parameters.
#'
#' @return the \code{SeuratObject} with the filtered cells removed.
#'
#' @export
#' @examples
#' data("so_mit_rbp")
#'
#' filtered_so <- automatic_filter(so_mit_rbp)
#'
#' ncol(so_mit_rbp)
#' ncol(filtered_so)
#'
setGeneric("automatic_filter", function(so, up_per = 0.9,
                                        low_per = 0.1,
                                        mit_per = 0.95,
                                        ...) {
  standardGeneric("automatic_filter")
})


#' Normalization of gene counts
#'
#' This function normalizes the counts inside a \code{SeuratObject} to
#' make comparison between different cell expression levels.
#' The default normalization is counts per ten thousand.
#'
#' @param so \code{SeuratObject} containing the counts to normalize
#' @param loga Boolean value telling the function if it should output the
#' logarithm of the normalized counts
#' @param scale_factor the scale factor that multiplies the division in order
#' to make it more human readable.
#' @param ... additional parameters
#'
#' @return the \code{SeuratObject} with the normalized counts
#'
#' @export
#' @examples
#' data("so_mit_rbp")
#'
#' norm_so <- sc_normalization(so_mit_rbp)
#' so_10k <- sc_normalization(so_mit_rbp, loga = FALSE)
#' so_1 <- sc_normalization(so_mit_rbp, loga = FALSE, scale_factor = 1)
#'
setGeneric("sc_normalization", function(so, loga = TRUE,
                                        scale_factor = 10000, ...) {
  standardGeneric("sc_normalization")
})


#' Function to find clusters using the most variable genes.
#'
#' This function first computes the most variable genes inside a \code{SeuratObject}
#' and then uses them for finding clusters using a graph in which the nodes (cells)
#' are connected to their K-nearest neighbors and. The edges are weighted using
#' Jaccadrd similarity.
#'
#' @param so \code{SeuratObject} containing the cells
#' @param n_var_genes the top number of variable genes used to compute the distances
#' for clustering
#' @param k_par number of nearest neighbors used for constructing the graph
#' @param res resolution. Higher values increase the number of clusters at the end.
#' it is suggested to keep this value between 0.4 and 1.2 for data set of around
#' 3K cells. The optimal number usually increases as the number of cells increases.
#'
#' @return the \code{SeuratObject} with the the clusters
#'
#' @export
#' @examples
#' data("so_mit_rbp")
#' norm_so <- sc_normalization(so_mit_rbp)
#'
#' sc_clustering(norm_so)
#'
setGeneric("sc_clustering", function(so, n_var_genes = 2000, k_par = 20,
                                     res = 0.5) {
  standardGeneric("sc_clustering")
})


#' Function to predict cell type using a deconvolution approach.
#'
#' This function use SingleR or Seurat, based on the type of reference database,
#' to predict the cell type of a query database using a reference database with
#' cell types already assigned. The reference database can be SeuratData,
#' SummarizedExperiment and SingleCellExperiment.
#'
#' @param so \code{SeuratObject} containing the cells
#' @param ref the reference database used to predict the cell types in so. It must
#' be one of: \code{SeuratObject}, \code{SummarizedExperiment} or \code{SingleCellExperiment}.
#' examples are \code{pbmc3k} from SeuratData or \code{MouseRNAseqData} from celldex.
#' @param diagnosis Boolean value to tell the function if it should also print
#' plots and information regarding predicted cell types confidence.
#' @param lab_name the name of the vector containing the names of the cell types
#' in the reference database.
#'
#' @return the \code{SeuratObject} with the predicted cell types in the meta.data
#' as "predictions".
#'
#' @export
#' @examples
#' if (requireNamespace("celldex", quietly = TRUE)) {
#'   suex_dataset <- celldex::MouseRNAseqData()
#'
#'   # SummarizedExperiment on SummarizedExperiment method
#'   sc_deconvolute(suex_dataset, suex_dataset)
#'
#'   if (requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
#'     e <- new.env()
#'     data("pbmc3k", package = "pbmc3k.SeuratData", envir = e)
#'     pbmc3k <- e$pbmc3k
#'     pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
#'     pbmc3k <- sc_normalization(pbmc3k)
#'     pbmc3k <- sc_clustering(pbmc3k)
#'
#'     # Seurat on SummarizedExperiment method
#'     sc_deconvolute(pbmc3k, suex_dataset)
#'
#'     # Seurat on Seurat method
#'     sc_deconvolute(pbmc3k, pbmc3k)
#'   }
#' }
setGeneric("sc_deconvolute", function(so, ref, diagnosis = TRUE, lab_name = "label.main") {
  standardGeneric("sc_deconvolute")
})


#' Function to plot information about the confidence of a cell type prediction
#'
#' This function produces plots and information from a \code{SeuratObject} that
#' was passed through sc_deconvolute to gain insight about the confidence of the
#' cell types predicted.
#'
#' @param so \code{SeuratObject} containig predicted cell types.
#'
#' @return plots and informations about the predicted cell types
#'
#' @export
#' @examples
#' if (requireNamespace("celldex", quietly = TRUE)) {
#'   suex_dataset <- celldex::MouseRNAseqData()
#'
#'   suex_pred <- sc_deconvolute(suex_dataset, suex_dataset, diagnosis = FALSE)
#'   prediction_diagnostics(suex_pred)
#' }
#'
setGeneric("prediction_diagnostics", function(so) {
  standardGeneric("prediction_diagnostics")
})


#' Function to compare the clusters to the predicted cell labels using
#' deconvolution
#'
#' This function is useful to compare the final results of the clusters obtained
#' from \code{sc_clustering} to the cell labels predicted by
#' \code{sc_deconvolute}.
#'
#' @param so \code{SeuratObject} containing the clusters found by
#' \code{sc_clustering} and the predicted cell labels found by
#' \code{sc_deconvolute}
#'
#' @return prints tables and a heat-map to compare the two results
#'
#' @export
#' @examples
#' if (requireNamespace("pbmc3k.SeuratData", quietly = TRUE)) {
#'   data("pbmc3k", package = "pbmc3k.SeuratData")
#'   pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
#'   pbmc3k <- sc_normalization(pbmc3k)
#'   pbmc3k <- sc_clustering(pbmc3k)
#'   pbmc3k <- sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE)
#'
#'   clu_dec_comparison(pbmc3k)
#' }
#'
setGeneric("clu_dec_comparison", function(so) {
  standardGeneric("clu_dec_comparison")
})
