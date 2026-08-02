#' Creates a Seurat object for downstream single-cell analysis from multiple file types
#'
#' @param sc_data Variable or path to file containing the single-cell sequencing data.
#' @param data_type The data or file type of \code{sc_data}.
#' @param project_name The name given to the Seurat object
#' @param min_cells Excludes genes not detected in at least this number of cells
#' @param min_features Excludes cells with less than this number of genes expression detected
#'
#' @importFrom Seurat CreateSeuratObject
#'
#' @return A Seurat object containing the single-cell data.
#' @export
#'
#' @examples
#' convert_to_so(sc_dataframe, "TestProject")
#' convert_to_so(sc_data.csv)
convert_to_so <- function(sc_data, data_type = "infer", project_name = "CreateSeuratObject",
                          min_cells = 3, min_features = 200){
  Seurat::CreateSeuratObject(counts = sc_data, project = project_name,
                             min.cells = min_cells, min.features = min_features)}
