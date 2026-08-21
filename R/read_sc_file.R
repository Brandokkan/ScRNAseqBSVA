#' @describeIn read_sc_file Method for creating the corresponding single-cell
#'  variable from a file or multiple files. For .mtx files, a gene and cell
#'  barcodes must also be submitted in .tsv format
#' @importFrom Seurat ReadMtx
#' @importFrom Matrix Matrix
#' @importFrom utils read.csv
#' @export
setMethod("read_sc_file", signature(path = "character"), function(path, bar_path, gene_path,
                                                                  ...) {
  if (grepl(".mtx$", path)) {
    if (is.character(bar_path) & is.character(gene_path)) {
      ReadMtx(path, bar_path, gene_path, feature.column = 1, ...)
    } else {
      stop("The suplied paths for the barcodes and/or gene names are not actually paths",
           call. = FALSE)
    }
  } else if (grepl(".csv$", path)) {
    df_file <- read.csv(path, row.names = 1, ...)
    Matrix(as.matrix(df_file), sparse = TRUE)
  } else if (grepl(".tsv$", path)) {
    df_file <- read.csv(path, sep = "\t", row.names = 1, ...)
    Matrix(as.matrix(df_file), sparse = TRUE)
  } else {
    stop("type of file not recognized or supported", call. = FALSE)
  }
})
