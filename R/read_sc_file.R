#' @describeIn read_sc_file Method for creating the corresponding single-cell
#'  variable from a file or multiple files. For .mtx files, a gene and cell
#'  barcodes must also be submitted in .tsv format
#' @param bar_path the path of the .tsv file containing the cells bar-codes
#' @param gene_path the path of the .tsv file containing the names of the genes
#' @importFrom Seurat ReadMtx
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
    read.csv(path, row.names = 1, ...)
  } else if (grepl(".tsv$", path)) {
    read.csv(path, sep = "\t", ...)
  } else {
    stop("type of file not recognized or supported", call. = FALSE)
  }
})
