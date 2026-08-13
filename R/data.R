#' Example single-cell object
#'
#' A subsetted Seurat object used in examples and tests.
#'
#' @format A `Seurat` object with 78 cells and 9337 features.
#' @source Derived from Panglao dataset.
"so_obj"

#' Example single-cell object without the ENSEMBL suffix
#'
#' A subsetted Seurat object used in examples and tests. The ENSEBLE suffix was
#' removed.
#'
#' @format A `Seurat` object with 78 cells and 9337 features.
#' @source Derived from Panglao dataset.
"so_no_ens"

#' Example single-cell object with percentage of mitochondrial and ribosomial
#' genes calculated.
#'
#' A subsetted Seurat object used in examples and tests. The percentage of
#' ribosomal and mitochondrial genes are already calculated for each cell.
#'
#' @format A `Seurat` object with 78 cells and 9337 features.
#' @source Derived from Panglao dataset.
"so_mit_rbp"

#' Example single-cell object with an additional normalized count layer
#'
#' This Seurat object contains an additional layer called "data" that contains
#' the normalized counts per ten thousand.
#'
#' @format A `Seurat` object with 78 cells and 9337 features.
#' @source Derived from Panglao dataset.
"so_norm"

#' Example single-cell object with an additional log-normalized count layer
#'
#' This Seurat object contains an additional layer called "data" that contains
#' the log-normalized counts per ten thousand.
#'
#' @format A `Seurat` object with 78 cells and 9337 features.
#' @source Derived from Panglao dataset.
"so_norm_log"
