mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
tsv_path <- system.file("extdata", "sc_test.tsv", package = "ScRNAseqBSVA")
empty_file_path <- system.file("extdata", "empty_test", package = "ScRNAseqBSVA")

invalid_file_str <- "type of file not recognized or supported"
invalid_paths_str <- "The suplied paths for the barcodes and/or gene names are not actually paths"

test_that("mtx files are read as a dgCmMatrix", {
  expect_s4_class(read_sc_file(mat_path, cell_path, gene_path), "dgCMatrix")
})

test_that("csv files are read correctlly", {
  expect_type(read_sc_file(csv_path), "list")
})

test_that("cell bar-codes and gene names are loaded correctlly into the
          dgCMatrix", {
  dgc_mat <- read_sc_file(mat_path, cell_path, gene_path)
  gene_name_mat <- as.vector(dgc_mat@Dimnames[[1]])
  bar_name_mat <- as.vector(dgc_mat@Dimnames[[2]])
  gene_name_file <- read.csv(gene_path, header = FALSE, sep = "\t")$"V1"
  bar_name_file <- read.csv(cell_path, header = FALSE, sep = "\t")$"V1"
  expect_true(all.equal.character(bar_name_mat, bar_name_file))
  expect_true(all.equal.character(gene_name_mat, gene_name_file))
})

test_that("an unrecognized file is not loaded", {
  expect_error(read_sc_file(empty_file_path), invalid_file_str)
})

test_that("the suplied paths for gene names and barcodes are not valid", {
  expect_error(read_sc_file(mat_path, 1, 2), invalid_paths_str)
})
