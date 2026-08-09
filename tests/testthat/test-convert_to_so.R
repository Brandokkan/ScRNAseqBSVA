mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
tsv_path <- system.file("extdata", "sc_test.tsv", package = "ScRNAseqBSVA")
empty_file_path <- system.file("extdata", "empty_test", package = "ScRNAseqBSVA")

# test <- convert_to_so(mat_path, bar_path = cell_path, gene_path = gene_path)
# test_csv <- convert_to_so(csv_path)

test_that("files are loaded and trasnformed into SO correctlly", {
  expect_s4_class(convert_to_so(mat_path, bar_path = cell_path,
                                gene_path = gene_path),
                  "Seurat")
  expect_s4_class(convert_to_so(csv_path), "Seurat")
  expect_s4_class(convert_to_so(tsv_path), "Seurat")
})
