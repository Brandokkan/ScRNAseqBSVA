## code to prepare `SeuratObject` dataset
mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")

so_obj <- Seurat::ReadMtx(mat_path, cell_path, gene_path, feature.column = 1)

usethis::use_data(so_obj, overwrite = TRUE, compress = "xz")
