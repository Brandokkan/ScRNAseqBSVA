mat_path <- system.file("extdata", "sc_test.mtx", package = "ScRNAseqBSVA")
gene_path <- system.file("extdata", "genes.tsv", package = "ScRNAseqBSVA")
cell_path <- system.file("extdata", "barcode.tsv", package = "ScRNAseqBSVA")
csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
tsv_path <- system.file("extdata", "sc_test.tsv", package = "ScRNAseqBSVA")
empty_file_path <- system.file("extdata", "empty_test", package = "ScRNAseqBSVA")

# test <- convert_to_so(mat_path, bar_path = cell_path, gene_path = gene_path)
# test_csv <- convert_to_so(csv_path)

test_that("files are loaded and trasnformed into SO correctlly", {
    expect_s4_class(
        convert_to_so(mat_path,
            bar_path = cell_path,
            gene_path = gene_path
        ),
        "Seurat"
    )
    expect_s4_class(convert_to_so(csv_path), "Seurat")
    expect_s4_class(convert_to_so(tsv_path), "Seurat")
})

test_that("files are converted correctlly", {
    data("so_obj")
    # so_obj_no_ens <- so_obj
    # rownames(so_obj_no_ens) <- make.unique(sub("-[^-]*$", "", rownames(so_obj_no_ens)))
    csv_df <- read.table(csv_path, header = TRUE, sep = ",", row.names = 1)
    tsv_df <- read.table(tsv_path, header = TRUE, sep = "\t", row.names = 1)
    csv_mat <- read_sc_file(csv_path)
    tsv_mat <- read_sc_file(tsv_path)
    mtx_mat <- read_sc_file(mat_path, bar_pat = cell_path, gene_path = gene_path)

    expect_equal(convert_to_so(csv_df, rem_ens_name = FALSE), so_obj)
    expect_equal(convert_to_so(tsv_df, rem_ens_name = FALSE), so_obj)
    expect_equal(convert_to_so(csv_mat, rem_ens_name = FALSE), so_obj)
    expect_equal(convert_to_so(tsv_mat, rem_ens_name = FALSE), so_obj)
    expect_equal(convert_to_so(mtx_mat, rem_ens_name = FALSE), so_obj)
})

test_that("unrecognized files and variables are not laoded", {
    expect_error(convert_to_so("test"))
    expect_error(convert_to_so(TRUE))
    expect_error(convert_to_so(5))
})
