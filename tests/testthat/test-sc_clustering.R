data("so_mit_rbp")
norm_so <- sc_normalization(so_mit_rbp)

test_that("the number of final variable genes are the same as the set one", {
    n_fea <- 100
    clu_so <- sc_clustering(norm_so, n_var_genes = n_fea)
    expect_equal(length(Seurat::VariableFeatures(clu_so)), n_fea)
})

test_that("all data used to for the clustering is present in the Surat Object", {
    clu_so <- sc_clustering(norm_so)
    expect_no_error(clu_so@reductions)
    expect_no_error(clu_so@graphs)
})

test_that("an error occurs with invalid values", {
    expect_error(sc_clustering("hello"))
    expect_error(sc_clustering(clu_so, n_var_genes = -20))
    expect_error(sc_clustering(clu_so, k_par = -9))
    expect_error(sc_clustering(clu_so, res = -18))
})
