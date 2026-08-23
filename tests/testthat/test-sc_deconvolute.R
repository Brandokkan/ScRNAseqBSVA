skip_if_not_installed("pbmc3k.SeuratData")
skip_if_not_installed("celldex")
skip_if_not_installed("scRNAseq")
skip_if_not_installed("scuttle")

data("so_mit_rbp")
so_complete <- sc_normalization(so_mit_rbp)
so_complete <- sc_clustering(so_complete)

e <- new.env()
data("pbmc3k", package = "pbmc3k.SeuratData", envir = e)
pbmc3k <- e$pbmc3k
pbmc3k <- Seurat::UpdateSeuratObject(pbmc3k)
pbmc3k <- sc_normalization(pbmc3k)
pbmc3k <- sc_clustering(pbmc3k)

suex_dataset <- celldex::MouseRNAseqData()

scex <- scRNAseq::ZeiselBrainData()
scex <- scuttle::logNormCounts(scex)
scex_not_norm <- scRNAseq::ZeiselBrainData()
scex_no_labels <- scex
scex_no_labels$level1class <- NULL

test_that("method for SummarizedExperiment and SingleCellExperiment works and gives
          sound results", {
    suex_labels <- SummarizedExperiment::colData(suex_dataset)$label.main
    suex_labels_p <- sc_deconvolute(suex_dataset, suex_dataset, diagnosis = FALSE)$predictions$labels
    concordance <- mean(suex_labels_p == suex_labels)
    expect_gt(concordance, 0.85)

    scex_labels <- scex$level1class
    scex_labels_p <- sc_deconvolute(scex, scex, diagnosis = FALSE, lab_name = "level1class")$predictions$labels
    concordance <- mean(scex_labels == scex_labels_p)
    expect_gt(concordance, 0.85)
})

test_that("method for SeuratObject works and gives sound results", {
    pbm_labels <- pbmc3k@meta.data$seurat_annotations
    pbm_labels_p <- sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE)@meta.data$predicted.id
    concordance <- mean(pbm_labels == pbm_labels_p, na.rm = TRUE)
    expect_gt(concordance, 0.85)
})

test_that("method for SeuratObject as query and SummarizedExperiment/SingleCellExperiment
          works", {
    suex_labels <- SummarizedExperiment::colData(suex_dataset)$label.main
    pbm_suex_labels_p <- sc_deconvolute(pbmc3k, suex_dataset,
        diagnosis = FALSE
    )@misc$predictions$labels
    expect_type(pbm_suex_labels_p, "character")
    expect_length(pbm_suex_labels_p, ncol(pbmc3k))

    scex_labels <- scex$level1class
    pbm_scex_labels_p <- sc_deconvolute(pbmc3k, scex,
        diagnosis = FALSE,
        lab_name = "level1class"
    )@misc$predictions$labels
    expect_type(pbm_scex_labels_p, "character")
    expect_length(pbm_scex_labels_p, ncol(pbmc3k))
})

test_that("function works with different database as query and reference", {
    expect_no_error(sc_deconvolute(so_complete, suex_dataset, diagnosis = FALSE))
})

test_that("not normalized query are rejected", {
    expect_error(sc_deconvolute(so_mit_rbp, pbmc3k, diagnosis = FALSE), "query was LogNormalized")
    expect_error(sc_deconvolute(so_mit_rbp, suex_dataset, diagnosis = FALSE), "query was LogNormalized")
    expect_error(sc_deconvolute(scex_not_norm, scex, diagnosis = FALSE), "query was LogNormalized")
})

test_that("reference datasets with no labels or with wrong lab_name asigned are rejected", {
    expect_error(sc_deconvolute(pbmc3k, so_mit_rbp, diagnosis = FALSE), "see what is the name of the vector")
    expect_error(sc_deconvolute(pbmc3k, scex_no_labels, diagnosis = FALSE), "see what is the name of the vector")
    expect_error(sc_deconvolute(scex, scex_no_labels, diagnosis = FALSE), "see what is the name of the vector")

    expect_error(
        sc_deconvolute(pbmc3k, pbmc3k, diagnosis = FALSE, lab_name = "wrong"),
        "see what is the name of the vector"
    )
    expect_error(
        sc_deconvolute(pbmc3k, suex_dataset, diagnosis = FALSE, lab_name = "wrong"),
        "see what is the name of the vector"
    )
    expect_error(
        sc_deconvolute(scex, scex, diagnosis = FALSE, lab_name = "wrong"),
        "see what is the name of the vector"
    )
})

test_that("sc_deconvolute plots a diagnostic graph when diagnosis = TRUE", {
    grDevices::pdf(NULL) # real device, writes nothing to disk
    grDevices::dev.control("enable") # make sure display-list recording is on
    on.exit(grDevices::dev.off(), add = TRUE)

    invisible(sc_deconvolute(suex_dataset, suex_dataset, diagnosis = TRUE))

    recorded <- grDevices::recordPlot()
    expect_gt(length(recorded[[1]]), 0) # display list is non-empty -> something was drawn
})

test_that("sc_deconvolute does not plot when diagnosis = FALSE", {
    grDevices::pdf(NULL)
    grDevices::dev.control("enable")
    on.exit(grDevices::dev.off(), add = TRUE)

    invisible(sc_deconvolute(suex_dataset, suex_dataset, diagnosis = FALSE))

    recorded <- grDevices::recordPlot()
    expect_length(recorded[[1]], 0) # nothing drawn
})
