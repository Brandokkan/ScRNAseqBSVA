
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ScRNAseqBSVA

<!-- badges: start -->

[![R-CMD-check](https://github.com/Brandokkan/ScRNAseqBSVA/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Brandokkan/ScRNAseqBSVA/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

ScRNAseqBSVA streamlines single-cell RNA-seq analysis: it loads the
data, preprocesses it, clusters the cells, assigns cell type labels by
deconvolution, and evaluates how well those two agree.

The package is built around [Seurat](https://satijalab.org/seurat/), and
collapses the most commonly used Seurat steps into single functions with
sensible defaults, so a standard workflow can be run end to end with a
handful of calls. It also interfaces with the Bioconductor world:
references from `SummarizedExperiment` / `SingleCellExperiment` objects
(e.g. `celldex` datasets) are accepted directly, and a `Seurat` object
can be handed to `Seurat::as.SingleCellExperiment()` at any point.

## Installation

The package depends on Bioconductor packages, so install `BiocManager`
first:

``` r
# install.packages("BiocManager")
BiocManager::install(c("SingleR", "SummarizedExperiment"))
```

You can then install the development version of ScRNAseqBSVA from
[GitHub](https://github.com/Brandokkan/ScRNAseqBSVA) with:

``` r
# install.packages("pak")
pak::pak("Brandokkan/ScRNAseqBSVA")
```

Some of the suggested packages used in the examples and in the vignette
(`SeuratData`, `pbmc3k.SeuratData`, `celldex`, `scRNAseq`, `scuttle`)
are not on CRAN and are only needed to reproduce them:

``` r
pak::pak("satijalab/seurat-data")
SeuratData::InstallData("pbmc3k")
BiocManager::install(c("celldex", "scRNAseq", "scuttle"))
```

R (\>= 4.6.0) is required.

## Overview

| Function | What it does |
|----|----|
| `read_sc_file()` | Reads a `.mtx`, `.csv` or `.tsv` count file into a `dgCMatrix` |
| `convert_to_so()` | Builds a `Seurat` object from a path, `data.frame` or `dgCMatrix`, filtering low-quality genes/cells and stripping the ENSEMBL suffix |
| `calculate_mt_rbp()` | Adds the per-cell percentage of mitochondrial and ribosomal protein genes to the metadata |
| `cell_qc_vis()` | Violin (`"vln"`) or scatter (`"scat"`) plots of the usual QC metrics |
| `automatic_filter()` | Percentile-based filtering on `nFeature_RNA` and `percent.mt` |
| `sc_normalization()` | Counts-per-10k normalization, log-transformed by default |
| `sc_clustering()` | Variable features, scaling, PCA, KNN graph and Louvain clustering in one call |
| `sc_deconvolute()` | Predicts cell types against a `Seurat`, `SummarizedExperiment` or `SingleCellExperiment` reference (Seurat label transfer or SingleR, depending on the reference) |
| `prediction_diagnostics()` | Score heatmap and summary of low-confidence calls for the predicted labels |
| `clu_dec_comparison()` | Compares the clusters against the predicted labels (tables + heatmap) |

## Example

A full workflow, from a count file to the comparison between clusters
and predicted cell types:

``` r
library(ScRNAseqBSVA)

# 1. load the counts and build a Seurat object
csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")
so <- convert_to_so(csv_path)

# 2. quality control
so <- calculate_mt_rbp(so, specie = "mouse")
cell_qc_vis(so, "vln", specie = "mouse")
cell_qc_vis(so, "scat", specie = "mouse")

# 3. filter, normalize and cluster
so <- automatic_filter(so)
so <- sc_normalization(so)
so <- sc_clustering(so)

# 4. label the cells against a reference, then inspect the calls
#    (the reference can be a Seurat, SummarizedExperiment or
#     SingleCellExperiment object)
ref <- celldex::MouseRNAseqData()
so <- sc_deconvolute(so, ref, diagnosis = FALSE)
prediction_diagnostics(so)

# 5. how well do the clusters and the predicted labels agree?
clu_dec_comparison(so)
```

The vignette walks through the same steps on the `pbmc3k` dataset, with
the plots and a discussion of the output:

``` r
vignette("working-with-ScRNAseqBSVA", package = "ScRNAseqBSVA")
```

## Example data

Small `Seurat` objects (78 cells, 9337 features, derived from the
Panglao dataset) are shipped with the package for examples and tests,
one per stage of the workflow: `so_obj`, `so_no_ens`, `so_mit_rbp`,
`so_norm` and `so_norm_log`. Raw count files in the supported formats
live in `inst/extdata`.

``` r
data("so_mit_rbp")
```

## Getting help

Every function has a help page (`?sc_deconvolute`); please report bugs
and feature requests at
<https://github.com/Brandokkan/ScRNAseqBSVA/issues>.

## License

AGPL (\>= 3). See [LICENSE.md](LICENSE.md).
