#' @describeIn sc_clustering method for computing the most variable
#' genes and use them for clustering
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat ScaleData
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat VariableFeatures
#' @importFrom utils head
#' @export
setMethod(
    "sc_clustering", signature(so = "Seurat"),
    function(so, n_var_genes = 2000, k_par = 20, res = 0.5) {
        so <- FindVariableFeatures(so,
            selection.method = "vst", nfeatures = n_var_genes
        )

        so <- ScaleData(so)

        so <- RunPCA(so, features = VariableFeatures(object = so))

        pc.sevfiv <- (so$pca@stdev)^2
        pc.sevfiv <- pc.sevfiv / sum(pc.sevfiv)
        pc.sevfiv <- head(cumsum(pc.sevfiv), 50)
        pc.sevfiv <- min(which(pc.sevfiv >= 0.75))

        so <- FindNeighbors(so, dims = seq_len(pc.sevfiv), k.param = k_par)

        so <- FindClusters(so, resolution = res)
        so
    }
)
