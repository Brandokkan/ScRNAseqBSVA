data("so_mit_rbp")

so_norm_log <- Seurat::NormalizeData(so_mit_rbp, normalization.method = "LogNormalize", scale.factor = 10000)
so_norm <- Seurat::NormalizeData(so_mit_rbp, normalization.method = "RC", scale.factor = 10000)

usethis::use_data(so_norm, overwrite = TRUE, compress = "xz")
usethis::use_data(so_norm_log, overwrite = TRUE, compress = "xz")
