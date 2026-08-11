csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")

so_no_ens <- convert_to_so(csv_path)

usethis::use_data(so_no_ens, overwrite = TRUE, compress = "xz")
