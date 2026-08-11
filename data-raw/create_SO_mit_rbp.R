csv_path <- system.file("extdata", "sc_test.csv", package = "ScRNAseqBSVA")

so_mit_rbp <- convert_to_so(csv_path)
so_mit_rbp <- calculate_mt_rbp(so_mit_rbp, specie = "mouse")

usethis::use_data(so_mit_rbp, overwrite = TRUE, compress = "xz")
