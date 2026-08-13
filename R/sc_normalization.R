#' @describeIn sc_normalization method to normalize the counts
#' @importFrom Seurat NormalizeData
#' @export
setMethod("sc_normalization", signature(so = "Seurat"), function(so, loga = TRUE,
                                                                 scale_factor = 10000, ...) {
  if (!is.numeric(scale_factor) | scale_factor <= 0) {
    stop("Error: scale_factor must be a positive number", call. = FALSE)
  }

  if (loga){
    NormalizeData(so, normalization.method = "LogNormalize", scale.factor = scale_factor, ...)
  } else {
    NormalizeData(so, normalization.method = "RC", scale.factor = scale_factor, ...)
  }
})
