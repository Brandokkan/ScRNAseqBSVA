#' Package-wide imports.
#'
#' \code{methods} provides the S4 machinery (\code{setGeneric} and
#' \code{setMethod}) used throughout the package. \code{standardGeneric} needs
#' no import: it is a primitive in \code{base}.
#'
#' \code{pheatmap} and \code{viridis} are never called here, but
#' \code{SingleR::plotScoreHeatmap} calls \code{pheatmap::pheatmap} and
#' \code{viridis::viridis} internally and only lists them under \code{Suggests},
#' so they are not guaranteed to be installed alongside SingleR. Declaring them
#' in \code{Imports} makes them hard requirements of this package, which is what
#' \code{\link{prediction_diagnostics}} needs to plot at all.
#'
#' @importFrom methods setGeneric setMethod
#' @importFrom pheatmap pheatmap
#' @importFrom viridis viridis
#' @noRd
NULL


# Metadata columns that Seurat's subset() resolves against the object rather
# than the calling frame. R CMD check's static analysis cannot see that, so it
# reports them as undefined globals; this declares them as known names.
utils::globalVariables(c("nFeature_RNA", "percent.mt"))


#' Draws the SingleR score heatmap on the device the caller was using.
#'
#' \code{plotScoreHeatmap} relies on \code{pheatmap}, which opens a temporary
#' \code{pdf(file = NULL)} device to measure the label sizes and then closes it
#' with \code{dev.off()}. \code{dev.off()} makes the \emph{next} device in the
#' list current rather than the one that was current before, so when more than
#' one device is open the heatmap is silently drawn on the wrong one. Building
#' the heatmap with \code{silent = TRUE} and restoring the original device
#' before drawing it keeps the plot where the caller expects it.
#'
#' The gtable is drawn explicitly rather than through \code{print()}, which
#' Bioconductor discourages outside \code{show} methods. This is exactly what
#' \code{pheatmap:::print.pheatmap} does, so the output is unchanged.
#'
#' @param predictions the prediction \code{DFrame} returned by SingleR
#'
#' @return the \code{pheatmap} object, invisibly
#'
#' @importFrom SingleR plotScoreHeatmap
#' @importFrom grDevices dev.cur dev.set
#' @importFrom grid grid.draw
#' @noRd
plot_prediction_scores <- function(predictions) {
    dev_before <- dev.cur()

    scores <- plotScoreHeatmap(predictions, silent = TRUE)

    # pheatmap's temporary device may have left another device current
    if (dev_before != 1L && dev.cur() != dev_before) {
        dev.set(dev_before)
    }

    grid.draw(scores$gtable)
    invisible(scores)
}
