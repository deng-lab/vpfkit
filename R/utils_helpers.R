#' Plot beta diversity using \code{mia} package

#' @description A utils function
#'
#' @param tse_obj a \code{TreeSummarizedExperiment} object
#' @param name name of the reduced dimension
#' @param NMDS if TRUE, run NMDS instead of MDS
#' @param scale if TRUE, scale the data
#' @param ... additional arguments to \code{vegdist} or \code{cmdscale}
#'
#' @return a \code{ggplot2} object
#' @importFrom SingleCellExperiment reducedDim
#' @noRd
plot_beta_diversity <- function(tse_obj, name, NMDS=FALSE, scale = F, ...) {
  if (ncol(tse_obj) < 2) return(NULL)
  if (NMDS == TRUE) {
    tse_obj <- mia::runNMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    xlab <- "Axis 1"
    ylab <- "Axis 2"
  } else {
    tse_obj <- mia::addMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    e <- attr(SingleCellExperiment::reducedDim(tse_obj, name), "eig")
    rel_eig <- 100 * e/sum(e[e>0])
    xlab <- paste("Axis 1 (", round(rel_eig[[1]], 2), "%)", sep = "")
    ylab <- paste("Axis 2 (", round(rel_eig[[2]], 2), "%)", sep = "")
  }
  cond <- if ("condition" %in% names(SummarizedExperiment::colData(tse_obj))) "condition" else NULL
  p <- scater::plotReducedDim(tse_obj, name, colour_by = cond) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(name) +
    ggplot2::theme_bw()
  if (!is.null(cond)) {
    p <- p + ggplot2::stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=colour_by), linetype=2)
  }
  return(p)
}
