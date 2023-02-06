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
#' @export
#'
#' @examples
#' @noRd
plot_beta_diversity <- function(tse_obj, name, NMDS=FALSE, scale = F, ...) {
  if (NMDS == TRUE) {
    tse_obj <- mia::runNMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    xlab <- "Axis 1"
    ylab <- "Axis 2"
  } else {
    tse_obj <- mia::runMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    e <- attr(SingleCellExperiment::reducedDim(tse_obj, name), "eig")
    rel_eig <- 100 * e/sum(e[e>0])
    xlab <- paste("Axis 1 (", round(rel_eig[[1]], 2), "%)", sep = "")
    ylab <- paste("Axis 2 (", round(rel_eig[[2]], 2), "%)", sep = "")
  }
  p <- scater::plotReducedDim(tse_obj, name, colour_by = "condition") +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=colour_by), linetype=2) +
    ggplot2::ggtitle(name) +
    ggplot2::theme_bw()
  return(p)
}
