#' Pick a colData column usable as a grouping variable
#'
#' A grouping variable is a categorical column with at least two levels that is
#' not simply a per-sample identifier. The previous implementation looked for a
#' column literally named `condition`; no ViroProfiler object has one, so it
#' silently produced an uncoloured ordination on every real dataset.
#'
#' @param tse_obj A `SummarizedExperiment`.
#' @return The column name, or `NULL` when there is no usable one.
#' @noRd
.vpf_guess_group_column <- function(tse_obj) {
  cd <- SummarizedExperiment::colData(tse_obj)
  if (ncol(cd) == 0L) return(NULL)
  n <- nrow(cd)
  for (nm in colnames(cd)) {
    v <- cd[[nm]]
    if (!(is.character(v) || is.factor(v) || is.logical(v))) next
    k <- length(unique(v[!is.na(v)]))
    ## Two or more levels, but not one level per sample: a column of unique
    ## sample names carries no group structure.
    if (k >= 2L && k < n) return(nm)
  }
  NULL
}

#' Plot beta diversity using the \code{mia} package
#'
#' @description Ordinates samples and returns a `ggplot2` object. Returns
#'   `NULL`, rather than failing, when the data cannot support an ordination.
#'
#' @param tse_obj a \code{TreeSummarizedExperiment} object
#' @param name name of the reduced dimension
#' @param NMDS if TRUE, run NMDS instead of MDS
#' @param colour_by name of a `colData` column used to colour the points.
#'   `NULL` (default) picks the first column that looks like a grouping
#'   variable; `NA` disables colouring.
#' @param ... additional arguments to \code{vegdist} or \code{cmdscale}
#'
#' @return a \code{ggplot2} object, or `NULL`
#' @importFrom SingleCellExperiment reducedDim
#' @noRd
plot_beta_diversity <- function(tse_obj, name, NMDS = FALSE, colour_by = NULL, ...) {
  ## Two samples give exactly one distance, which has no ordination geometry:
  ## every ordination of two points is the same line segment.
  if (ncol(tse_obj) < 3L) return(NULL)

  if (is.null(colour_by)) {
    colour_by <- .vpf_guess_group_column(tse_obj)
  } else if (length(colour_by) == 1L && is.na(colour_by)) {
    colour_by <- NULL
  } else if (!colour_by %in% colnames(SummarizedExperiment::colData(tse_obj))) {
    warning(sprintf("Column '%s' is not in colData; points will not be coloured.", colour_by),
            call. = FALSE)
    colour_by <- NULL
  }

  if (NMDS) {
    tse_obj <- mia::runNMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    xlab <- "Axis 1"
    ylab <- "Axis 2"
  } else {
    tse_obj <- mia::addMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    e <- attr(SingleCellExperiment::reducedDim(tse_obj, name), "eig")
    if (is.null(e) || !any(e > 0)) {
      xlab <- "Axis 1"
      ylab <- "Axis 2"
    } else {
      rel_eig <- 100 * e / sum(e[e > 0])
      xlab <- paste0("Axis 1 (", round(rel_eig[[1]], 2), "%)")
      ylab <- paste0("Axis 2 (", round(rel_eig[[2]], 2), "%)")
    }
  }

  p <- scater::plotReducedDim(tse_obj, name, colour_by = colour_by) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(name) +
    ggplot2::theme_bw()

  ## ggplot2's default t-distribution ellipse needs more than three points per
  ## group; with fewer it drops the group and warns every time the plot is
  ## drawn, which in a Shiny app means on every reactive update.
  if (!is.null(colour_by)) {
    grp <- SummarizedExperiment::colData(tse_obj)[[colour_by]]
    if (length(grp) && min(table(grp[!is.na(grp)])) >= 4L) {
      p <- p + ggplot2::stat_ellipse(
        ggplot2::aes(fill = .data$colour_by),
        geom = "polygon", alpha = 0.1, linetype = 2
      )
    }
  }
  p
}
