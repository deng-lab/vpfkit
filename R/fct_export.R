#' Export a TSE object to RDS
#'
#' @param tse TreeSummarizedExperiment object
#' @param file Output file path
#' @return Invisible NULL
#' @export
export_vpftse <- function(tse, file) {
  saveRDS(tse, file)
  invisible(NULL)
}

#' Export abundance matrix from TSE
#'
#' @param tse TreeSummarizedExperiment object
#' @param file Output file path
#' @param assay.type Which assay to export (default: "counts")
#' @param format Output format: "csv" or "xlsx"
#' @return Invisible NULL
#' @export
#' @importFrom openxlsx write.xlsx
export_abundance <- function(tse, file, assay.type = "counts", format = c("csv", "xlsx")) {
  format <- match.arg(format)
  df <- SummarizedExperiment::assay(tse, assay.type) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Contig")
  if (format == "csv") {
    utils::write.csv(df, file, row.names = FALSE)
  } else {
    openxlsx::write.xlsx(df, file)
  }
  invisible(NULL)
}

#' Export contig annotations from TSE
#'
#' @param tse TreeSummarizedExperiment object
#' @param file Output file path
#' @param format Output format: "tsv" or "csv"
#' @return Invisible NULL
#' @export
export_annotations <- function(tse, file, format = c("tsv", "csv")) {
  format <- match.arg(format)
  df <- SummarizedExperiment::rowData(tse) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Contig")
  if (format == "tsv") {
    utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)
  } else {
    utils::write.csv(df, file, row.names = FALSE)
  }
  invisible(NULL)
}
