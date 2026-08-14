## Export helpers.
##
## All three exporters share the same contract: validate first, write into a
## directory that is known to exist and be writable, and return the path that
## was actually written so a caller can act on it. None of them may report
## success for a file that is not on disk.

#' Flatten columns a delimited text file cannot represent
#'
#' `rowData()` may carry list, `List` or `DataFrame` columns. `write.table()`
#' fails on those with "unimplemented type 'list' in 'EncodeElement'", which
#' does not say which column caused it.
#'
#' @param df A data frame.
#' @param collapse Separator used to join multi-valued entries.
#' @return `df` with every column atomic.
#' @noRd
.vpf_flatten_columns <- function(df, collapse = ";") {
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.atomic(col) && is.null(dim(col))) next
    df[[nm]] <- vapply(
      seq_len(nrow(df)),
      function(i) {
        v <- tryCatch(as.character(unlist(col[i], use.names = FALSE)), error = function(e) NA_character_)
        if (length(v) == 0L) NA_character_ else paste(v, collapse = collapse)
      },
      character(1L)
    )
  }
  df
}

#' Turn a matrix-like assay into a data frame with an explicit ID column
#'
#' @param mat Matrix or matrix-like object.
#' @param id_column Name of the identifier column to prepend.
#' @return A data frame whose first column holds the row identifiers.
#' @noRd
.vpf_matrix_to_df <- function(mat, id_column = "Contig") {
  mat <- as.matrix(mat)
  ids <- rownames(mat)
  if (is.null(ids)) {
    warning(
      "Matrix has no rownames; the identifier column will hold row numbers.",
      call. = FALSE
    )
    ids <- as.character(seq_len(nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("V", seq_len(ncol(mat)))
  }
  df <- as.data.frame(mat, stringsAsFactors = FALSE, optional = TRUE)
  rownames(df) <- NULL
  out <- cbind(stats::setNames(data.frame(ids, stringsAsFactors = FALSE), id_column), df)
  out
}

#' Make a string usable as an Excel worksheet name
#'
#' Excel rejects worksheet names longer than 31 characters, names containing
#' any of `[ ] : * ? / \`, empty names, and names wrapped in apostrophes.
#' `openxlsx` passes such names straight through and the resulting workbook
#' will not open.
#'
#' @param name Desired sheet name.
#' @param fallback Name used when `name` sanitizes to nothing.
#' @return A valid worksheet name.
#' @noRd
.vpf_safe_sheet_name <- function(name, fallback = "Sheet1") {
  if (!is.character(name) || length(name) != 1L || is.na(name)) name <- fallback
  clean <- gsub("[\\[\\]:*?/\\\\]", "_", name, perl = TRUE)
  clean <- gsub("^'+|'+$", "", clean)
  clean <- trimws(clean)
  if (!nzchar(clean)) clean <- fallback
  if (nchar(clean) > 31L) clean <- substr(clean, 1L, 31L)
  ## "History" is reserved by Excel.
  if (identical(tolower(clean), "history")) clean <- substr(paste0(clean, "_"), 1L, 31L)
  clean
}

#' Write a data frame as delimited text
#'
#' Tab-separated output is written unquoted, because that is what the tools
#' downstream of a ViroProfiler run expect. Unquoted output is only safe if no
#' field contains the separator, so embedded tabs, carriage returns and
#' newlines are replaced by a single space first and the affected columns are
#' named in a warning. Without that step one tab inside a free-text annotation
#' shifts every later column of that row and nothing reports an error.
#'
#' Comma-separated output goes through `utils::write.csv()`, which quotes
#' properly and needs no sanitizing.
#'
#' @param df Data frame to write.
#' @param file Destination path.
#' @param sep Field separator.
#' @return `file`, invisibly.
#' @noRd
.vpf_write_delim <- function(df, file, sep) {
  if (identical(sep, "\t")) {
    df <- .vpf_sanitize_delims(df)
    utils::write.table(df, file, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE, na = "NA")
  } else {
    utils::write.csv(df, file, row.names = FALSE)
  }
  if (!file.exists(file)) {
    stop(sprintf("Failed to write '%s'.", file), call. = FALSE)
  }
  invisible(file)
}

#' Replace characters that would break an unquoted delimited file
#'
#' @param df A data frame.
#' @return `df` with tabs, carriage returns and newlines replaced by spaces in
#'   every character or factor column.
#' @noRd
.vpf_sanitize_delims <- function(df) {
  affected <- character(0)
  for (nm in names(df)) {
    col <- df[[nm]]
    if (!is.character(col) && !is.factor(col)) next
    chr <- as.character(col)
    if (!any(grepl("[\t\r\n]", chr), na.rm = TRUE)) next
    affected <- c(affected, nm)
    df[[nm]] <- gsub("[\t\r\n]+", " ", chr)
  }
  if (length(affected)) {
    warning(
      sprintf(
        "Replaced embedded tabs/newlines with spaces in column(s): %s.",
        paste(affected, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  df
}

#' Export a TSE object to RDS
#'
#' @param tse A `SummarizedExperiment` (usually a `TreeSummarizedExperiment`).
#' @param file Output file path. The parent directory is created when missing.
#' @return The path written, invisibly.
#' @export
#' @examples
#' \dontrun{
#' export_vpftse(tse, "results/viroprofiler_output.rds")
#' }
export_vpftse <- function(tse, file) {
  vpf_assert_se(tse)
  file <- vpf_prepare_outfile(file)
  saveRDS(tse, file)
  if (!file.exists(file)) {
    stop(sprintf("Failed to write '%s'.", file), call. = FALSE)
  }
  invisible(file)
}

#' Export an abundance matrix from a TSE object
#'
#' Writes one row per contig and one column per sample, with the contig
#' identifiers in a leading `Contig` column so that they survive formats that
#' have no concept of row names.
#'
#' Legacy assay names are resolved automatically: asking for `trimmed_mean`
#' finds an assay stored under the older name `tmm`, and vice versa. Note that
#' `tmm` here has always meant CoverM's trimmed mean of per-base coverage
#' depth, not edgeR's Trimmed Mean of M-values.
#'
#' @param tse A `SummarizedExperiment`.
#' @param file Output file path. The parent directory is created when missing.
#' @param assay.type Assay to export (default `"counts"`).
#' @param format Output format, `"csv"` or `"xlsx"`.
#' @param sheet_name Worksheet name for `format = "xlsx"`. Defaults to the
#'   assay name. Sanitized to Excel's rules (31 characters, no `[ ] : * ? / \`).
#' @return The path written, invisibly.
#' @export
#' @importFrom openxlsx write.xlsx
#' @examples
#' \dontrun{
#' export_abundance(tse, "abundance_counts.csv")
#' export_abundance(tse, "abundance.xlsx", assay.type = "tpm", format = "xlsx")
#' }
export_abundance <- function(tse, file, assay.type = "counts",
                             format = c("csv", "xlsx"), sheet_name = NULL) {
  vpf_assert_se(tse)
  format <- match.arg(format)
  ## Objects written before the rename store the trimmed mean as `tmm`.
  ## Normalizing here means the exported worksheet is named for the quantity
  ## it holds, whichever spelling the object arrived with.
  tse <- normalize_assay_names(tse, quiet = TRUE)
  stored <- vpf_match_assay(tse, assay.type)
  file <- vpf_prepare_outfile(file)

  df <- .vpf_matrix_to_df(SummarizedExperiment::assay(tse, stored), id_column = "Contig")

  if (format == "csv") {
    return(.vpf_write_delim(df, file, sep = ","))
  }

  ## Excel's hard limits. Hitting them produces a workbook that will not open,
  ## which is worse than refusing to write one.
  if (nrow(df) + 1L > 1048576L) {
    stop(sprintf(
      "Cannot write %d rows to xlsx (Excel's limit is 1,048,576). Use format = \"csv\".",
      nrow(df)
    ), call. = FALSE)
  }
  if (ncol(df) > 16384L) {
    stop(sprintf(
      "Cannot write %d columns to xlsx (Excel's limit is 16,384). Use format = \"csv\".",
      ncol(df)
    ), call. = FALSE)
  }

  sheet <- .vpf_safe_sheet_name(if (is.null(sheet_name)) stored else sheet_name)
  openxlsx::write.xlsx(stats::setNames(list(df), sheet), file = file)
  if (!file.exists(file)) {
    stop(sprintf("Failed to write '%s'.", file), call. = FALSE)
  }
  invisible(file)
}

#' Export contig annotations from a TSE object
#'
#' Writes `rowData()` with the contig identifiers in a leading `Contig` column.
#' Multi-valued columns (`List`, list, nested `DataFrame`) are collapsed with
#' `;` so that they survive a flat text format.
#'
#' @param tse A `SummarizedExperiment`.
#' @param file Output file path. The parent directory is created when missing.
#' @param format Output format, `"tsv"` or `"csv"`.
#' @return The path written, invisibly.
#' @export
#' @examples
#' \dontrun{
#' export_annotations(tse, "contig_annotations.tsv")
#' }
export_annotations <- function(tse, file, format = c("tsv", "csv")) {
  vpf_assert_se(tse)
  format <- match.arg(format)
  file <- vpf_prepare_outfile(file)

  rd <- as.data.frame(SummarizedExperiment::rowData(tse), optional = TRUE)
  rd <- .vpf_flatten_columns(rd)
  ids <- rownames(SummarizedExperiment::rowData(tse))
  if (is.null(ids)) {
    ids <- rownames(tse)
  }
  if (is.null(ids)) {
    warning("TSE has no rownames; the Contig column will hold row numbers.", call. = FALSE)
    ids <- as.character(seq_len(nrow(rd)))
  }
  rownames(rd) <- NULL
  df <- cbind(data.frame(Contig = ids, stringsAsFactors = FALSE), rd)

  .vpf_write_delim(df, file, sep = if (format == "tsv") "\t" else ",")
}
