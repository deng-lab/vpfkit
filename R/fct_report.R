## Quarto report generation.

#' Is a usable Quarto installation available?
#'
#' Split out from [generate_report()] so that the failure path can be tested
#' without uninstalling anything.
#'
#' @return `TRUE` when both the `quarto` R package and the Quarto CLI are
#'   available.
#' @noRd
.vpf_quarto_available <- function() {
  if (!requireNamespace("quarto", quietly = TRUE)) return(FALSE)
  path <- tryCatch(quarto::quarto_path(), error = function(e) NULL)
  !is.null(path) && nzchar(path) && file.exists(path)
}

#' Message explaining how to make Quarto usable
#'
#' @return A single string.
#' @noRd
.vpf_quarto_hint <- function() {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    return(paste0(
      "Package 'quarto' is required for report generation. ",
      "Install it with install.packages('quarto'), and install the Quarto CLI ",
      "from https://quarto.org/docs/get-started/."
    ))
  }
  paste0(
    "The Quarto CLI was not found. Install it from ",
    "https://quarto.org/docs/get-started/, or point R at an existing ",
    "installation with Sys.setenv(QUARTO_PATH = '/path/to/quarto'). ",
    "RStudio ships its own copy, which R only sees when RStudio's ",
    "'quarto' is on PATH."
  )
}

#' Generate a ViroProfiler quality report
#'
#' Renders a self-contained HTML report summarizing the contents of a TSE
#' object: dataset composition, CheckV quality, taxonomy, coverage breadth,
#' diversity, abundance and gene annotations. Every section degrades to an
#' explanatory note when the data it needs is absent, so a minimal object
#' still produces a complete report rather than a render failure.
#'
#' The report is rendered in a private staging directory and only then copied
#' to `output_file`, so no intermediate file is ever written next to the
#' caller's data and two concurrent calls cannot overwrite each other.
#'
#' @param tse A `SummarizedExperiment` (usually a `TreeSummarizedExperiment`).
#' @param output_file Output HTML file path. The parent directory is created
#'   when missing. Relative paths are resolved against the current working
#'   directory before rendering starts.
#' @param covfrac_threshold Optional minimum covered fraction used to mask
#'   abundances before the report computes anything from them. `NULL` (the
#'   default) leaves the assays exactly as stored and the report says so.
#'   When set, masking is applied to raw counts *before* normalization and the
#'   threshold is echoed in the report.
#' @param assay.type Assay used for the abundance and diversity sections.
#'   Legacy assay names are resolved automatically.
#' @param title Report title.
#' @param quiet Suppress Quarto's progress output (default `TRUE`).
#' @return The path of the report that was written, invisibly.
#' @export
#' @examples
#' \dontrun{
#' generate_report(tse, "reports/viroprofiler_report.html")
#' generate_report(tse, "reports/masked.html", covfrac_threshold = 0.75)
#' }
generate_report <- function(tse, output_file, covfrac_threshold = NULL,
                            assay.type = "counts",
                            title = "ViroProfiler Quality Report",
                            quiet = TRUE) {
  vpf_assert_se(tse)
  if (!.vpf_quarto_available()) {
    stop(.vpf_quarto_hint(), call. = FALSE)
  }
  if (!is.null(covfrac_threshold)) {
    if (!is.numeric(covfrac_threshold) || length(covfrac_threshold) != 1L ||
        is.na(covfrac_threshold) || covfrac_threshold < 0 || covfrac_threshold > 1) {
      stop("`covfrac_threshold` must be NULL or a single number in [0, 1].", call. = FALSE)
    }
  }

  ## Objects written before the rename store the trimmed mean as `tmm`; the
  ## report should name the quantity the same way regardless.
  tse <- normalize_assay_names(tse, quiet = TRUE)

  ## Resolve the destination before Quarto runs, so that a bad path fails now
  ## rather than after a render that then has nowhere to go.
  output_file <- vpf_prepare_outfile(output_file)

  ## A mistyped assay name would otherwise only surface as a missing section
  ## halfway down the finished report. Warn rather than stop: an object that
  ## genuinely lacks an abundance assay should still get its other sections.
  if (is.na(vpf_match_assay(tse, assay.type, required = FALSE))) {
    warning(
      sprintf(
        "Assay '%s' not found; abundance sections will be omitted. Available assays: %s.",
        assay.type,
        if (length(SummarizedExperiment::assayNames(tse)))
          paste(sprintf("'%s'", SummarizedExperiment::assayNames(tse)), collapse = ", ")
        else "none"
      ),
      call. = FALSE
    )
  }

  template <- system.file("report_template.qmd", package = "vpfkit")
  if (!nzchar(template)) {
    stop("Report template not found in the installed vpfkit package.", call. = FALSE)
  }

  ## A private staging directory per call. Quarto writes its output next to
  ## the input document, so a shared directory would let concurrent renders
  ## collide on the same output name.
  stage <- tempfile("vpfkit-report-")
  if (!dir.create(stage, recursive = TRUE)) {
    stop("Could not create a staging directory for the report.", call. = FALSE)
  }
  on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)
  stage <- normalizePath(stage, winslash = "/", mustWork = TRUE)

  staged_qmd <- file.path(stage, "report.qmd")
  staged_rds <- file.path(stage, "tse.rds")
  staged_html <- file.path(stage, "report.html")

  if (!isTRUE(file.copy(template, staged_qmd))) {
    stop("Could not copy the report template into the staging directory.", call. = FALSE)
  }
  saveRDS(tse, staged_rds)

  ## Absolute paths only: chunks execute in the staging directory, so a
  ## relative path would be resolved against the wrong root.
  params <- list(
    tse_path = normalizePath(staged_rds, winslash = "/", mustWork = TRUE),
    assay_type = assay.type,
    title = title,
    covfrac_threshold = if (is.null(covfrac_threshold)) "none" else format(covfrac_threshold)
  )

  quarto::quarto_render(
    input = staged_qmd,
    output_format = "html",
    output_file = "report.html",
    execute_params = params,
    execute_dir = stage,
    quiet = quiet,
    as_job = FALSE
  )

  ## quarto_render() signals most failures, but a render that produces no file
  ## while still returning normally would otherwise be reported as success.
  if (!file.exists(staged_html)) {
    stop(
      sprintf("Quarto returned without error but produced no report at '%s'.", staged_html),
      call. = FALSE
    )
  }
  if (file.size(staged_html) == 0L) {
    stop("Quarto produced an empty report.", call. = FALSE)
  }
  if (!isTRUE(file.copy(staged_html, output_file, overwrite = TRUE)) ||
      !file.exists(output_file)) {
    stop(sprintf("The report was rendered but could not be copied to '%s'.", output_file),
         call. = FALSE)
  }

  invisible(output_file)
}
