#' Generate a ViroProfiler quality report
#'
#' Renders a Quarto HTML report summarizing the contents of a TSE object,
#' including quality metrics, taxonomy, diversity, and gene annotations.
#'
#' @param tse TreeSummarizedExperiment object
#' @param output_file Output HTML file path
#' @param quiet Suppress Quarto output (default: TRUE)
#' @return Invisible path to the generated report
#' @export
generate_report <- function(tse, output_file, quiet = TRUE) {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    stop("Package 'quarto' is required for report generation. Install with: install.packages('quarto')", call. = FALSE)
  }

  # Save TSE to temp file for the template to read
  tse_tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tse_tmp), add = TRUE)
  saveRDS(tse, tse_tmp)

  # Find template
  template <- system.file("report_template.qmd", package = "vpfkit")
  if (template == "") stop("Report template not found in vpfkit package", call. = FALSE)

  # Copy template to a temp location (quarto renders in-place)
  tmp_qmd <- tempfile(fileext = ".qmd")
  on.exit(unlink(tmp_qmd), add = TRUE)
  file.copy(template, tmp_qmd)

  # Render
  quarto::quarto_render(
    input = tmp_qmd,
    output_file = basename(output_file),
    execute_params = list(tse_path = tse_tmp),
    quiet = quiet
  )

  # Move rendered file to desired location
  rendered <- sub("\\.qmd$", ".html", tmp_qmd)
  if (file.exists(rendered)) {
    file.copy(rendered, output_file, overwrite = TRUE)
    unlink(rendered)
  }

  invisible(output_file)
}
