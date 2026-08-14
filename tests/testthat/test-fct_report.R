## Rendering a report costs several seconds, so the full renders below are kept
## to the cases that have actually broken: a destination outside the session
## temp directory, and objects with the shape a real ViroProfiler run produces.

make_report_tse <- function(n_contigs = 6, n_samples = 3,
                            assay_names = c("counts", "tpm", "covfrac"),
                            with_group = FALSE, with_genes = FALSE,
                            with_annotations = TRUE) {
  set.seed(7)
  contigs <- paste0("NODE_", seq_len(n_contigs), "_length_5000_cov_9.1")
  samples <- paste0("S", seq_len(n_samples))
  mk <- function(nm) {
    m <- matrix(
      if (nm == "covfrac") stats::runif(n_contigs * n_samples)
      else as.numeric(sample.int(500, n_contigs * n_samples, replace = TRUE)),
      nrow = n_contigs, dimnames = list(contigs, samples)
    )
    m
  }
  rd <- if (with_annotations) {
    S4Vectors::DataFrame(
      checkv_quality = rep(c("High-quality", "Low-quality"), length.out = n_contigs),
      checkv_contig_length = seq_len(n_contigs) * 2000L,
      checkv_completeness = seq(20, 95, length.out = n_contigs),
      Family = rep(c("Peduoviridae", NA), length.out = n_contigs),
      row.names = contigs
    )
  } else {
    S4Vectors::DataFrame(row.names = contigs)
  }
  cd <- S4Vectors::DataFrame(sample_name = samples, row.names = samples)
  if (with_group) cd$group <- rep(c("case", "control"), length.out = n_samples)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = stats::setNames(lapply(assay_names, mk), assay_names),
    rowData = rd, colData = cd
  )
  if (with_genes) {
    S4Vectors::metadata(se)$gene_annotations <- data.frame(
      gene = paste0("g", 1:4),
      dramv_category = c("AMG", "viral", "viral", NA),
      pharokka_card = c("NA", "aac(6')", NA, "NA"),
      stringsAsFactors = FALSE
    )
  }
  se
}

skip_without_quarto <- function() {
  skip_if_not(.vpf_quarto_available(), "Quarto CLI or the quarto R package is unavailable")
}

# --------------------------------------------------------------------------
# Output location
# --------------------------------------------------------------------------

test_that("generate_report writes to the requested path outside the temp directory", {
  skip_without_quarto()
  tse <- make_report_tse()
  ## The destination must not be a plain tempfile(): the previous
  ## implementation only appeared to work because Quarto wrote its output into
  ## the session temp directory, where a tempfile() destination happened to
  ## point as well.
  root <- withr::local_tempdir()
  out <- file.path(root, "reports", "deeply", "nested", "quality_report.html")

  result <- generate_report(tse, out)

  expect_equal(result, out)
  expect_true(file.exists(out))
  expect_gt(file.size(out), 10000)
  head_lines <- readLines(out, n = 5, warn = FALSE)
  expect_true(any(grepl("<!DOCTYPE|<html", head_lines, ignore.case = TRUE)))
})

test_that("generate_report leaves no intermediate files beside the output", {
  skip_without_quarto()
  tse <- make_report_tse()
  root <- withr::local_tempdir()
  out <- file.path(root, "report.html")
  generate_report(tse, out)
  expect_equal(list.files(root, recursive = TRUE), "report.html")
})

test_that("the rendered report embeds its figures instead of linking them", {
  skip_without_quarto()
  tse <- make_report_tse()
  out <- file.path(withr::local_tempdir(), "report.html")
  generate_report(tse, out)

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  n_img <- lengths(regmatches(html, gregexpr("<img", html, fixed = TRUE)))
  expect_gt(n_img, 0)
  ## Every image must be a data URI; a self-contained report cannot reference
  ## files from a staging directory that no longer exists.
  n_data <- lengths(regmatches(html, gregexpr("<img[^>]*src=\"data:image", html)))
  expect_equal(n_img, n_data)
})

# --------------------------------------------------------------------------
# Graceful degradation
# --------------------------------------------------------------------------

test_that("a report renders for the shape a real ViroProfiler run produces", {
  skip_without_quarto()
  ## Two samples, no grouping variable, no gene annotations, and the legacy
  ## trimmed-mean assay name.
  tse <- make_report_tse(n_samples = 2, assay_names = c("counts", "tmm", "covfrac"))
  out <- file.path(withr::local_tempdir(), "report.html")
  expect_no_error(generate_report(tse, out))

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "at least three")
  expect_match(html, "No grouping variable")
  expect_match(html, "No gene annotations")
})

test_that("degenerate objects still produce a report", {
  skip_without_quarto()
  root <- withr::local_tempdir()

  cases <- list(
    one_sample   = make_report_tse(n_samples = 1),
    one_contig   = make_report_tse(n_contigs = 1),
    no_rowdata   = make_report_tse(with_annotations = FALSE),
    no_covfrac   = make_report_tse(assay_names = "counts"),
    covfrac_only = make_report_tse(assay_names = "covfrac")
  )
  for (nm in names(cases)) {
    out <- file.path(root, paste0(nm, ".html"))
    ## covfrac_only legitimately warns that it has no abundance assay.
    expect_no_error(suppressWarnings(generate_report(cases[[nm]], out)))
    expect_true(file.exists(out), info = nm)
  }
})

test_that("the report refuses to compute diversity from a covered fraction", {
  skip_without_quarto()
  tse <- make_report_tse(assay_names = c("covfrac", "counts"))
  out <- file.path(withr::local_tempdir(), "report.html")
  generate_report(tse, out, assay.type = "covfrac")

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "breadth of coverage rather than abundance")
})

test_that("the report echoes the coverage-breadth threshold it was given", {
  skip_without_quarto()
  tse <- make_report_tse()
  root <- withr::local_tempdir()

  masked <- file.path(root, "masked.html")
  generate_report(tse, masked, covfrac_threshold = 0.75)
  expect_match(paste(readLines(masked, warn = FALSE), collapse = "\n"),
               "covered fraction &gt;= 0.75|covered fraction >= 0.75")

  unmasked <- file.path(root, "unmasked.html")
  generate_report(tse, unmasked)
  expect_match(paste(readLines(unmasked, warn = FALSE), collapse = "\n"),
               "not applied")
})

test_that("gene annotations are summarized when present", {
  skip_without_quarto()
  tse <- make_report_tse(with_genes = TRUE, with_group = TRUE)
  out <- file.path(withr::local_tempdir(), "report.html")
  generate_report(tse, out)

  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "Genes annotated")
  expect_match(html, "AMG genes")
  expect_match(html, "Candidate grouping variable")
})

# --------------------------------------------------------------------------
# Input validation
# --------------------------------------------------------------------------

test_that("generate_report validates its arguments before rendering", {
  tse <- make_report_tse()
  out <- file.path(withr::local_tempdir(), "report.html")

  expect_error(generate_report(data.frame(x = 1), out), "must be a SummarizedExperiment")
  skip_without_quarto()
  expect_error(generate_report(tse, out, covfrac_threshold = 1.5), "\\[0, 1\\]")
  expect_error(generate_report(tse, out, covfrac_threshold = "a"), "\\[0, 1\\]")
  expect_error(generate_report(tse, out, covfrac_threshold = c(0.1, 0.2)), "\\[0, 1\\]")
  expect_error(generate_report(tse, ""), "non-empty file path")
  expect_error(generate_report(tse, dirname(out)), "existing directory")
})

test_that("generate_report warns about a mistyped assay but still renders", {
  skip_without_quarto()
  tse <- make_report_tse()
  out <- file.path(withr::local_tempdir(), "report.html")
  expect_warning(generate_report(tse, out, assay.type = "nope"), "Assay 'nope' not found")
  expect_true(file.exists(out))
})

test_that("generate_report explains how to install Quarto when it is missing", {
  ## The availability check is a separate function precisely so that this path
  ## can be tested without uninstalling anything.
  local_mocked_bindings(.vpf_quarto_available = function() FALSE)
  tse <- make_report_tse()
  out <- file.path(withr::local_tempdir(), "report.html")

  expect_error(generate_report(tse, out), "quarto\\.org|install\\.packages")
})

test_that("the Quarto availability check agrees with the CLI on this machine", {
  available <- .vpf_quarto_available()
  expect_type(available, "logical")
  if (requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))) {
    expect_true(available)
  }
})
