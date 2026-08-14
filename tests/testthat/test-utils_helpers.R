make_helper_tse <- function(n_contigs = 8, n_samples = 6, group = TRUE,
                            assay_names = c("counts", "covfrac")) {
  set.seed(11)
  contigs <- paste0("c", seq_len(n_contigs))
  samples <- paste0("S", seq_len(n_samples))
  mk <- function(nm) {
    matrix(
      if (nm == "covfrac") stats::runif(n_contigs * n_samples)
      else as.numeric(sample.int(200, n_contigs * n_samples, replace = TRUE)),
      nrow = n_contigs, dimnames = list(contigs, samples)
    )
  }
  cd <- S4Vectors::DataFrame(sample_name = samples, row.names = samples)
  if (group) cd$group <- rep(c("case", "control"), length.out = n_samples)
  ## A TreeSummarizedExperiment, not a bare SummarizedExperiment: mia's
  ## ordination writes into reducedDims, which only the richer class has.
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = stats::setNames(lapply(assay_names, mk), assay_names),
    colData = cd
  )
}

# --------------------------------------------------------------------------
# Assay name resolution
# --------------------------------------------------------------------------

test_that("legacy and current assay names resolve to each other", {
  expect_equal(vpf_canonical_assay("tmm"), "trimmed_mean")
  expect_equal(vpf_canonical_assay("trimmed_mean"), "trimmed_mean")
  expect_equal(vpf_canonical_assay("covered_fraction"), "covfrac")
  expect_equal(vpf_canonical_assay("something_else"), "something_else")

  expect_equal(vpf_match_assay(c("counts", "tmm"), "trimmed_mean"), "tmm")
  expect_equal(vpf_match_assay(c("counts", "trimmed_mean"), "tmm"), "trimmed_mean")
  expect_equal(vpf_match_assay(c("counts", "covfrac"), "covered_fraction"), "covfrac")
  expect_equal(vpf_match_assay(c("counts"), "counts"), "counts")
})

test_that("an unresolvable assay name errors or returns NA on request", {
  expect_error(vpf_match_assay(c("counts", "tpm"), "nope"), "Assay 'nope' not found")
  expect_true(is.na(vpf_match_assay(c("counts", "tpm"), "nope", required = FALSE)))
  expect_error(vpf_match_assay(c("counts"), NA_character_), "single non-empty assay name")
  expect_error(vpf_match_assay(c("counts"), c("a", "b")), "single non-empty assay name")
})

test_that("assay descriptions distinguish trimmed mean from edgeR's TMM", {
  expect_match(vpf_assay_description("tmm"), "coverage depth")
  expect_match(vpf_assay_description("trimmed_mean"), "coverage depth")
  expect_match(vpf_assay_description("covfrac"), "breadth")
  expect_match(vpf_assay_description("counts"), "read count")
  expect_match(vpf_assay_description("mystery"), "unknown")
})

# --------------------------------------------------------------------------
# Coverage-breadth masking
# --------------------------------------------------------------------------

test_that("refind_abundance zeroes only what falls below the threshold", {
  tse <- make_helper_tse(n_contigs = 3, n_samples = 2, group = FALSE)
  SummarizedExperiment::assay(tse, "counts") <-
    matrix(c(10, 20, 30, 1, 2, 3), nrow = 3,
           dimnames = dimnames(SummarizedExperiment::assay(tse, "counts")))
  SummarizedExperiment::assay(tse, "covfrac") <-
    matrix(c(0.9, 0.4, 0.5, 0.2, 1.0, 0.49), nrow = 3,
           dimnames = dimnames(SummarizedExperiment::assay(tse, "covfrac")))

  out <- SummarizedExperiment::assay(refind_abundance(tse, "counts", "covfrac", 0.5), "counts")
  ## A contig sitting exactly at the threshold is kept.
  expect_equal(as.vector(out), c(10, 0, 30, 0, 2, 0))
})

test_that("refind_abundance treats a missing covered fraction as absent", {
  tse <- make_helper_tse(n_contigs = 2, n_samples = 1, group = FALSE)
  SummarizedExperiment::assay(tse, "counts") <- matrix(c(10, 20), nrow = 2,
    dimnames = dimnames(SummarizedExperiment::assay(tse, "counts")))
  SummarizedExperiment::assay(tse, "covfrac") <- matrix(c(NA_real_, 0.9), nrow = 2,
    dimnames = dimnames(SummarizedExperiment::assay(tse, "covfrac")))

  out <- SummarizedExperiment::assay(refind_abundance(tse, "counts", "covfrac", 0.5), "counts")
  ## NA must not propagate into the abundance matrix.
  expect_false(anyNA(out))
  expect_equal(as.vector(out), c(0, 20))
})

test_that("refind_abundance rejects a threshold outside [0, 1]", {
  tse <- make_helper_tse()
  expect_error(refind_abundance(tse, "counts", "covfrac", 1.5), "between 0 and 1")
  expect_error(refind_abundance(tse, "counts", "covfrac", -0.1), "between 0 and 1")
  expect_error(refind_abundance(tse, "counts", "covfrac", NA), "non-missing number")
  expect_error(refind_abundance(tse, "counts", "covfrac", "0.5"), "non-missing number")
  expect_error(refind_abundance(tse, "counts", "covfrac", c(0.1, 0.2)), "non-missing number")
})

test_that("refind_abundance is a no-op at threshold 0 and keeps only full coverage at 1", {
  tse <- make_helper_tse(n_contigs = 3, n_samples = 1, group = FALSE)
  SummarizedExperiment::assay(tse, "counts") <- matrix(c(5, 6, 7), nrow = 3,
    dimnames = dimnames(SummarizedExperiment::assay(tse, "counts")))
  SummarizedExperiment::assay(tse, "covfrac") <- matrix(c(0, 0.5, 1), nrow = 3,
    dimnames = dimnames(SummarizedExperiment::assay(tse, "covfrac")))

  expect_equal(as.vector(SummarizedExperiment::assay(
    refind_abundance(tse, "counts", "covfrac", 0), "counts")), c(5, 6, 7))
  expect_equal(as.vector(SummarizedExperiment::assay(
    refind_abundance(tse, "counts", "covfrac", 1), "counts")), c(0, 0, 7))
})

test_that("refind_abundance warns when asked to mask an already-normalized assay", {
  tse <- make_helper_tse(assay_names = c("tpm", "covfrac"))
  expect_warning(refind_abundance(tse, "tpm", "covfrac", 0.5), "already normalized")
})

test_that("refind_abundance accepts legacy assay names and a missing covfrac assay", {
  tse <- make_helper_tse(assay_names = c("counts", "tmm"))
  expect_warning(out <- refind_abundance(tse, "trimmed_mean", "covfrac", 0.5),
                 "not found")
  expect_equal(SummarizedExperiment::assay(out, "tmm"),
               SummarizedExperiment::assay(tse, "tmm"))
})

# --------------------------------------------------------------------------
# Output path preparation
# --------------------------------------------------------------------------

test_that("vpf_prepare_outfile creates parents and returns an absolute path", {
  root <- withr::local_tempdir()
  p <- vpf_prepare_outfile(file.path(root, "x", "y", "out.csv"))
  expect_true(dir.exists(file.path(root, "x", "y")))
  expect_true(startsWith(p, "/") || grepl("^[A-Za-z]:", p))
})

test_that("vpf_prepare_outfile rejects unusable destinations", {
  root <- withr::local_tempdir()
  expect_error(vpf_prepare_outfile(root), "existing directory")
  expect_error(vpf_prepare_outfile(""), "non-empty file path")
  expect_error(vpf_prepare_outfile(NA_character_), "non-empty file path")
})

# --------------------------------------------------------------------------
# Grouping-variable detection
# --------------------------------------------------------------------------

test_that("a grouping column is found when one exists", {
  expect_equal(.vpf_guess_group_column(make_helper_tse(group = TRUE)), "group")
})

test_that("no grouping column is reported when every column is an identifier", {
  ## This is the shape of a real ViroProfiler object: colData holds only
  ## sample_name, which is unique per sample and therefore not a grouping
  ## variable.
  expect_null(.vpf_guess_group_column(make_helper_tse(group = FALSE)))
})

test_that("a constant column is not treated as a grouping variable", {
  tse <- make_helper_tse(group = FALSE)
  SummarizedExperiment::colData(tse)$batch <- rep("b1", ncol(tse))
  expect_null(.vpf_guess_group_column(tse))
})

# --------------------------------------------------------------------------
# Beta diversity
# --------------------------------------------------------------------------

test_that("plot_beta_diversity returns NULL rather than ordinating fewer than three samples", {
  expect_null(plot_beta_diversity(make_helper_tse(n_samples = 2), "bd",
                                  NMDS = TRUE, method = "bray", assay.type = "counts"))
  expect_null(plot_beta_diversity(make_helper_tse(n_samples = 1), "bd",
                                  NMDS = TRUE, method = "bray", assay.type = "counts"))
})

test_that("plot_beta_diversity colours by a detected grouping column", {
  skip_if_not_installed("scater")
  skip_if_not_installed("vegan")
  tse <- make_helper_tse(group = TRUE)
  p <- suppressWarnings(suppressMessages(
    plot_beta_diversity(tse, "bd", NMDS = TRUE, method = "bray", assay.type = "counts")))
  expect_s3_class(p, "ggplot")
  ## Building the plot is what catches an aes() that refers to a column
  ## scater no longer produces.
  expect_no_error(suppressWarnings(ggplot2::ggplot_build(p)))
})

test_that("plot_beta_diversity warns about an unknown colour column but still plots", {
  skip_if_not_installed("scater")
  skip_if_not_installed("vegan")
  tse <- make_helper_tse(group = TRUE)
  expect_warning(
    p <- suppressMessages(plot_beta_diversity(tse, "bd", NMDS = TRUE, method = "bray",
                                              assay.type = "counts", colour_by = "nope")),
    "not in colData"
  )
  expect_s3_class(p, "ggplot")
})
