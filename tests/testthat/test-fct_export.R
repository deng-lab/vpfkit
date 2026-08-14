## Test objects are built directly rather than through create_vpftse(), so that
## these tests exercise the exporters and nothing else.

make_tse <- function(n_contigs = 5, n_samples = 3, assay_names = c("counts", "tpm", "covfrac"),
                     rowdata = TRUE, rownames = TRUE) {
  set.seed(42)
  contigs <- paste0("NODE_", seq_len(n_contigs), "_length_1000_cov_5.5")
  samples <- paste0("S", seq_len(n_samples))
  mk <- function(nm) {
    m <- matrix(
      if (nm == "covfrac") stats::runif(n_contigs * n_samples)
      else as.numeric(sample.int(1000, n_contigs * n_samples, replace = TRUE)),
      nrow = n_contigs, ncol = n_samples
    )
    if (rownames) rownames(m) <- contigs
    colnames(m) <- samples
    m
  }
  assays <- stats::setNames(lapply(assay_names, mk), assay_names)
  rd <- if (rowdata) {
    S4Vectors::DataFrame(
      checkv_quality = rep(c("High-quality", "Low-quality"), length.out = n_contigs),
      checkv_contig_length = seq_len(n_contigs) * 1000L,
      checkv_completeness = seq(10, 90, length.out = n_contigs),
      Family = rep(c("Peduoviridae", NA), length.out = n_contigs),
      row.names = if (rownames) contigs else NULL
    )
  } else {
    S4Vectors::DataFrame(row.names = if (rownames) contigs else NULL)
  }
  SummarizedExperiment::SummarizedExperiment(
    assays = assays, rowData = rd,
    colData = S4Vectors::DataFrame(sample_name = samples, row.names = samples)
  )
}

# --------------------------------------------------------------------------
# export_vpftse
# --------------------------------------------------------------------------

test_that("export_vpftse round-trips a TSE", {
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".rds")
  expect_equal(export_vpftse(tse, tmp), tmp)
  expect_true(file.exists(tmp))
  loaded <- readRDS(tmp)
  expect_s4_class(loaded, "SummarizedExperiment")
  expect_equal(dim(loaded), dim(tse))
  expect_equal(SummarizedExperiment::assay(loaded, "counts"),
               SummarizedExperiment::assay(tse, "counts"))
  expect_equal(as.data.frame(SummarizedExperiment::rowData(loaded)),
               as.data.frame(SummarizedExperiment::rowData(tse)))
})

test_that("export_vpftse rejects objects that are not SummarizedExperiments", {
  tmp <- withr::local_tempfile(fileext = ".rds")
  expect_error(export_vpftse(data.frame(a = 1), tmp), "must be a SummarizedExperiment")
})

# --------------------------------------------------------------------------
# export_abundance
# --------------------------------------------------------------------------

test_that("export_abundance CSV round-trips values and contig names", {
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".csv")
  export_abundance(tse, tmp, assay.type = "counts", format = "csv")

  back <- utils::read.csv(tmp, check.names = FALSE)
  orig <- SummarizedExperiment::assay(tse, "counts")
  expect_equal(colnames(back), c("Contig", colnames(orig)))
  expect_equal(back$Contig, rownames(orig))
  expect_equal(as.matrix(back[, -1]), unname(orig), ignore_attr = TRUE)
})

test_that("export_abundance XLSX round-trips and names the sheet after the assay", {
  skip_if_not_installed("openxlsx")
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  export_abundance(tse, tmp, assay.type = "tpm", format = "xlsx")

  expect_equal(openxlsx::getSheetNames(tmp), "tpm")
  back <- openxlsx::read.xlsx(tmp)
  orig <- SummarizedExperiment::assay(tse, "tpm")
  expect_equal(nrow(back), nrow(orig))
  expect_equal(back$Contig, rownames(orig))
  expect_equal(as.matrix(back[, -1]), unname(orig), ignore_attr = TRUE)
})

test_that("export_abundance sanitizes worksheet names to Excel's rules", {
  skip_if_not_installed("openxlsx")
  tse <- make_tse()

  tmp <- withr::local_tempfile(fileext = ".xlsx")
  export_abundance(tse, tmp, format = "xlsx", sheet_name = strrep("a", 60))
  expect_equal(nchar(openxlsx::getSheetNames(tmp)), 31L)

  tmp2 <- withr::local_tempfile(fileext = ".xlsx")
  export_abundance(tse, tmp2, format = "xlsx", sheet_name = "bad[]:*?/\\name")
  sheet <- openxlsx::getSheetNames(tmp2)
  expect_false(grepl("[\\[\\]:*?/\\\\]", sheet, perl = TRUE))

  tmp3 <- withr::local_tempfile(fileext = ".xlsx")
  export_abundance(tse, tmp3, format = "xlsx", sheet_name = "   ")
  expect_true(nzchar(openxlsx::getSheetNames(tmp3)))
})

test_that("export_abundance resolves legacy and current assay names in both directions", {
  legacy <- make_tse(assay_names = c("counts", "tmm"))
  current <- make_tse(assay_names = c("counts", "trimmed_mean"))

  f1 <- withr::local_tempfile(fileext = ".csv")
  expect_silent(export_abundance(legacy, f1, assay.type = "trimmed_mean"))
  f2 <- withr::local_tempfile(fileext = ".csv")
  expect_silent(export_abundance(current, f2, assay.type = "tmm"))

  expect_equal(utils::read.csv(f1)$Contig, rownames(legacy))
  expect_equal(utils::read.csv(f2)$Contig, rownames(current))
})

test_that("export_abundance reports a missing assay by name and lists the alternatives", {
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".csv")
  expect_error(export_abundance(tse, tmp, assay.type = "nonexistent"),
               "Assay 'nonexistent' not found")
  expect_error(export_abundance(tse, tmp, assay.type = "nonexistent"), "counts")
})

test_that("export_abundance warns rather than silently numbering unnamed rows", {
  tse <- make_tse(rownames = FALSE)
  tmp <- withr::local_tempfile(fileext = ".csv")
  expect_warning(export_abundance(tse, tmp), "no rownames")
})

# --------------------------------------------------------------------------
# export_annotations
# --------------------------------------------------------------------------

test_that("export_annotations TSV round-trips every rowData column", {
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".tsv")
  export_annotations(tse, tmp, format = "tsv")

  back <- utils::read.delim(tmp, check.names = FALSE)
  rd <- as.data.frame(SummarizedExperiment::rowData(tse))
  expect_equal(colnames(back), c("Contig", colnames(rd)))
  expect_equal(nrow(back), nrow(rd))
  expect_equal(back$Contig, rownames(tse))
  expect_equal(back$checkv_quality, rd$checkv_quality)
  expect_equal(back$checkv_contig_length, rd$checkv_contig_length)
  ## NA must survive as NA, not as the string "NA".
  expect_equal(is.na(back$Family), is.na(rd$Family))
})

test_that("export_annotations CSV round-trips every rowData column", {
  tse <- make_tse()
  tmp <- withr::local_tempfile(fileext = ".csv")
  export_annotations(tse, tmp, format = "csv")
  back <- utils::read.csv(tmp, check.names = FALSE)
  rd <- as.data.frame(SummarizedExperiment::rowData(tse))
  expect_equal(colnames(back), c("Contig", colnames(rd)))
  expect_equal(nrow(back), nrow(rd))
})

test_that("an embedded tab does not shift the columns of a TSV export", {
  tse <- make_tse()
  SummarizedExperiment::rowData(tse)$note <-
    c("plain", "contains\ttab", "contains\nnewline", "carriage\rreturn", "fine")

  tmp <- withr::local_tempfile(fileext = ".tsv")
  expect_warning(export_annotations(tse, tmp, format = "tsv"), "note")

  back <- utils::read.delim(tmp, check.names = FALSE)
  ## Without sanitizing, the newline would add a row and the tab a column.
  expect_equal(nrow(back), nrow(tse))
  expect_equal(ncol(back), ncol(SummarizedExperiment::rowData(tse)) + 1L)
  expect_false(any(grepl("[\t\r\n]", back$note)))
})

test_that("export_annotations flattens list columns instead of failing", {
  tse <- make_tse()
  SummarizedExperiment::rowData(tse)$hits <-
    IRanges::CharacterList(list(c("a", "b"), "c", character(0), "d", c("e", "f")))

  tmp <- withr::local_tempfile(fileext = ".tsv")
  expect_silent(export_annotations(tse, tmp, format = "tsv"))
  back <- utils::read.delim(tmp, check.names = FALSE)
  expect_equal(back$hits[1], "a;b")
  expect_equal(nrow(back), nrow(tse))
})

test_that("export_annotations handles rowData with no columns", {
  tse <- make_tse(rowdata = FALSE)
  tmp <- withr::local_tempfile(fileext = ".tsv")
  export_annotations(tse, tmp)
  back <- utils::read.delim(tmp, check.names = FALSE)
  expect_equal(colnames(back), "Contig")
  expect_equal(nrow(back), nrow(tse))
})

# --------------------------------------------------------------------------
# Path handling, shared by all three exporters
# --------------------------------------------------------------------------

test_that("exporters create a missing output directory", {
  tse <- make_tse()
  root <- withr::local_tempdir()
  nested <- file.path(root, "a", "b", "c")

  expect_true(file.exists(export_abundance(tse, file.path(nested, "ab.csv"))))
  expect_true(file.exists(export_annotations(tse, file.path(nested, "anno.tsv"))))
  expect_true(file.exists(export_vpftse(tse, file.path(nested, "tse.rds"))))
})

test_that("exporters handle paths containing spaces", {
  tse <- make_tse()
  root <- withr::local_tempdir()
  dir.create(file.path(root, "dir with space"))
  out <- file.path(root, "dir with space", "file with space.csv")
  expect_true(file.exists(export_abundance(tse, out)))
  expect_equal(nrow(utils::read.csv(out)), nrow(tse))
})

test_that("exporters reject unwritable destinations with an actionable message", {
  skip_on_os("windows")
  skip_if(as.integer(Sys.info()[["effective_user"]] == "root"), "running as root")
  tse <- make_tse()
  root <- withr::local_tempdir()
  ro <- file.path(root, "readonly")
  dir.create(ro)
  withr::defer(Sys.chmod(ro, "0755"))
  Sys.chmod(ro, "0555")
  skip_if(file.access(ro, mode = 2L) == 0L, "directory permissions not enforced")

  expect_error(export_abundance(tse, file.path(ro, "a.csv")), "not writable")
  expect_error(export_annotations(tse, file.path(ro, "a.tsv")), "not writable")
  expect_error(export_vpftse(tse, file.path(ro, "a.rds")), "not writable")
})

test_that("exporters reject a directory or an empty path", {
  tse <- make_tse()
  root <- withr::local_tempdir()
  expect_error(export_abundance(tse, root), "existing directory")
  expect_error(export_abundance(tse, ""), "non-empty file path")
  expect_error(export_abundance(tse, NA), "non-empty file path")
  expect_error(export_annotations(tse, c("a.tsv", "b.tsv")), "non-empty file path")
})
