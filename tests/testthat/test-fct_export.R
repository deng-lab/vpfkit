# Helper to create a test TSE from fixtures
.make_test_tse <- function() {
  fix <- test_path("fixtures")
  create_vpftse(
    fin_abcount = file.path(fix, "coverm_abundance_count.tsv"),
    fin_abtpm = file.path(fix, "coverm_abundance_tpm.tsv"),
    fin_abtmm = file.path(fix, "coverm_abundance_tmm.tsv"),
    fin_abcov = file.path(fix, "coverm_abundance_covfrac.tsv"),
    fin_taxa = file.path(fix, "mmseqs2_taxa.tsv"),
    fin_checkv = file.path(fix, "checkv_quality_summary.tsv"),
    fin_virsorter2 = file.path(fix, "virsorter2_final_viral_score.tsv"),
    fin_vibrant = file.path(fix, "vibrant_quality.tsv"),
    fin_dvf = file.path(fix, "dvf_virus.tsv"),
    fin_replicyc = file.path(fix, "bacphlip.bacphlip")
  )
}

test_that("export_vpftse writes and reads back correctly", {
  tse <- .make_test_tse()
  tmp <- withr::local_tempfile(fileext = ".rds")
  export_vpftse(tse, tmp)
  expect_true(file.exists(tmp))
  loaded <- readRDS(tmp)
  expect_s4_class(loaded, "TreeSummarizedExperiment")
  expect_equal(nrow(loaded), 5)
})

test_that("export_abundance writes CSV", {
  tse <- .make_test_tse()
  tmp <- withr::local_tempfile(fileext = ".csv")
  export_abundance(tse, tmp, assay.type = "counts", format = "csv")
  expect_true(file.exists(tmp))
  df <- read.csv(tmp)
  expect_equal(nrow(df), 5)
  expect_true("Contig" %in% colnames(df))
})

test_that("export_abundance writes XLSX", {
  tse <- .make_test_tse()
  tmp <- withr::local_tempfile(fileext = ".xlsx")
  export_abundance(tse, tmp, assay.type = "counts", format = "xlsx")
  expect_true(file.exists(tmp))
  df <- openxlsx::read.xlsx(tmp)
  expect_equal(nrow(df), 5)
})

test_that("export_annotations writes TSV", {
  tse <- .make_test_tse()
  tmp <- withr::local_tempfile(fileext = ".tsv")
  export_annotations(tse, tmp, format = "tsv")
  expect_true(file.exists(tmp))
  df <- read.delim(tmp)
  expect_equal(nrow(df), 5)
  expect_true("checkv_quality" %in% colnames(df))
})

test_that("export_annotations writes CSV", {
  tse <- .make_test_tse()
  tmp <- withr::local_tempfile(fileext = ".csv")
  export_annotations(tse, tmp, format = "csv")
  expect_true(file.exists(tmp))
  df <- read.csv(tmp)
  expect_equal(nrow(df), 5)
})
