test_that("generate_report creates HTML output", {
  skip_if_not(requireNamespace("quarto", quietly = TRUE), "quarto package not available")
  skip_if_not(nzchar(Sys.which("quarto")), "quarto CLI not installed")

  fix <- test_path("fixtures")
  tse <- create_vpftse(
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

  tmp_html <- withr::local_tempfile(fileext = ".html")
  result <- generate_report(tse, tmp_html)
  expect_true(file.exists(tmp_html))
  expect_equal(result, tmp_html)
  # Check it's valid HTML
  content <- readLines(tmp_html, n = 5)
  expect_true(any(grepl("<html|<!DOCTYPE", content, ignore.case = TRUE)))
})

test_that("generate_report errors without quarto package", {
  # This test just verifies the error message is helpful
  skip("Skipping: cannot unload quarto if installed")
})
