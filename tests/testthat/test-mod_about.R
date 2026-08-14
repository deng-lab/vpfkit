test_that("module ui works", {
  ui <- mod_about_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_about_ui)))
})

test_that("the tutorial ships with the package and is found via app_sys", {
  # Regression: the commented-out tab used includeMarkdown("inst/app/www/...").
  # `inst/` is stripped at install time, so that path resolves only in a source
  # checkout and the tab would have been empty in an installed package.
  path <- app_sys("app/www/tutorial.md")
  expect_true(nzchar(path))
  expect_true(file.exists(path))
  txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(txt, "ViroProfiler-viewer")
  # The tutorial has to describe the app that exists now.
  for (tab in c("Data", "Filter", "Taxonomy", "Diversity", "Contigs",
                "Host & lifestyle", "Genes", "Compare", "Export")) {
    expect_match(txt, tab, fixed = TRUE)
  }
  expect_match(txt, "covfrac", fixed = TRUE)
  expect_match(txt, "not edgeR's TMM", fixed = TRUE)
})

test_that("the footer ships with the package", {
  path <- app_sys("app/www/footer.html")
  expect_true(nzchar(path))
  expect_true(file.exists(path))
})

test_that("the tutorial and footer render", {
  testServer(mod_about_server, {
    expect_false(is.null(output$tutorial))
    expect_false(is.null(output$footer))
    expect_match(output$versions, "vpfkit")
    expect_match(output$versions, "R version")
  })
})

test_that("markdown rendering degrades rather than failing", {
  f <- tempfile(fileext = ".md")
  writeLines(c("# Title", "", "Some *text*."), f)
  on.exit(unlink(f), add = TRUE)
  out <- vpf_render_markdown(f)
  expect_true(inherits(out, "html") || inherits(out, "shiny.tag"))
})
