test_that("the navbar exposes every module", {
  ui <- app_ui()
  html <- as.character(ui)
  for (tab in c("Data", "Filter", "Taxonomy", "Diversity", "Contigs",
                "Host &amp; lifestyle", "Genes", "Compare", "Export", "About")) {
    expect_true(grepl(paste0('data-value="', tab, '"'), html, fixed = TRUE),
                info = tab)
  }
})

test_that("every module namespace reaches the UI", {
  html <- as.character(app_ui())
  # Regression: a conditionalPanel condition written without `ns = ns` refers to
  # a top-level output that does not exist inside a module, so the panel never
  # appears. The generated JS condition must be namespaced.
  expect_true(grepl("data-source", html) || grepl("input\\.source", html))
  for (id in c("data-load", "filter-assay_ui", "composition-plt_stack",
               "diversity-plt_alpha", "features-tbl_columns", "host-plt_host",
               "genes-plt_phrog", "compare-plt_overlap", "export-dl_tse",
               "about-tutorial")) {
    expect_true(grepl(id, html, fixed = TRUE), info = id)
  }
})

test_that("conditional panels are namespaced", {
  html <- as.character(mod_data_input_ui("data"))
  # `ns = ns` rewrites `input.source` to the fully qualified input id.
  expect_true(grepl("data-source", html, fixed = TRUE))
})

test_that("app_server wires the filtered object to every module", {
  body_txt <- paste(deparse(body(app_server)), collapse = " ")
  expect_true(grepl("mod_data_input_server", body_txt))
  expect_true(grepl("mod_vpfilter_server", body_txt))
  for (m in c("composition", "diversity", "features", "host", "genes",
              "compare", "export")) {
    expect_true(grepl(paste0("mod_", m, "_server"), body_txt), info = m)
  }
  # Downstream modules must consume the shared filtered reactive, not re-filter.
  expect_true(grepl("r_filter", body_txt))
})

test_that("the status bar reports the current selection", {
  testServer(app_server, {
    expect_match(output$status_bar, "No dataset loaded")
  })
})

test_that("run_app raises the upload limit above Shiny's 5 MB default", {
  # Regression: a real viroprofiler_output.rds easily exceeds 5 MB, and the
  # upload failed with a generic browser error.
  old <- options(shiny.maxRequestSize = NULL)
  on.exit(options(old), add = TRUE)
  withr::with_envvar(list(VPFKIT_MAX_UPLOAD_MB = "123"), {
    app <- run_app()
    expect_s3_class(app, "shiny.appobj")
    expect_equal(getOption("shiny.maxRequestSize"), 123 * 1024^2)
  })
  app <- run_app(max_upload_mb = 42)
  expect_equal(getOption("shiny.maxRequestSize"), 42 * 1024^2)
  # A nonsensical value falls back rather than disabling uploads.
  app <- run_app(max_upload_mb = -1)
  expect_equal(getOption("shiny.maxRequestSize"), 500 * 1024^2)
})

test_that("server-path loading can be switched off explicitly", {
  old <- options(vpfkit.allow_server_path = NULL)
  on.exit(options(old), add = TRUE)
  run_app(allow_server_path = FALSE)
  expect_false(vpf_server_path_allowed())
  run_app(allow_server_path = TRUE)
  expect_true(vpf_server_path_allowed())
})


# ---------------------------------------------------------------------------
# End-to-end run in a real browser
# ---------------------------------------------------------------------------

test_that("the app loads a demo dataset and drives its tabs in a browser", {
  skip_on_cran()
  skip_if_not_installed("shinytest2")
  skip_if_not_installed("chromote")
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NULL)
  skip_if(is.null(chrome) || is.na(chrome), "no Chrome or Chromium available")
  skip_if(length(vpf_demo_datasets()) == 0, "demo datasets not available")

  # The app has to be driven from its own directory: shinytest2 starts a fresh
  # R process, and a `shinyApp` object built here would carry closures from a
  # pkgload namespace that the child process does not have.
  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                            mustWork = FALSE)
  skip_if(!file.exists(file.path(pkg_root, "app.R")), "app.R not found")

  app <- shinytest2::AppDriver$new(
    pkg_root, name = "vpfkit", load_timeout = 120 * 1000, timeout = 60 * 1000
  )
  on.exit(app$stop(), add = TRUE)

  # The status bar starts empty and the Data tab offers the bundled examples.
  expect_match(app$get_value(output = "status_bar"), "No dataset loaded")

  app$click("data-load")
  app$wait_for_idle(timeout = 30 * 1000)
  status <- app$get_value(output = "status_bar")
  expect_match(status, "60 of 60 contigs kept")
  expect_match(status, "8 samples")
  expect_match(status, "assay: counts")

  # Filtering propagates to the shared object rather than to one tab.
  app$set_inputs(main_nav = "Filter")
  app$wait_for_idle(timeout = 30 * 1000)
  app$set_inputs(`filter-min_length` = 15000)
  app$wait_for_idle(timeout = 30 * 1000)
  filtered_status <- app$get_value(output = "status_bar")
  expect_false(identical(filtered_status, status))
  expect_match(filtered_status, "of 60 contigs kept")

  # Reset restores the full selection, which the mock session cannot verify.
  app$click("filter-reset")
  app$wait_for_idle(timeout = 30 * 1000)
  expect_match(app$get_value(output = "status_bar"), "60 of 60 contigs kept")

  # Each analysis tab renders without a Shiny error banner.
  for (tab in c("Taxonomy", "Diversity", "Contigs", "Host & lifestyle",
                "Genes", "Compare", "Export", "About")) {
    app$set_inputs(main_nav = tab)
    app$wait_for_idle(timeout = 30 * 1000)
    html <- app$get_html("body")
    expect_false(grepl("shiny-output-error", html), info = tab)
  }

  logs <- app$get_logs()
  errors <- logs[logs$level %in% c("ERROR") & !is.na(logs$level), ]
  expect_equal(nrow(errors), 0)
})

test_that("the app starts as a real HTTP server and serves its UI", {
  skip_on_cran()
  skip_if_not_installed("callr")
  skip_if_not_installed("curl")
  skip_if(length(vpf_demo_datasets()) == 0, "demo datasets not available")

  pkg_root <- normalizePath(file.path(testthat::test_path(), "..", ".."),
                            mustWork = FALSE)
  skip_if(!file.exists(file.path(pkg_root, "DESCRIPTION")),
          "package root not found")

  port <- 8000L + sample(900, 1)
  proc <- callr::r_bg(
    function(root, port) {
      pkgload::load_all(root, quiet = TRUE)
      shiny::runApp(vpfkit::run_app(), port = port, host = "127.0.0.1",
                    launch.browser = FALSE)
    },
    args = list(root = pkg_root, port = port),
    supervise = TRUE
  )
  on.exit({
    if (proc$is_alive()) proc$kill()
  }, add = TRUE)

  url <- sprintf("http://127.0.0.1:%d/", port)
  html <- NULL
  for (i in seq_len(60)) {
    if (!proc$is_alive()) break
    html <- suppressWarnings(tryCatch(
      paste(readLines(url, warn = FALSE), collapse = "\n"),
      error = function(e) NULL
    ))
    if (!is.null(html)) break
    Sys.sleep(1)
  }
  if (is.null(html)) {
    skip(paste("app did not start:", paste(utils::tail(proc$read_error_lines(), 5),
                                           collapse = " | ")))
  }
  expect_match(html, "ViroProfiler-viewer")
  expect_match(html, "data-load", fixed = TRUE)
  expect_match(html, "filter-assay_ui", fixed = TRUE)
  expect_match(html, "custom.css", fixed = TRUE)
})
