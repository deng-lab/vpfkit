vpf_host_filter <- function(tse, assay = "counts", raw = tse) {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive("test")
  )
}

vpf_host_demo <- function(pattern = "Gut virome") {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_host_real <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("module ui works", {
  ui <- mod_host_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_host_ui)))
})

test_that("no dataset renders placeholders", {
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(NULL)), {
    expect_null(host_column())
    expect_null(cycle_column())
    expect_false(is.null(output$host_status))
    expect_false(is.null(output$cycle_status))
    expect_false(is.null(output$plt_host))
    expect_false(is.null(output$plt_cycle))
    expect_false(is.null(output$plt_provirus))
  })
})

test_that("iPHoP predictions are summarized when present", {
  tse <- vpf_host_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(tse)), {
    expect_equal(host_column(), "iphop_genus")
    hd <- host_data()
    expect_true(nrow(hd) > 0)
    expect_true("Score" %in% colnames(hd))
    expect_false(is.null(output$plt_host))
    expect_false(is.null(output$plt_host_score))
    expect_false(is.null(output$tbl_host))
    expect_false(is.null(output$host_controls))
  })
})

test_that("the host score threshold filters the table", {
  tse <- vpf_host_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(tse)), {
    session$setInputs(min_host_score = 0)
    n_all <- nrow(host_data())
    session$setInputs(min_host_score = 99)
    expect_lte(nrow(host_data()), n_all)
  })
})

test_that("replication cycle is summarized and weighted by abundance", {
  tse <- vpf_host_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(tse)), {
    expect_equal(cycle_column(), "bacphlip_replicyc")
    cd <- cycle_data()
    expect_true("Abundance" %in% colnames(cd))
    expect_equal(nrow(cd), nrow(tse))
    expect_false(is.null(output$plt_cycle))
    expect_false(is.null(output$plt_cycle_abundance))
    expect_false(is.null(output$tbl_cycle))
  })
})

test_that("a dataset without host prediction says which option produces it", {
  tse <- vpf_host_real()
  skip_if(is.null(tse), "real ViroProfiler output not available")
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(tse)), {
    # The real run had no iPHoP output.
    expect_null(host_column())
    expect_false(is.null(output$host_status))
    expect_false(is.null(output$plt_host))
    expect_false(is.null(output$plt_host_score))
    expect_false(is.null(output$tbl_host))
    # But it does carry BACPHLIP and VIBRANT lifestyle calls.
    expect_equal(cycle_column(), "bacphlip_replicyc")
    expect_false(is.null(output$plt_cycle))
    # And CheckV provirus status.
    expect_equal(provirus_column(), "checkv_provirus")
    expect_false(is.null(output$plt_provirus))
  })
})

test_that("a dataset without lifestyle or provirus columns degrades cleanly", {
  tse <- vpf_host_demo()
  skip_if(is.null(tse), "demo dataset not available")
  keep <- setdiff(colnames(SummarizedExperiment::rowData(tse)),
                  c("bacphlip_replicyc", "replidec_replicyc", "vibrant_replicyc",
                    "checkv_provirus", "genomad_topology"))
  SummarizedExperiment::rowData(tse) <-
    SummarizedExperiment::rowData(tse)[, keep, drop = FALSE]
  testServer(mod_host_server, args = list(r_filter = vpf_host_filter(tse)), {
    expect_null(cycle_column())
    expect_null(provirus_column())
    expect_false(is.null(output$cycle_status))
    expect_false(is.null(output$provirus_status))
    expect_false(is.null(output$plt_cycle))
    expect_false(is.null(output$plt_cycle_abundance))
    expect_false(is.null(output$plt_provirus))
    expect_false(is.null(output$tbl_cycle))
  })
})

test_that("zero surviving contigs does not error", {
  tse <- vpf_host_demo()
  skip_if(is.null(tse), "demo dataset not available")
  empty <- tse[integer(0), ]
  testServer(mod_host_server,
             args = list(r_filter = vpf_host_filter(empty, raw = tse)), {
    expect_null(rows())
    expect_false(is.null(output$plt_host))
    expect_false(is.null(output$plt_cycle))
    expect_false(is.null(output$plt_provirus))
  })
})
