vpf_feat_filter <- function(tse, assay = "counts", raw = tse) {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive("test")
  )
}

vpf_feat_demo <- function(pattern = "Gut virome") {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_feat_real <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("module ui works", {
  ui <- mod_features_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_features_ui)))
})

test_that("no dataset renders placeholders", {
  testServer(mod_features_server, args = list(r_filter = vpf_feat_filter(NULL)), {
    expect_null(annotation_table())
    expect_false(is.null(output$status))
    expect_false(is.null(output$tbl_columns))
    expect_false(is.null(output$tbl_contig))
    expect_false(is.null(output$plt_contig))
  })
})

test_that("the annotation table carries every rowData column plus the contig id", {
  tse <- vpf_feat_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_features_server, args = list(r_filter = vpf_feat_filter(tse)), {
    session$setInputs(include_abundance = FALSE)
    df <- annotation_table()
    expect_equal(nrow(df), nrow(tse))
    expect_equal(colnames(df)[1], "Contig")
    expect_true(all(colnames(SummarizedExperiment::rowData(tse)) %in% colnames(df)))
    expect_false(is.null(output$tbl_columns))
  })
})

test_that("per-sample abundance columns can be appended", {
  tse <- vpf_feat_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_features_server, args = list(r_filter = vpf_feat_filter(tse)), {
    session$setInputs(include_abundance = TRUE)
    df <- annotation_table()
    expect_true(any(grepl("\\[counts\\]$", colnames(df))))
    expect_equal(sum(grepl("\\[counts\\]$", colnames(df))), ncol(tse))
  })
})

test_that("the single-contig view resolves a default and its abundance", {
  tse <- vpf_feat_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_features_server, args = list(r_filter = vpf_feat_filter(tse)), {
    session$setInputs(include_abundance = FALSE)
    expect_equal(selected_contig(), rownames(tse)[[1]])
    expect_equal(length(contig_ids()), nrow(tse))
    expect_false(is.null(output$plt_contig))
    expect_false(is.null(output$tbl_contig))
    expect_false(is.null(output$contig_summary))
    expect_false(is.null(output$contig_hint))
    # A stale selection from a previous dataset must not crash the module.
    session$setInputs(contig = "not_a_contig")
    expect_equal(selected_contig(), rownames(tse)[[1]])
    # Nor does the empty string selectize sends before it is populated.
    session$setInputs(contig = "")
    expect_equal(selected_contig(), rownames(tse)[[1]])
  })
})

test_that("real ViroProfiler output browses without error", {
  tse <- vpf_feat_real()
  skip_if(is.null(tse), "real ViroProfiler output not available")
  testServer(mod_features_server, args = list(r_filter = vpf_feat_filter(tse)), {
    session$setInputs(include_abundance = TRUE)
    df <- annotation_table()
    expect_equal(nrow(df), 18)
    expect_true("checkv_quality" %in% colnames(df))
    expect_false(is.null(output$tbl_columns))
    expect_false(is.null(output$plt_contig))
    # Columns that are entirely missing must be visible as such.
    expect_false(is.null(output$tbl_contig))
  })
})

test_that("zero surviving contigs does not error", {
  tse <- vpf_feat_demo()
  skip_if(is.null(tse), "demo dataset not available")
  empty <- tse[integer(0), ]
  testServer(mod_features_server,
             args = list(r_filter = vpf_feat_filter(empty, raw = tse)), {
    expect_null(annotation_table())
    expect_null(selected_contig())
    expect_length(contig_ids(), 0)
    expect_false(is.null(output$status))
    expect_false(is.null(output$tbl_columns))
    expect_false(is.null(output$plt_contig))
    expect_false(is.null(output$contig_hint))
  })
})
