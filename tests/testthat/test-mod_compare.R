vpf_cmp_filter <- function(tse, assay = "counts", raw = tse, name = "Dataset 1") {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive(name)
  )
}

vpf_cmp_demo <- function(pattern) {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_cmp_demo_path <- function(pattern) {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

test_that("module ui works", {
  ui <- mod_compare_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_compare_ui)))
})

test_that("nothing loaded renders placeholders", {
  testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(NULL)), {
    expect_null(second())
    expect_null(overlap())
    expect_false(is.null(output$summary))
    expect_false(is.null(output$plt_overlap))
    expect_false(is.null(output$plt_alpha))
    expect_false(is.null(output$plt_composition))
    expect_false(is.null(output$contig_overlap))
  })
})

test_that("a second dataset loads and is compared at a shared rank", {
  a <- vpf_cmp_demo("Gut virome")
  b_path <- vpf_cmp_demo_path("Environmental")
  skip_if(is.null(a) || is.na(b_path), "demo datasets not available")
  testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(a)), {
    session$setInputs(source = "demo",
                      demo_choice = "Environmental virome, 40 contigs x 6 samples (Ocean vs Soil)",
                      alpha_index = "shannon")
    session$setInputs(load = 1)
    expect_s4_class(second(), "TreeSummarizedExperiment")
    expect_true(length(common_ranks()) > 0)
    session$setInputs(rank = "Family")
    ov <- overlap()
    expect_equal(ov$rank, "Family")
    expect_true(length(ov$shared) + length(ov$only_a) + length(ov$only_b) > 0)
    expect_false(is.null(output$plt_overlap))
    expect_false(is.null(output$tbl_overlap))
    expect_false(is.null(output$plt_alpha))
    expect_false(is.null(output$plt_composition))
    expect_false(is.null(output$summary))
  })
})

test_that("independent runs sharing no contig identifier are explained", {
  a <- vpf_cmp_demo("Gut virome")
  b_path <- vpf_cmp_demo_path("Environmental")
  skip_if(is.null(a) || is.na(b_path), "demo datasets not available")
  testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(a)), {
    session$setInputs(source = "demo",
                      demo_choice = "Environmental virome, 40 contigs x 6 samples (Ocean vs Soil)")
    session$setInputs(load = 1)
    expect_length(intersect(rownames(a), rownames(second())), 0)
    expect_false(is.null(output$contig_overlap))
  })
})

test_that("shared contig identifiers are listed when they exist", {
  a <- vpf_cmp_demo("Gut virome")
  skip_if(is.null(a), "demo dataset not available")
  path <- vpf_cmp_demo_path("Gut virome")
  testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(a)), {
    session$setInputs(source = "demo",
                      demo_choice = "Gut virome, 60 contigs x 8 samples (Healthy vs Disease)")
    session$setInputs(load = 1)
    expect_equal(length(intersect(rownames(a), rownames(second()))), nrow(a))
    expect_false(is.null(output$contig_overlap))
  })
})

test_that("a bad second file is refused without crashing", {
  a <- vpf_cmp_demo("Gut virome")
  skip_if(is.null(a), "demo dataset not available")
  bad <- tempfile(fileext = ".rds")
  saveRDS(1:3, bad)
  on.exit(unlink(bad), add = TRUE)
  testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(a)), {
    session$setInputs(source = "upload")
    session$setInputs(upload = data.frame(name = "bad.rds", size = 1, type = "",
                                          datapath = bad, stringsAsFactors = FALSE))
    session$setInputs(load = 1)
    expect_null(second())
    expect_false(is.null(output$load_status))
  })
})

test_that("two objects with no shared rank are reported", {
  a <- vpf_cmp_demo("Gut virome")
  skip_if(is.null(a), "demo dataset not available")
  stripped <- a
  SummarizedExperiment::rowData(stripped) <-
    SummarizedExperiment::rowData(stripped)[, "checkv_quality", drop = FALSE]
  path <- tempfile(fileext = ".rds")
  saveRDS(a, path)
  on.exit(unlink(path), add = TRUE)
  withr::with_options(list(vpfkit.allow_server_path = TRUE), {
    testServer(mod_compare_server, args = list(r_filter = vpf_cmp_filter(stripped)), {
      session$setInputs(source = "path", path = path)
      session$setInputs(load = 1)
      expect_length(common_ranks(), 0)
      expect_null(overlap())
      expect_false(is.null(output$summary))
      expect_false(is.null(output$plt_overlap))
      expect_false(is.null(output$plt_composition))
    })
  })
})
