vpf_fake_filter <- function(tse, assay = "counts", raw = tse, name = "test") {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive(name)
  )
}

vpf_demo_tse <- function(pattern) {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_real_tse <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("module ui works", {
  ui <- mod_composition_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_composition_ui)))
})

test_that("no dataset renders placeholders rather than errors", {
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(NULL)), {
    expect_false(is.null(output$plt_stack))
    expect_false(is.null(output$plt_ranks))
    expect_false(is.null(output$status))
    expect_equal(output$stack_caption, "")
  })
})

test_that("composition aggregates by rank and keeps unclassified contigs", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Family", scale = "relative", top_n = 12,
                      group = "", drop_unclassified = FALSE)
    expect_true("Unclassified" %in% rownames(composition()$matrix))
    expect_equal(composition()$rank, "Family")
    expect_false(is.null(output$plt_stack))
    expect_false(is.null(output$plt_ranks))
    expect_match(output$stack_caption, "filtered subset")
    # Relative abundance must sum to 1 per sample.
    pd <- plot_data()$df
    sums <- tapply(pd$Value, pd$Sample, sum)
    expect_true(all(abs(sums - 1) < 1e-8))
  })
})

test_that("hiding unclassified contigs is announced", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Family", scale = "relative", top_n = 12,
                      group = "", drop_unclassified = TRUE)
    expect_false("Unclassified" %in% rownames(composition()$matrix))
    expect_false(is.null(output$unclassified_warning))
  })
})

test_that("top-n pooling produces an Other category", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Genus", scale = "relative", top_n = 3,
                      group = "", drop_unclassified = FALSE)
    expect_true("Other" %in% levels(plot_data()$df$Taxon))
  })
})

test_that("grouping adds a facet column", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Family", scale = "absolute", top_n = 8,
                      group = "group", drop_unclassified = FALSE)
    pd <- plot_data()$df
    expect_true("Group" %in% colnames(pd))
    expect_setequal(unique(pd$Group), c("Healthy", "Disease"))
    expect_false(is.null(output$plt_stack))
  })
})

test_that("blank taxonomy strings become Unclassified, not an empty taxon", {
  # The real pipeline writes "" rather than NA when a rank is unassigned.
  tse <- vpf_real_tse()
  skip_if(is.null(tse), "real ViroProfiler output not available")
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Species", scale = "relative", top_n = 12,
                      group = "", drop_unclassified = FALSE)
    taxa <- rownames(composition()$matrix)
    expect_false(any(trimws(taxa) == ""))
    expect_true("Unclassified" %in% taxa)
    expect_false(is.null(output$plt_stack))
    expect_false(is.null(output$plt_ranks))
  })
})

test_that("a dataset with no taxonomy explains itself", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  SummarizedExperiment::rowData(tse) <-
    SummarizedExperiment::rowData(tse)[, "checkv_quality", drop = FALSE]
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    expect_length(ranks(), 0)
    expect_null(composition())
    expect_false(is.null(output$status))
    expect_false(is.null(output$plt_stack))
  })
})

test_that("zero surviving contigs does not error", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  empty <- tse[integer(0), ]
  testServer(mod_composition_server,
             args = list(r_filter = vpf_fake_filter(empty, raw = tse)), {
    expect_null(composition())
    expect_false(is.null(output$plt_stack))
    expect_false(is.null(output$plt_ranks))
  })
})

test_that("samples with zero total abundance do not produce NaN bars", {
  tse <- vpf_demo_tse("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  SummarizedExperiment::assay(tse, "counts")[, 1] <- 0L
  testServer(mod_composition_server, args = list(r_filter = vpf_fake_filter(tse)), {
    session$setInputs(rank = "Family", scale = "relative", top_n = 12,
                      group = "", drop_unclassified = FALSE)
    pd <- plot_data()$df
    expect_false(anyNA(pd$Value))
    expect_false(is.null(output$plt_stack))
  })
})
