vpf_exp_filter <- function(tse, assay = "counts", raw = tse, name = "demo.rds") {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(data.frame(
      Step = c("Loaded", "Contig length"),
      Applied = c("-", "yes"),
      Removed = c(0L, 2L),
      Remaining = c(10L, 8L),
      Detail = c("", ">= 3,000 bp"),
      stringsAsFactors = FALSE
    )),
    name = shiny::reactive(name)
  )
}

vpf_exp_demo <- function(pattern = "Gut virome") {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

test_that("module ui works", {
  ui <- mod_export_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_export_ui)))
})

test_that("report availability explains missing Quarto installations", {
  local_mocked_bindings(
    .vpf_quarto_available = function() FALSE,
    .vpf_quarto_hint = function() "Install Quarto to render reports."
  )
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(NULL)), {
    expect_match(output$report_availability, "Install Quarto to render reports")
  })
})

test_that("provenance describes the selection even with nothing loaded", {
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(NULL)), {
    txt <- provenance_text()
    expect_match(txt, "No dataset loaded")
    expect_match(output$selection_summary, "ViroProfiler-viewer export provenance")
  })
})

test_that("provenance records the assay semantics and every filter step", {
  tse <- vpf_exp_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(tse)), {
    txt <- provenance_text()
    expect_match(txt, "Contigs after filtering")
    expect_match(txt, "Active abundance assay:   counts")
    expect_match(txt, "Contig length")
    expect_match(txt, "CoverM read count", fixed = TRUE)
    expect_match(txt, "R version")
  })
})

test_that("download handlers write the files they promise", {
  tse <- vpf_exp_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(tse)), {
    rds <- tempfile(fileext = ".rds")
    write_tse(rds)
    expect_true(file.exists(rds))
    expect_s4_class(readRDS(rds), "TreeSummarizedExperiment")
    expect_match(name_tse(), "^demo_filtered_.*\\.rds$")

    csv <- tempfile(fileext = ".csv")
    write_abundance(csv, "csv")
    ab <- utils::read.csv(csv)
    expect_equal(nrow(ab), nrow(tse))
    expect_equal(colnames(ab)[1], "Contig")
    expect_match(name_abundance_csv(), "abundance_counts")

    tsv <- tempfile(fileext = ".tsv")
    write_annotations(tsv, "tsv")
    an <- utils::read.delim(tsv)
    expect_equal(nrow(an), nrow(tse))

    meta <- tempfile(fileext = ".csv")
    write_metadata(meta)
    md <- utils::read.csv(meta)
    expect_equal(nrow(md), ncol(tse))
    expect_equal(colnames(md)[1], "Sample")

    prov <- tempfile(fileext = ".txt")
    write_provenance(prov)
    expect_true(any(grepl("provenance", readLines(prov), ignore.case = TRUE)))

    xlsx <- tempfile(fileext = ".xlsx")
    write_abundance(xlsx, "xlsx")
    expect_true(file.exists(xlsx))
    expect_gt(file.size(xlsx), 0)

    unlink(c(rds, csv, tsv, meta, prov, xlsx))
  })
})

test_that("the gene-annotation download explains its own absence", {
  tse <- vpf_exp_demo()
  skip_if(is.null(tse), "demo dataset not available")
  ga <- vpf_gene_annotations(tse)
  skip_if(is.null(ga), "demo dataset has no gene annotations")
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(tse)), {
    f <- tempfile(fileext = ".tsv")
    write_genes(f)
    expect_equal(nrow(utils::read.delim(f)), nrow(ga))
    unlink(f)
  })

  stripped <- tse
  S4Vectors::metadata(stripped) <- list()
  testServer(mod_export_server, args = list(r_filter = vpf_exp_filter(stripped)), {
    f <- tempfile(fileext = ".tsv")
    expect_error(write_genes(f), "gene-level annotations")
    unlink(f)
  })
})

test_that("downloads refuse an empty selection with a readable message", {
  tse <- vpf_exp_demo()
  skip_if(is.null(tse), "demo dataset not available")
  empty <- tse[integer(0), ]
  testServer(mod_export_server,
             args = list(r_filter = vpf_exp_filter(empty, raw = tse)), {
    f <- tempfile(fileext = ".rds")
    expect_error(write_tse(f), "No contigs pass")
    unlink(f)
  })
})

test_that("file names are sanitized and dated", {
  tse <- vpf_exp_demo()
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_export_server,
             args = list(r_filter = vpf_exp_filter(tse, name = "my run 2024.rds")), {
    expect_equal(stem(), "my_run_2024")
    expect_match(name_annotations_csv(),
                 paste0("^my_run_2024_annotations_", Sys.Date(), "\\.csv$"))
  })
})
