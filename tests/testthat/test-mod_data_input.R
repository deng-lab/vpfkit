vpf_demo_path <- function(pattern) {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

vpf_real_path <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }
  candidate <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  if (file.exists(candidate)) candidate else NA_character_
}

test_that("module ui works", {
  ui <- mod_data_input_ui(id = "test")
  golem::expect_shinytaglist(ui)
  fmls <- formals(mod_data_input_ui)
  for (i in c("id")) {
    expect_true(i %in% names(fmls))
  }
})

test_that("nothing is loaded before the load button is pressed", {
  testServer(mod_data_input_server, {
    session$setInputs(source = "demo")
    expect_null(session$getReturned()$tse())
    # Every output must still render a placeholder rather than erroring.
    expect_false(is.null(output$tbl_annotations))
    expect_false(is.null(output$tbl_assays))
    expect_false(is.null(output$tbl_coldata))
    expect_false(is.null(output$plt_libsize))
    expect_false(is.null(output$overview))
  })
})

test_that("a bundled demo dataset loads and is described", {
  demo <- vpf_demo_path("Gut virome")
  skip_if(is.na(demo), "demo dataset not available")
  testServer(mod_data_input_server, {
    session$setInputs(source = "demo",
                      demo_choice = "Gut virome, 60 contigs x 8 samples (Healthy vs Disease)")
    session$setInputs(load = 1)
    tse <- session$getReturned()$tse()
    expect_s4_class(tse, "TreeSummarizedExperiment")
    expect_equal(dim(tse), c(60L, 8L))
    expect_false(is.null(output$overview))
    expect_false(is.null(output$tbl_annotations))
    expect_false(is.null(output$tbl_assays))
    expect_false(is.null(output$tbl_coldata))
    expect_false(is.null(output$plt_libsize))
  })
})

test_that("a missing upload produces a message rather than an error", {
  testServer(mod_data_input_server, {
    session$setInputs(source = "upload")
    session$setInputs(load = 1)
    expect_null(session$getReturned()$tse())
    expect_false(is.null(output$load_status))
  })
})

test_that("a file that is not a TSE is rejected with an explanation", {
  bad <- tempfile(fileext = ".rds")
  saveRDS(data.frame(a = 1), bad)
  on.exit(unlink(bad), add = TRUE)
  testServer(mod_data_input_server, {
    session$setInputs(source = "upload")
    session$setInputs(upload = data.frame(name = "bad.rds", size = 1,
                                          type = "", datapath = bad,
                                          stringsAsFactors = FALSE))
    session$setInputs(load = 1)
    expect_null(session$getReturned()$tse())
  })
})

test_that("server-path loading honours the production switch", {
  demo <- vpf_demo_path("Gut virome")
  skip_if(is.na(demo), "demo dataset not available")
  withr::with_options(list(vpfkit.allow_server_path = FALSE), {
    testServer(mod_data_input_server, {
      session$setInputs(source = "path", path = demo)
      session$setInputs(load = 1)
      expect_null(session$getReturned()$tse())
    })
  })
  withr::with_options(list(vpfkit.allow_server_path = TRUE), {
    testServer(mod_data_input_server, {
      session$setInputs(source = "path", path = demo)
      session$setInputs(load = 1)
      expect_s4_class(session$getReturned()$tse(), "TreeSummarizedExperiment")
    })
  })
})

test_that("sample metadata is joined onto colData and unlocks grouping", {
  real <- vpf_real_path()
  skip_if(is.na(real), "real ViroProfiler output not available")
  meta <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(sample = c("HT02", "UC20"), condition = c("healthy", "colitis"),
               stringsAsFactors = FALSE),
    meta, row.names = FALSE
  )
  on.exit(unlink(meta), add = TRUE)

  withr::with_options(list(vpfkit.allow_server_path = TRUE), {
    testServer(mod_data_input_server, {
      session$setInputs(source = "path", path = real)
      session$setInputs(load = 1)
      tse <- session$getReturned()$tse()
      expect_equal(dim(tse), c(18L, 2L))
      # The real object carries only sample_name, so no grouping is possible.
      expect_length(vpf_group_candidates(tse), 0)

      session$setInputs(metadata = data.frame(name = "meta.csv", size = 1,
                                              type = "", datapath = meta,
                                              stringsAsFactors = FALSE))
      session$setInputs(apply_metadata = 1)
      joined <- session$getReturned()$tse()
      expect_true("condition" %in%
                    colnames(SummarizedExperiment::colData(joined)))
      expect_equal(vpf_group_candidates(joined), "condition")

      session$setInputs(clear_metadata = 1)
      expect_false("condition" %in%
                     colnames(SummarizedExperiment::colData(session$getReturned()$tse())))
    })
  })
})

test_that("metadata whose identifiers do not match is refused", {
  demo <- vpf_demo_path("Gut virome")
  skip_if(is.na(demo), "demo dataset not available")
  meta <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = c("nope_1", "nope_2"), g = c("a", "b")),
                   meta, row.names = FALSE)
  on.exit(unlink(meta), add = TRUE)
  testServer(mod_data_input_server, {
    session$setInputs(source = "demo",
                      demo_choice = "Gut virome, 60 contigs x 8 samples (Healthy vs Disease)")
    session$setInputs(load = 1)
    before <- ncol(SummarizedExperiment::colData(session$getReturned()$tse()))
    session$setInputs(metadata = data.frame(name = "meta.csv", size = 1, type = "",
                                            datapath = meta, stringsAsFactors = FALSE))
    session$setInputs(apply_metadata = 1)
    after <- ncol(SummarizedExperiment::colData(session$getReturned()$tse()))
    expect_equal(before, after)
  })
})

vpf_enriched_path <- function() {
  p <- Sys.getenv("VPFKIT_ENRICHED_TSE", unset = "")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }
  candidate <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output_enriched.rds"
  if (file.exists(candidate)) candidate else NA_character_
}

test_that("the enriched real object is read with its canonical assay names", {
  path <- vpf_enriched_path()
  skip_if(is.na(path), "enriched ViroProfiler output not available")
  withr::with_options(list(vpfkit.allow_server_path = TRUE), {
    testServer(mod_data_input_server, {
      session$setInputs(source = "path", path = path)
      session$setInputs(load = 1)
      tse <- session$getReturned()$tse()
      expect_equal(dim(tse), c(18L, 2L))
      expect_true("trimmed_mean" %in% SummarizedExperiment::assayNames(tse))
      expect_false("tmm" %in% SummarizedExperiment::assayNames(tse))
      # Nothing had to be renamed, so no legacy-name warning is shown.
      expect_length(renamed_assays(), 0)
      # This object ships a grouping variable, unlike the earlier one.
      expect_true("group" %in% vpf_group_candidates(tse))
      # Assay descriptions come from the object itself.
      info <- vpf_assay_info(tse)
      expect_match(info$description[info$assay == "trimmed_mean"], "CoverM")
      for (o in c("overview", "tbl_annotations", "tbl_assays", "tbl_coldata",
                  "plt_libsize")) {
        expect_false(is.null(output[[o]]), info = o)
      }
    })
  })
})

test_that("the legacy real object is renamed on read and the user is told", {
  path <- vpf_real_path()
  skip_if(is.na(path), "real ViroProfiler output not available")
  withr::with_options(list(vpfkit.allow_server_path = TRUE), {
    testServer(mod_data_input_server, {
      session$setInputs(source = "path", path = path)
      session$setInputs(load = 1)
      tse <- session$getReturned()$tse()
      expect_true("trimmed_mean" %in% SummarizedExperiment::assayNames(tse))
      expect_equal(unname(renamed_assays()), "trimmed_mean")
      expect_equal(names(renamed_assays()), "tmm")
      expect_false(is.null(output$assay_help))
    })
  })
})
