vpf_fake_input <- function(tse, name = "test") {
  list(tse = shiny::reactive(tse), name = shiny::reactive(name))
}

vpf_filter_tse <- function(n_row = 8, n_col = 4) {
  set.seed(42)
  counts <- matrix(as.integer(stats::rpois(n_row * n_col, 30)), nrow = n_row)
  rownames(counts) <- paste0("contig_", seq_len(n_row))
  colnames(counts) <- paste0("S", seq_len(n_col))
  covfrac <- matrix(seq(0, 1, length.out = n_row * n_col), nrow = n_row,
                    dimnames = dimnames(counts))
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts, covfrac = covfrac),
    rowData = S4Vectors::DataFrame(
      checkv_contig_length = seq(1000, by = 2000, length.out = n_row),
      checkv_quality = factor(rep(c("Complete", "Low-quality"), length.out = n_row),
                              levels = c("Complete", "Low-quality")),
      checkv_completeness = c(NA_real_, seq(10, 100, length.out = n_row - 1)),
      genomad_score = seq(0.1, 0.95, length.out = n_row),
      virsorter2_max_score = seq(0.2, 0.99, length.out = n_row),
      Family = rep(c("Siphoviridae", NA_character_), length.out = n_row),
      row.names = rownames(counts)
    ),
    colData = S4Vectors::DataFrame(
      sample_name = colnames(counts),
      group = rep(c("A", "B"), length.out = n_col),
      row.names = colnames(counts)
    )
  )
}

vpf_real_tse_path <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (nzchar(p) && file.exists(p)) {
    return(p)
  }
  candidate <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  if (file.exists(candidate)) candidate else NA_character_
}

test_that("module ui works", {
  ui <- mod_vpfilter_ui(id = "test")
  golem::expect_shinytaglist(ui)
  fmls <- formals(mod_vpfilter_ui)
  for (i in c("id")) {
    expect_true(i %in% names(fmls))
  }
})

test_that("no dataset yields no object and no error", {
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(NULL)), {
    ret <- session$getReturned()
    expect_null(ret$tse())
    expect_null(ret$assay())
    expect_null(ret$audit())
    expect_false(is.null(output$summary_box))
    expect_false(is.null(output$tbl_audit))
    expect_false(is.null(output$plt_length))
    expect_false(is.null(output$plt_quality))
    expect_false(is.null(output$plt_completeness))
    expect_false(is.null(output$plt_score))
    expect_false(is.null(output$plt_prevalence))
  })
})

test_that("default settings keep every contig", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"),
                      min_completeness = 0, completeness_keep_na = TRUE,
                      min_genomad = 0, genomad_keep_na = TRUE,
                      min_virsorter2 = 0, virsorter2_keep_na = TRUE,
                      min_prevalence = 0, min_total_abundance = 0)
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), nrow(tse))
    expect_equal(ret$assay(), "counts")
  })
})

test_that("filters that have no annotation are reported as not applied", {
  # Regression: reading a missing rowData column returned NULL, and comparing
  # NULL to a threshold gives logical(0), which silently removed every contig
  # while looking like a legitimate empty result.
  tse <- vpf_filter_tse()
  SummarizedExperiment::rowData(tse) <-
    SummarizedExperiment::rowData(tse)[, c("Family"), drop = FALSE]
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_prevalence = 0,
                      min_total_abundance = 0)
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), nrow(tse))
    audit <- ret$audit()
    expect_true(all(audit$Applied[audit$Step %in%
      c("Contig length", "CheckV quality", "CheckV completeness",
        "geNomad score", "VirSorter2 score")] == "no"))
    expect_match(audit$Detail[audit$Step == "Contig length"], "absent")
  })
})

test_that("a missing covered-fraction assay does not crash the filter", {
  # Regression: assay(tse, "covfrac") was called unconditionally.
  tse <- vpf_filter_tse()
  tse <- tse[, , drop = FALSE]
  SummarizedExperiment::assays(tse) <-
    SummarizedExperiment::assays(tse)["counts"]
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_length = 0, min_prevalence = 0,
                      min_total_abundance = 0,
                      quality = c("Complete", "Low-quality"))
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), nrow(tse))
    audit <- ret$audit()
    expect_equal(audit$Applied[audit$Step == "Coverage masking"], "no")
  })
})

test_that("an empty quality selection keeps every tier", {
  # Regression: an empty multi-select made `%in%` return all FALSE, so the AND
  # branch removed every contig.
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = character(0), min_prevalence = 0,
                      min_total_abundance = 0)
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), nrow(tse))
    audit <- ret$audit()
    expect_match(audit$Detail[audit$Step == "CheckV quality"], "all kept")
  })
})

test_that("missing inputs fall back to permissive defaults", {
  # Regression: switch(NULL, ...) raised "EXPR must be a length 1 vector", and a
  # NULL threshold made the length comparison return logical(0).
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    ret <- session$getReturned()
    expect_equal(ret$assay(), "counts")
    expect_equal(nrow(ret$tse()), nrow(tse))
  })
})

test_that("NA annotations do not break subsetting and are kept or dropped on request", {
  # Regression: an NA in a logical subscript is an error for
  # SummarizedExperiment; checkv_completeness has one NA in this fixture.
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"),
                      min_completeness = 50, completeness_keep_na = TRUE,
                      min_prevalence = 0, min_total_abundance = 0)
    with_na <- nrow(session$getReturned()$tse())
    session$setInputs(completeness_keep_na = FALSE)
    without_na <- nrow(session$getReturned()$tse())
    expect_equal(with_na - without_na, 1)
  })
})

test_that("the length threshold is inclusive", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 3000,
                      quality = c("Complete", "Low-quality"),
                      min_prevalence = 0, min_total_abundance = 0)
    kept <- session$getReturned()$tse()
    lens <- SummarizedExperiment::rowData(kept)$checkv_contig_length
    expect_true(all(lens >= 3000))
    expect_true(3000 %in% lens)
  })
})

test_that("coverage masking zeroes abundance without removing contigs", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"),
                      min_prevalence = 0, min_total_abundance = 0)
    unmasked <- sum(SummarizedExperiment::assay(session$getReturned()$tse(), "counts"))
    session$setInputs(min_covfrac = 0.9)
    masked_tse <- session$getReturned()$tse()
    expect_equal(nrow(masked_tse), nrow(tse))
    expect_lt(sum(SummarizedExperiment::assay(masked_tse, "counts")), unmasked)
  })
})

test_that("zero surviving contigs is handled everywhere", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 1e9,
                      quality = c("Complete", "Low-quality"),
                      min_prevalence = 0, min_total_abundance = 0)
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), 0)
    expect_false(is.null(output$summary_box))
    expect_false(is.null(output$tbl_audit))
    expect_false(is.null(output$plt_length))
    expect_false(is.null(output$plt_quality))
    expect_false(is.null(output$plt_prevalence))
  })
})

test_that("prevalence and total-abundance filters work on the masked assay", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0.5, min_length = 0,
                      quality = c("Complete", "Low-quality"),
                      min_prevalence = 1, min_total_abundance = 0)
    ret <- session$getReturned()
    audit <- ret$audit()
    expect_equal(audit$Applied[audit$Step == "Prevalence"], "yes")
    kept <- ret$tse()
    if (nrow(kept) > 0) {
      expect_true(all(vpf_prevalence(kept, "counts") >= 1))
    }
  })
})

test_that("the audit table accounts for every contig", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0.3, min_length = 3000,
                      quality = "Complete", min_completeness = 20,
                      completeness_keep_na = FALSE, min_genomad = 0.2,
                      genomad_keep_na = FALSE, min_virsorter2 = 0.3,
                      virsorter2_keep_na = FALSE, min_prevalence = 0,
                      min_total_abundance = 0)
    ret <- session$getReturned()
    audit <- ret$audit()
    expect_equal(audit$Remaining[[1]], nrow(tse))
    expect_equal(audit$Remaining[[nrow(audit)]], nrow(ret$tse()))
    expect_equal(sum(audit$Removed), nrow(tse) - nrow(ret$tse()))
  })
})

test_that("reset sends permissive values back to every control", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 5000,
                      quality = "Complete", min_prevalence = 0,
                      min_total_abundance = 0)
    expect_lt(nrow(session$getReturned()$tse()), nrow(tse))
    # A mock session neither records nor replays update messages, so what can
    # be asserted here is that the observer runs cleanly against a dataset
    # whose optional filters are only partly present. The effect on the
    # controls is covered by the browser test in test-app.R.
    expect_no_error(session$setInputs(reset = 1))
    session$flushReact()
    expect_s4_class(session$getReturned()$tse(), "TreeSummarizedExperiment")

    # The same observer must survive a dataset with no optional filters at all.
    expect_no_error(session$setInputs(reset = 2))
  })
})

test_that("covered fraction is not offered as an abundance assay", {
  tse <- vpf_filter_tse()
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    # Even an explicit request for the detection assay falls back to a valid one.
    session$setInputs(assay = "covfrac", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"),
                      min_prevalence = 0, min_total_abundance = 0)
    expect_equal(session$getReturned()$assay(), "counts")
  })
})

test_that("real ViroProfiler output filters without error", {
  path <- vpf_real_tse_path()
  skip_if(is.na(path), "real ViroProfiler output not available")
  tse <- readRDS(path)
  testServer(mod_vpfilter_server,
             args = list(r_input = vpf_fake_input(tse, "viroprofiler_output.rds")), {
    session$setInputs(assay = "counts", min_covfrac = 0.5, min_length = 0,
                      quality = c("Complete", "High-quality", "Medium-quality",
                                  "Low-quality", "Not-determined"),
                      min_completeness = 0, completeness_keep_na = TRUE,
                      min_genomad = 0, genomad_keep_na = TRUE,
                      min_virsorter2 = 0, virsorter2_keep_na = TRUE,
                      min_prevalence = 0, min_total_abundance = 0)
    ret <- session$getReturned()
    expect_equal(nrow(ret$tse()), 18)
    expect_equal(ncol(ret$tse()), 2)
    for (o in c("summary_box", "tbl_audit", "plt_length", "plt_quality",
                "plt_completeness", "plt_score", "plt_prevalence")) {
      expect_false(is.null(output[[o]]), info = o)
    }
    # This object stores the CoverM trimmed mean under the misleading name.
    session$setInputs(assay = "tmm")
    expect_equal(ret$assay(), "tmm")
    expect_false(is.null(output$assay_description))
  })
})

# ---------------------------------------------------------------------------
# Viral-evidence votes
# ---------------------------------------------------------------------------

vpf_voted_tse <- function() {
  tse <- vpf_filter_tse(n_row = 8, n_col = 4)
  rd <- SummarizedExperiment::rowData(tse)
  rd$viral_vote_taxonomy <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  rd$viral_vote_checkv <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)
  rd$viral_vote_genomad <- c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE)
  rd$viral_vote_vibrant <- c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
  rd$viral_vote_n <- as.integer(rd$viral_vote_taxonomy + rd$viral_vote_checkv +
                                  rd$viral_vote_genomad + rd$viral_vote_vibrant)
  rd$viral_selected <- rd$viral_vote_n > 0
  SummarizedExperiment::rowData(tse) <- rd
  S4Vectors::metadata(tse)$viral_selection <- list(
    rule = "vote",
    combination = "OR (permissive union, not a consensus)",
    votes_used = c("taxonomy", "checkv", "genomad", "vibrant"),
    votes_missing = c("virsorter2", "dvf"),
    n_total = 8L, n_selected = 7L
  )
  tse
}

test_that("viral-evidence votes are exposed as a filter", {
  tse <- vpf_voted_tse()
  expect_equal(vpf_viral_votes(tse), c(4L, 2L, 2L, 2L, 1L, 1L, 1L, 0L))
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"), min_prevalence = 0,
                      min_total_abundance = 0, min_votes = 1)
    expect_false(is.null(filterable()$votes))
    expect_equal(nrow(session$getReturned()$tse()), 8)

    # Excluding contigs admitted on a single detector is the point of the filter.
    session$setInputs(min_votes = 2)
    expect_equal(nrow(session$getReturned()$tse()), 4)
    audit <- session$getReturned()$audit()
    expect_equal(audit$Applied[audit$Step == "Viral-evidence votes"], "yes")

    session$setInputs(min_votes = 4)
    expect_equal(nrow(session$getReturned()$tse()), 1)

    expect_false(is.null(output$plt_votes))
    expect_false(is.null(output$plt_vote_sources))
    expect_false(is.null(output$evidence_status))
  })
})

test_that("a dataset with no vote columns says which step records them", {
  tse <- vpf_filter_tse()
  expect_null(vpf_viral_votes(tse))
  testServer(mod_vpfilter_server, args = list(r_input = vpf_fake_input(tse)), {
    session$setInputs(assay = "counts", min_covfrac = 0, min_length = 0,
                      quality = c("Complete", "Low-quality"), min_prevalence = 0,
                      min_total_abundance = 0)
    expect_null(filterable()$votes)
    audit <- session$getReturned()$audit()
    expect_equal(audit$Applied[audit$Step == "Viral-evidence votes"], "no")
    expect_match(audit$Detail[audit$Step == "Viral-evidence votes"],
                 "annotate_viral_votes")
    expect_false(is.null(output$evidence_status))
    expect_false(is.null(output$plt_votes))
    expect_false(is.null(output$plt_vote_sources))
  })
})

test_that("legacy tmm assays are renamed on read", {
  # Regression: CoverM's trimmed mean was stored under the name `tmm`, which
  # reads as edgeR's unrelated TMM normalization.
  tse <- vpf_filter_tse()
  names(SummarizedExperiment::assays(tse))[1] <- "tmm"
  f <- tempfile(fileext = ".rds")
  saveRDS(tse, f)
  on.exit(unlink(f), add = TRUE)
  res <- vpf_read_tse(f)
  expect_null(res$error)
  expect_true("trimmed_mean" %in% SummarizedExperiment::assayNames(res$tse))
  expect_false("tmm" %in% SummarizedExperiment::assayNames(res$tse))
  expect_equal(unname(res$renamed_assays), "trimmed_mean")
})

test_that("real viral-evidence votes drive the filter", {
  p <- Sys.getenv("VPFKIT_ENRICHED_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output_enriched.rds"
  }
  skip_if(!file.exists(p), "enriched ViroProfiler output not available")
  tse <- readRDS(p)
  votes <- vpf_viral_votes(tse)
  expect_false(is.null(votes))
  testServer(mod_vpfilter_server,
             args = list(r_input = vpf_fake_input(tse, "enriched.rds")), {
    session$setInputs(assay = "trimmed_mean", min_covfrac = 0.5, min_length = 0,
                      quality = c("Complete", "High-quality", "Medium-quality",
                                  "Low-quality", "Not-determined"),
                      min_completeness = 0, completeness_keep_na = TRUE,
                      min_genomad = 0, genomad_keep_na = TRUE,
                      min_virsorter2 = 0, virsorter2_keep_na = TRUE,
                      min_prevalence = 0, min_total_abundance = 0, min_votes = 1)
    ret <- session$getReturned()
    expect_equal(ret$assay(), "trimmed_mean")
    expect_equal(nrow(ret$tse()), 18)

    # One contig entered on VIBRANT alone; requiring two detectors drops it.
    session$setInputs(min_votes = 2)
    expect_equal(nrow(ret$tse()), sum(votes >= 2))
    expect_lt(nrow(ret$tse()), 18)
    session$setInputs(min_votes = 4)
    expect_equal(nrow(ret$tse()), sum(votes >= 4))

    for (o in c("evidence_status", "plt_votes", "plt_vote_sources",
                "summary_box", "tbl_audit", "plt_prevalence")) {
      expect_false(is.null(output[[o]]), info = o)
    }
  })
})
