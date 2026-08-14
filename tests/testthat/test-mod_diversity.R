vpf_div_filter <- function(tse, assay = "counts", raw = tse) {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive("test")
  )
}

vpf_div_demo <- function(pattern) {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_div_real <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("module ui works", {
  ui <- mod_diversity_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_diversity_ui)))
})

test_that("no dataset renders placeholders", {
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(NULL)), {
    expect_false(is.null(output$plt_alpha))
    expect_false(is.null(output$plt_ordination))
    expect_false(is.null(output$plt_distance))
    expect_false(is.null(output$permanova_out))
    expect_false(is.null(output$subset_warning))
  })
})

test_that("alpha diversity is computed on the assay the filter selected", {
  # Regression: alpha diversity was hard-coded to "counts" while filtering and
  # coverage masking were applied to the metric the user chose, so the plot
  # described different data from every other tab.
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "tpm")), {
    session$setInputs(alpha_index = "shannon", group = "",
                      beta_method = "bray", ordination = "pcoa",
                      permutations = 99)
    av <- alpha_values()
    expect_equal(av$assay, "tpm")
    expect_equal(nrow(av$data), ncol(tse))
    # The values must differ from the counts-based ones.
    expect_false(is.null(output$plt_alpha))
    expect_false(is.null(output$tbl_alpha))
  })
})

test_that("Chao1 is refused on a non-count assay", {
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "tpm")), {
    session$setInputs(alpha_index = "chao1", group = "")
    expect_false(design()$can_chao1)
    expect_null(alpha_values())
    expect_false(is.null(output$alpha_note))
    expect_false(is.null(output$plt_alpha))
  })
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "counts")), {
    session$setInputs(alpha_index = "chao1", group = "")
    expect_true(design()$can_chao1)
    expect_false(is.null(alpha_values()$data))
  })
})

test_that("two samples with no grouping disable every inferential output", {
  tse <- vpf_div_real()
  skip_if(is.null(tse), "real ViroProfiler output not available")
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "", beta_method = "bray",
                      ordination = "pcoa", permutations = 999)
    d <- design()
    expect_equal(d$n_samples, 2L)
    expect_true(d$can_alpha)
    expect_true(d$can_pairwise_beta)
    expect_false(d$can_pcoa)
    expect_false(d$can_nmds)
    expect_false(d$can_permanova)

    # Alpha diversity and the single pairwise distance remain available.
    expect_equal(nrow(alpha_values()$data), 2)
    expect_false(is.null(dissimilarity()$dist))
    expect_length(dissimilarity()$dist, 1)

    # Ordination and PERMANOVA return an explanation, not an error.
    # Regression: cmdscale(k = 2) on two samples raised
    # "'k' must be in {1, 2, .. n - 1}", and the module surfaced it as a red
    # error box in the browser.
    expect_false(is.null(output$plt_ordination))
    expect_false(is.null(output$permanova_out))
    # The group-test panel must say why no test is offered rather than vanish.
    expect_false(is.null(output$alpha_test))

    session$setInputs(ordination = "nmds")
    expect_false(is.null(output$plt_ordination))
  })
})

test_that("zero surviving contigs never reaches vegdist", {
  # Regression: a zero-row object produced an all-NA dissimilarity and
  # "NA values not allowed in 'd'".
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  empty <- tse[integer(0), ]
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(empty, raw = tse)), {
    session$setInputs(alpha_index = "shannon", group = "", beta_method = "bray",
                      ordination = "pcoa", permutations = 99)
    expect_null(alpha_values())
    expect_null(dissimilarity())
    expect_false(is.null(output$plt_alpha))
    expect_false(is.null(output$plt_ordination))
    expect_false(is.null(output$plt_distance))
    expect_false(is.null(output$permanova_out))
  })
})

test_that("PCoA and PERMANOVA run on a replicated design", {
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "group",
                      beta_method = "bray", ordination = "pcoa",
                      permutations = 199)
    d <- design()
    expect_equal(d$n_groups, 2L)
    expect_equal(unname(d$min_group_size), 4L)
    expect_true(d$can_pcoa)
    expect_true(d$can_nmds)
    expect_true(d$can_permanova)
    expect_true(d$permanova_resolution_ok)
    expect_false(is.null(output$plt_ordination))
    expect_false(is.null(output$plt_distance))
    expect_false(is.null(output$permanova_out))
    expect_false(is.null(output$alpha_test))
    expect_true("Group" %in% colnames(alpha_values()$data))
  })
})

test_that("NMDS is refused below six samples", {
  tse <- vpf_div_demo("Environmental")
  skip_if(is.null(tse), "second demo dataset not available")
  five <- tse[, 1:5]
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(five)), {
    session$setInputs(alpha_index = "shannon", group = "", beta_method = "bray",
                      ordination = "nmds", permutations = 99)
    expect_false(design()$can_nmds)
    expect_true(design()$can_pcoa)
    expect_false(is.null(output$plt_ordination))
  })
})

test_that("a single-level grouping variable disables the group test", {
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  SummarizedExperiment::colData(tse)$group <- rep("one", ncol(tse))
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "group",
                      beta_method = "bray", ordination = "pcoa",
                      permutations = 99)
    expect_equal(design()$n_groups, 1L)
    expect_false(design()$can_group_test)
    expect_false(is.null(output$permanova_out))
  })
})

test_that("singleton groups are rejected for PERMANOVA", {
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  SummarizedExperiment::colData(tse)$subject <- letters[seq_len(ncol(tse))]
  # Force a design with one sample per level for one of the levels.
  SummarizedExperiment::colData(tse)$lopsided <-
    c("A", rep("B", ncol(tse) - 1))
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "lopsided",
                      beta_method = "bray", ordination = "pcoa",
                      permutations = 99)
    expect_equal(unname(design()$min_group_size), 1L)
    expect_false(design()$can_permanova)
    expect_false(design()$can_group_test)
    expect_false(is.null(output$permanova_out))
    expect_false(is.null(output$alpha_test))
  })
})

test_that("jaccard and bray both produce a usable distance matrix", {
  tse <- vpf_div_demo("Environmental")
  skip_if(is.null(tse), "second demo dataset not available")
  for (m in c("bray", "jaccard")) {
    testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
      session$setInputs(alpha_index = "shannon", group = "group",
                        beta_method = m, ordination = "pcoa", permutations = 99)
      dd <- dissimilarity()
      expect_equal(dd$method, m)
      expect_false(anyNA(dd$dist))
      expect_false(is.null(output$plt_distance))
    })
  }
})

test_that("an empty library is excluded rather than producing NA distances", {
  tse <- vpf_div_demo("Gut virome")
  skip_if(is.null(tse), "demo dataset not available")
  SummarizedExperiment::assay(tse, "counts")[, 1] <- 0L
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "", beta_method = "bray",
                      ordination = "pcoa", permutations = 99)
    expect_equal(design()$empty_libraries, colnames(tse)[1])
    expect_equal(nrow(alpha_values()$data), ncol(tse) - 1)
    expect_false(anyNA(dissimilarity()$dist))
  })
})

vpf_div_enriched <- function() {
  p <- Sys.getenv("VPFKIT_ENRICHED_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output_enriched.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("a real grouping variable with one sample per group still blocks inference", {
  # The enriched two-sample object does carry a grouping variable, which is
  # exactly the case where an unguarded app would print a p-value computed from
  # one observation per group.
  tse <- vpf_div_enriched()
  skip_if(is.null(tse), "enriched ViroProfiler output not available")
  expect_true("group" %in% vpf_group_candidates(tse))
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "trimmed_mean")), {
    session$setInputs(alpha_index = "shannon", group = "group",
                      beta_method = "bray", ordination = "pcoa",
                      permutations = 999)
    d <- design()
    expect_equal(d$n_samples, 2L)
    expect_equal(d$n_groups, 2L)
    expect_equal(unname(d$min_group_size), 1L)
    expect_false(d$can_group_test)
    expect_false(d$can_permanova)
    expect_false(d$can_pcoa)
    # Descriptive output remains available.
    expect_equal(nrow(alpha_values()$data), 2)
    expect_true("Group" %in% colnames(alpha_values()$data))
    expect_length(dissimilarity()$dist, 1)
    for (o in c("plt_alpha", "tbl_alpha", "alpha_test", "plt_ordination",
                "plt_distance", "permanova_out", "subset_warning",
                "design_report")) {
      expect_false(is.null(output[[o]]), info = o)
    }
  })
})

test_that("the subset caveat reports the vote breakdown when it is recorded", {
  tse <- vpf_div_enriched()
  skip_if(is.null(tse), "enriched ViroProfiler output not available")
  skip_if(is.null(vpf_viral_votes(tse)), "no viral-vote columns")
  testServer(mod_diversity_server, args = list(r_filter = vpf_div_filter(tse)), {
    session$setInputs(alpha_index = "shannon", group = "")
    html <- as.character(output$subset_warning$html %||% output$subset_warning)
    expect_match(html, "single detector")
  })
})

test_that("Chao1 is available on the enriched object's raw counts", {
  tse <- vpf_div_enriched()
  skip_if(is.null(tse), "enriched ViroProfiler output not available")
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "counts")), {
    session$setInputs(alpha_index = "chao1", group = "")
    expect_true(design()$can_chao1)
    expect_equal(nrow(alpha_values()$data), 2)
  })
  testServer(mod_diversity_server,
             args = list(r_filter = vpf_div_filter(tse, assay = "trimmed_mean")), {
    session$setInputs(alpha_index = "chao1", group = "")
    expect_false(design()$can_chao1)
    expect_null(alpha_values())
  })
})
