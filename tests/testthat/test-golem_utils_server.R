test_that("not_in works", {
  expect_true(1 %not_in% 2:10)
  expect_false(1 %not_in% 1:10)
})

test_that("not_null works", {
  expect_true(not_null(1))
  expect_false(not_null(NULL))
})

test_that("not_na works", {
  expect_true(not_na(1))
  expect_false(not_na(NA))
})

test_that("drop_nulls works", {
  expect_equal(
    drop_nulls(
      list(x = NULL, y = 2)
    ),
    list(y = 2)
  )
})

test_that("%||% works", {
  expect_equal(
    NULL %||% 1,
    1
  )
  expect_equal(
    2 %||% 1,
    2
  )
})

test_that("%|NA|% works", {
  expect_equal(
    NA %|NA|% 1,
    1
  )
  expect_equal(
    2 %|NA|% 1,
    2
  )
})

test_that("rv and rvtl work", {
  expect_true(
    inherits(rv, "function")
  )
  expect_true(
    inherits(rvtl, "function")
  )
})


# ---------------------------------------------------------------------------
# Viewer helpers
# ---------------------------------------------------------------------------

vpf_fixture <- function(name) {
  paths <- vpf_demo_datasets()
  hit <- paths[grepl(name, names(paths), ignore.case = TRUE)]
  if (length(hit) == 0) {
    return(NULL)
  }
  readRDS(hit[[1]])
}

vpf_toy_tse <- function(n_row = 6, n_col = 4, seed = 1) {
  set.seed(seed)
  counts <- matrix(as.integer(stats::rpois(n_row * n_col, 20)), nrow = n_row)
  rownames(counts) <- paste0("contig_", seq_len(n_row))
  colnames(counts) <- paste0("S", seq_len(n_col))
  covfrac <- matrix(stats::runif(n_row * n_col), nrow = n_row,
                    dimnames = dimnames(counts))
  TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = counts, covfrac = covfrac),
    rowData = S4Vectors::DataFrame(
      checkv_contig_length = seq(1000, by = 1000, length.out = n_row),
      checkv_quality = factor(rep(c("Complete", "Low-quality"), length.out = n_row)),
      checkv_completeness = seq(0, 100, length.out = n_row),
      genomad_score = seq(0.1, 0.9, length.out = n_row),
      Family = c("Siphoviridae", "Myoviridae", rep(NA_character_, n_row - 2)),
      row.names = rownames(counts)
    ),
    colData = S4Vectors::DataFrame(
      sample_name = colnames(counts),
      group = rep(c("A", "B"), length.out = n_col),
      row.names = colnames(counts)
    )
  )
}

test_that("vpf_blank_na turns pipeline placeholders into NA", {
  expect_equal(
    vpf_blank_na(c("Caudoviricetes", "", " ", "NA", "unknown", "-")),
    c("Caudoviricetes", NA, NA, NA, NA, NA)
  )
  expect_equal(vpf_blank_na(factor(c("a", ""))), c("a", NA))
})

test_that("assay dictionary separates abundance from detection", {
  tse <- vpf_toy_tse()
  info <- vpf_assay_info(tse)
  expect_setequal(info$assay, c("counts", "covfrac"))
  expect_equal(info$role[info$assay == "covfrac"], "detection")
  # covered fraction must never be offered as an abundance: it is a breadth
  # measure and using it for composition or diversity is a silent error.
  expect_false("covfrac" %in% vpf_abundance_assay_choices(tse))
  expect_equal(vpf_default_assay(tse), "counts")
  expect_equal(vpf_covfrac_assay(tse), "covfrac")
})

test_that("both spellings of the CoverM trimmed mean are recognized", {
  dict <- vpf_assay_dictionary()
  expect_equal(dict$tmm$role, "abundance")
  expect_equal(dict$trimmed_mean$role, "abundance")
  expect_match(dict$tmm$description, "NOT edgeR's TMM", fixed = TRUE)
  expect_match(dict$trimmed_mean$description, "NOT edgeR's TMM", fixed = TRUE)
  expect_false(dict$tmm$integer_counts)
  expect_false(dict$trimmed_mean$integer_counts)
})

test_that("unknown assay names degrade instead of failing", {
  tse <- vpf_toy_tse()
  SummarizedExperiment::assay(tse, "mystery") <-
    SummarizedExperiment::assay(tse, "counts")
  info <- vpf_assay_info(tse)
  expect_false(info$known[info$assay == "mystery"])
  expect_true("mystery" %in% vpf_abundance_assay_choices(tse))
})

test_that("vpf_is_count_assay only accepts raw integer counts", {
  tse <- vpf_toy_tse()
  expect_true(vpf_is_count_assay(tse, "counts"))
  expect_false(vpf_is_count_assay(tse, "covfrac"))
  expect_false(vpf_is_count_assay(tse, "absent"))
  expect_false(vpf_is_count_assay(NULL, "counts"))
})

test_that("sample identifiers prefer sample_name", {
  tse <- vpf_toy_tse()
  expect_equal(vpf_sample_names(tse), c("S1", "S2", "S3", "S4"))
  expect_equal(vpf_sample_names(NULL), character(0))
})

test_that("available ranks ignore all-blank columns", {
  tse <- vpf_toy_tse()
  expect_equal(vpf_available_ranks(tse), "Family")
  SummarizedExperiment::rowData(tse)$Genus <- rep("", nrow(tse))
  expect_false("Genus" %in% vpf_available_ranks(tse))
})

test_that("a sample identifier is never a grouping variable", {
  tse <- vpf_toy_tse()
  expect_equal(vpf_group_candidates(tse), "group")
  # All-unique columns would create one group per sample.
  SummarizedExperiment::colData(tse)$replicate <- letters[seq_len(ncol(tse))]
  expect_false("replicate" %in% vpf_group_candidates(tse))
  expect_false("sample_name" %in% vpf_group_candidates(tse))
  # Recognized by name even when it happens to repeat values.
  SummarizedExperiment::colData(tse)$library_id <- rep(c("L1", "L2"), 2)
  expect_false("library_id" %in% vpf_group_candidates(tse))
})

test_that("a two-sample run may still carry a labelling variable", {
  # With n = 2 any two-level variable has one sample per level, so refusing it
  # outright would make an uploaded metadata table useless on a pilot run. The
  # design check is what refuses the tests.
  tse <- vpf_toy_tse(n_col = 2)
  SummarizedExperiment::colData(tse)$condition <- c("healthy", "colitis")
  expect_true("condition" %in% vpf_group_candidates(tse))
  d <- vpf_design_check(tse, "condition", "counts")
  expect_equal(d$n_groups, 2L)
  expect_equal(unname(d$min_group_size), 1L)
  expect_false(d$can_group_test)
  expect_false(d$can_permanova)
  expect_false(d$can_differential)
})

test_that("permutation counts follow the multinomial formula", {
  # 4 vs 4 has choose(8,4) = 70 labelled allocations; equally sized groups are
  # interchangeable so 35 partitions remain.
  expect_equal(vpf_permutation_count(c(4, 4)), 35)
  expect_equal(vpf_permutation_count(c(3, 4)), 35)
  expect_equal(vpf_permutation_count(c(2, 2)), 3)
  expect_equal(vpf_permutation_count(5), 1)
})

test_that("design guards match the mathematics of each procedure", {
  # PCoA needs n >= 3 because classical MDS has at most n - 1 axes.
  two <- vpf_toy_tse(n_col = 2)
  d2 <- vpf_design_check(two, NULL, "counts")
  expect_equal(d2$n_samples, 2L)
  expect_true(d2$can_alpha)
  expect_true(d2$can_pairwise_beta)
  expect_false(d2$can_pcoa)
  expect_false(d2$can_nmds)
  expect_false(d2$can_permanova)
  expect_false(d2$can_group_test)

  three <- vpf_toy_tse(n_col = 3)
  d3 <- vpf_design_check(three, NULL, "counts")
  expect_true(d3$can_pcoa)
  expect_true(d3$pcoa_saturated)
  expect_false(d3$can_nmds)

  # NMDS in 2D is saturated below n = 6, where stress is trivially near zero.
  six <- vpf_toy_tse(n_col = 6)
  expect_false(vpf_design_check(vpf_toy_tse(n_col = 5), NULL, "counts")$can_nmds)
  expect_true(vpf_design_check(six, NULL, "counts")$can_nmds)

  d_grp <- vpf_design_check(six, "group", "counts")
  expect_equal(d_grp$n_groups, 2L)
  expect_equal(unname(d_grp$min_group_size), 3L)
  expect_true(d_grp$can_group_test)
  expect_true(d_grp$can_permanova)
  # 3 vs 3 gives 10 distinguishable permutations, so p < 0.05 is unreachable.
  expect_false(d_grp$permanova_resolution_ok)
  # 4 vs 4 gives 35, which is still below the 20-permutation floor? No: 35 > 20.
  eight <- vpf_toy_tse(n_col = 8)
  expect_true(vpf_design_check(eight, "group", "counts")$permanova_resolution_ok)
})

test_that("a single-level grouping variable cannot pass the guards", {
  tse <- vpf_toy_tse(n_col = 6)
  SummarizedExperiment::colData(tse)$group <- rep("only", 6)
  d <- vpf_design_check(tse, "group", "counts")
  expect_equal(d$n_groups, 1L)
  expect_false(d$can_group_test)
  expect_false(d$can_permanova)
  expect_false(d$can_differential)
})

test_that("empty libraries are detected and excluded", {
  tse <- vpf_toy_tse(n_col = 4)
  SummarizedExperiment::assay(tse, "counts")[, 2] <- 0L
  d <- vpf_design_check(tse, NULL, "counts")
  expect_equal(d$empty_libraries, "S2")
})

test_that("coverage masking zeroes low-breadth abundances and treats NA as absent", {
  tse <- vpf_toy_tse(n_row = 3, n_col = 2)
  SummarizedExperiment::assay(tse, "counts")[] <- 10L
  cf <- matrix(c(0.9, 0.1, NA, 0.6, 0.4, 0.8), nrow = 3)
  dimnames(cf) <- dimnames(SummarizedExperiment::assay(tse, "counts"))
  SummarizedExperiment::assay(tse, "covfrac") <- cf
  out <- refind_abundance(tse, "counts", "covfrac", 0.5)
  expect_equal(
    unname(as.numeric(SummarizedExperiment::assay(out, "counts"))),
    c(10, 0, 0, 10, 0, 10)
  )
  # A missing assay is a no-op rather than an error.
  expect_warning(expect_identical(refind_abundance(tse, "counts", "nope", 0.5), tse))
  expect_error(refind_abundance(tse, "counts", "covfrac", NA_real_))
})

test_that("row subsetting resolves NA in the subscript", {
  tse <- vpf_toy_tse(n_row = 4)
  keep <- c(TRUE, NA, FALSE, TRUE)
  expect_equal(nrow(vpf_subset_rows(tse, keep, na_keeps = FALSE)), 2)
  expect_equal(nrow(vpf_subset_rows(tse, keep, na_keeps = TRUE)), 3)
  # A length mismatch would otherwise recycle silently.
  expect_equal(nrow(vpf_subset_rows(tse, c(TRUE, FALSE))), 4)
})

test_that("numeric coercion refuses columns that cannot be compared", {
  expect_equal(vpf_as_numeric_column(c(1, 2)), c(1, 2))
  expect_equal(vpf_as_numeric_column(factor(c("10", "2"))), c(10, 2))
  expect_null(vpf_as_numeric_column(c("high", "low")))
  expect_null(vpf_as_numeric_column(rep(NA, 3)))
  expect_null(vpf_as_numeric_column(NULL))
})

test_that("prevalence counts detections, not abundance", {
  tse <- vpf_toy_tse(n_row = 2, n_col = 4)
  SummarizedExperiment::assay(tse, "counts")[1, ] <- c(5L, 0L, 3L, 0L)
  SummarizedExperiment::assay(tse, "counts")[2, ] <- c(0L, 0L, 0L, 0L)
  expect_equal(unname(vpf_prevalence(tse, "counts")), c(0.5, 0))
  expect_equal(vpf_prevalence(tse, "absent"), numeric(0))
})

test_that("vpf_read_tse reports what is wrong instead of failing", {
  bad <- tempfile(fileext = ".rds")
  saveRDS(list(a = 1), bad)
  on.exit(unlink(bad), add = TRUE)
  expect_match(vpf_read_tse(bad)$error, "TreeSummarizedExperiment")
  expect_match(vpf_read_tse("/nonexistent/file.rds")$error, "not found")
  expect_match(vpf_read_tse(NULL)$error, "No file selected")

  good <- tempfile(fileext = ".rds")
  saveRDS(vpf_toy_tse(), good)
  on.exit(unlink(good), add = TRUE)
  res <- vpf_read_tse(good)
  expect_null(res$error)
  expect_s4_class(res$tse, "TreeSummarizedExperiment")
})

test_that("sample metadata joins by identifier and reports mismatches", {
  tse <- vpf_toy_tse(n_col = 4)
  meta <- data.frame(
    sample = c("S1", "S2", "S3", "S9"),
    treatment = c("ctrl", "ctrl", "drug", "drug"),
    stringsAsFactors = FALSE
  )
  res <- vpf_join_sample_metadata(tse, meta)
  expect_null(res$error)
  expect_equal(res$matched, 3L)
  expect_equal(res$unmatched_samples, "S4")
  expect_equal(res$unused_rows, "S9")
  expect_true("treatment" %in% colnames(SummarizedExperiment::colData(res$tse)))
  expect_true("treatment" %in% vpf_group_candidates(res$tse))

  # A clashing column name must not overwrite existing colData.
  meta2 <- data.frame(sample = c("S1", "S2", "S3", "S4"),
                      group = c("x", "x", "y", "y"), stringsAsFactors = FALSE)
  res2 <- vpf_join_sample_metadata(tse, meta2)
  expect_true("group_meta" %in% colnames(SummarizedExperiment::colData(res2$tse)))
  expect_true("group" %in% colnames(SummarizedExperiment::colData(res2$tse)))

  # No overlap at all is an explicit, diagnosable failure.
  meta3 <- data.frame(sample = c("X1", "X2"), g = c("a", "b"), stringsAsFactors = FALSE)
  res3 <- vpf_join_sample_metadata(tse, meta3)
  expect_match(res3$error, "None of the sample identifiers matched")
})

test_that("metadata files are parsed from csv, tsv and xlsx", {
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = c("S1", "S2"), g = c("a", "b")),
                   csv, row.names = FALSE)
  on.exit(unlink(csv), add = TRUE)
  expect_null(vpf_read_metadata_file(csv, "meta.csv")$error)

  tsv <- tempfile(fileext = ".tsv")
  utils::write.table(data.frame(sample = c("S1", "S2"), g = c("a", "b")),
                     tsv, sep = "\t", row.names = FALSE, quote = FALSE)
  on.exit(unlink(tsv), add = TRUE)
  parsed <- vpf_read_metadata_file(tsv, "meta.tsv")
  expect_null(parsed$error)
  expect_equal(colnames(parsed$data), c("sample", "g"))

  one_col <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(sample = c("S1", "S2")), one_col, row.names = FALSE)
  on.exit(unlink(one_col), add = TRUE)
  expect_match(vpf_read_metadata_file(one_col, "one.csv")$error,
               "at least one metadata column")
})

test_that("annotation status names the option that produces each family", {
  tse <- vpf_toy_tse()
  st <- vpf_annotation_status(tse)
  expect_true(st$available[st$key == "checkv"])
  expect_true(st$available[st$key == "genomad"])
  expect_false(st$available[st$key == "iphop"])
  expect_match(st$param[st$key == "iphop"], "--use_iphop")
  expect_match(st$param[st$key == "replicyc"], "--replicyc")
  expect_false(st$available[st$key == "gene_annotations"])
})

test_that("gene annotations are read from metadata when present", {
  tse <- vpf_toy_tse()
  expect_null(vpf_gene_annotations(tse))
  S4Vectors::metadata(tse)$gene_annotations <- data.frame(
    Contig = "contig_1", gene_id = "g1", stringsAsFactors = FALSE
  )
  expect_equal(nrow(vpf_gene_annotations(tse)), 1)
  S4Vectors::metadata(tse)$gene_annotations <- data.frame()
  expect_null(vpf_gene_annotations(tse))
})

test_that("message plots are real plotly objects", {
  p <- vpf_message_plot("nothing to show")
  expect_s3_class(p, "plotly")
  broken <- vpf_ggplotly("not a ggplot")
  expect_s3_class(broken, "plotly")
})

test_that("bundled demo datasets are discoverable", {
  demos <- vpf_demo_datasets()
  skip_if(length(demos) == 0, "demo datasets not present in this checkout")
  expect_true(all(file.exists(unlist(demos))))
  tse <- readRDS(demos[[1]])
  expect_s4_class(tse, "TreeSummarizedExperiment")
})
