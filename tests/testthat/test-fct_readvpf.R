# Fixtures in tests/testthat/fixtures/ reproduce the header line of the real
# ViroProfiler output they stand for, verbatim. Where the real header was
# discovered from a reference run, the file it came from is named in the
# comment above the corresponding test.
#
# Reference run:
#   /home/allen/data2/testdata/viroprofiler_real_full  (2 samples, 22 contigs)

fx <- function(...) test_path("fixtures", ...)

REAL_VP <- "/home/allen/data2/testdata/viroprofiler_real_full"

skip_without_real_data <- function() {
  testthat::skip_if_not(dir.exists(REAL_VP),
                        "reference ViroProfiler run is not available")
}

# Assembles the fixture TSE that several tests share.
build_fixture_tse <- function(...) {
  create_vpftse(
    fin_abcount    = fx("coverm_abundance_count.tsv"),
    fin_abtpm      = fx("coverm_abundance_tpm.tsv"),
    fin_abtmm      = fx("coverm_abundance_tmm.tsv"),
    fin_abcov      = fx("coverm_abundance_covfrac.tsv"),
    fin_taxa       = fx("taxonomy_tse.tsv"),
    fin_checkv     = fx("checkv_quality_summary.tsv"),
    fin_virsorter2 = fx("virsorter2_final_viral_score.tsv"),
    fin_vibrant    = fx("vibrant_quality.tsv"),
    fin_dvf        = fx("dvf_virus.tsv"),
    fin_replicyc   = fx("bacphlip.bacphlip"),
    ...
  )
}


# --- read_checkv ---------------------------------------------------------
# Header from checkv/quality_summary.tsv of the reference run.
test_that("read_checkv parses the CheckV quality summary", {
  df <- read_checkv(fx("checkv_quality_summary.tsv"))

  expect_s3_class(df, "data.frame")
  expect_true(all(c("Contig", "checkv_quality", "checkv_contig_length",
                    "checkv_completeness", "checkv_provirus", "checkv_proviral_length",
                    "checkv_kmer_freq", "checkv_warnings") %in% colnames(df)))
  expect_equal(nrow(df), 5)
  expect_true(is.factor(df$checkv_quality))
  expect_equal(levels(df$checkv_quality),
               c("Complete", "High-quality", "Medium-quality", "Low-quality",
                 "Not-determined"))
  # "NA" in the completeness column is a real missing value, not the string.
  expect_true(is.numeric(df$checkv_completeness))
  expect_true(is.na(df$checkv_completeness[df$Contig == "ctg005"]))
})

test_that("read_checkv warns about quality values outside the CheckV vocabulary", {
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("contig_id\tcheckv_quality\tcompleteness",
               "ctg001\tPretty-good\t80"), bad)
  expect_warning(read_checkv(bad), "unrecognized checkv_quality")
})

test_that("read_checkv errors on missing file and missing columns", {
  expect_error(read_checkv("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_checkv(bad), "missing columns")
})

test_that("read_checkv errors on a zero-byte file rather than parsing nothing", {
  empty <- withr::local_tempfile(fileext = ".tsv")
  file.create(empty)
  expect_error(read_checkv(empty), "empty")
})


# --- read_virsorter2 -----------------------------------------------------
# Header from virsorter2/out_vs2/final-viral-score.tsv of the reference run,
# which has one score column because ViroProfiler passes a single
# --include-groups value.
test_that("read_virsorter2 parses the single-group score table", {
  df <- read_virsorter2(fx("virsorter2_final_viral_score.tsv"))
  expect_true(all(c("Contig", "virsorter2_max_score", "virsorter2_max_score_group",
                    "virsorter2_viral", "virsorter2_cellular") %in% colnames(df)))
  expect_equal(nrow(df), 5)
})

test_that("read_virsorter2 does not assume a fixed set of group columns", {
  df <- read_virsorter2(fx("virsorter2_final_viral_score_multigroup.tsv"))
  expect_true("virsorter2_NCLDV" %in% colnames(df))
  expect_equal(nrow(df), 3)
  expect_equal(df$virsorter2_max_score_group[df$Contig == "ctg003"], "NCLDV")
})

test_that("read_virsorter2 strips the || decoration VirSorter2 adds by default", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("seqname\tmax_score\tmax_score_group",
               "ctg001||full\t0.99\tdsDNAphage",
               "ctg002||0_partial\t0.80\tdsDNAphage"), f)
  df <- read_virsorter2(f)
  expect_equal(sort(df$Contig), c("ctg001", "ctg002"))
})

test_that("read_virsorter2 errors on missing file and missing columns", {
  expect_error(read_virsorter2("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_virsorter2(bad), "missing columns")
})


# --- read_vibrant --------------------------------------------------------
# Header from vibrant/.../VIBRANT_genome_quality_contigs.tsv of the reference run.
test_that("read_vibrant parses quality and orders the levels best to worst", {
  df <- suppressWarnings(read_vibrant(fx("vibrant_quality.tsv")))
  expect_true(all(c("Contig", "vibrant_quality", "vibrant_replicyc") %in% colnames(df)))
  expect_true(is.factor(df$vibrant_quality))
  # Not the alphabetical order: medium ranks above low.
  expect_equal(levels(df$vibrant_quality),
               c("complete circular", "high quality draft", "medium quality draft",
                 "low quality draft"))
})

test_that("read_vibrant collapses prophage fragments keeping the best quality", {
  # ctg001 appears three times: itself (low), _fragment_1 (high), _fragment_2 (low).
  # The previous implementation sorted by the string literal "vibrant_quality",
  # which sorted nothing, and kept whichever row came first.
  df <- suppressWarnings(read_vibrant(fx("vibrant_quality.tsv")))
  expect_equal(sum(df$Contig == "ctg001"), 1)
  expect_equal(as.character(df$vibrant_quality[df$Contig == "ctg001"]),
               "high quality draft")
})

test_that("read_vibrant errors on missing file", {
  expect_error(read_vibrant("/nonexistent/file.tsv"), "not found")
})


# --- read_dvf ------------------------------------------------------------
# DeepVirFinder writes `name len score pvalue qvalue`; the same header is used
# by viroprofiler/assets/no_dvf_scores.tsv.
test_that("read_dvf applies its three thresholds and keeps the raw values", {
  df <- read_dvf(fx("dvf_virus.tsv"))
  expect_true(all(c("Contig", "dvf_score", "dvf_pvalue", "dvf_qvalue") %in% colnames(df)))
  expect_equal(nrow(df), 2)
  expect_setequal(df$Contig, c("ctg001", "ctg004"))
})

test_that("read_dvf can return every row unfiltered", {
  df <- read_dvf(fx("dvf_virus.tsv"), filter = FALSE)
  expect_equal(nrow(df), 5)
})

test_that("read_dvf errors on missing file", {
  expect_error(read_dvf("/nonexistent/file.tsv"), "not found")
})


# --- read_coverm ---------------------------------------------------------
# Header from abundance/abundance_contigs_count.tsv.gz of the reference run:
# `Contig` plus one column per sample, named with the real sample names.
test_that("read_coverm parses an abundance table", {
  df <- read_coverm(fx("coverm_abundance_count.tsv"))
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 5)
  expect_equal(colnames(df), c("HT02", "UC20", "HT03"))
  expect_true("ctg001" %in% rownames(df))
})

test_that("read_coverm leaves sample names alone unless asked to rename them", {
  # Earlier versions rewrote any column matching "ds10Ms" to "Sample_",
  # renaming samples in datasets that had nothing to do with that project.
  df <- read_coverm(fx("coverm_abundance_count.tsv"))
  expect_false(any(grepl("^Sample_", colnames(df))))
  df2 <- read_coverm(fx("coverm_abundance_count.tsv"),
                     sample_rename = c("HT" = "Case"))
  expect_true("Case02" %in% colnames(df2))
})

test_that("read_coverm refuses a table whose sample column is not numeric", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("Contig\tS1", "ctg001\tnot_a_number"), f)
  expect_error(read_coverm(f), "non-numeric")
})

test_that("read_coverm errors on missing file and on a missing Contig column", {
  expect_error(read_coverm("/nonexistent/file.tsv"), "not found")
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("Genome\tS1", "ctg001\t5"), f)
  expect_error(read_coverm(f), "missing 'Contig' column")
})


# --- read_taxonomy2 ------------------------------------------------------
# Header from taxonomy/taxonomy_tse.tsv of the reference run, written by
# viroprofiler/bin/merge_taxonomy.py.
test_that("read_taxonomy2 parses the merged taxonomy and drops unassigned rows", {
  df <- read_taxonomy2(fx("taxonomy_tse.tsv"))
  expect_true(all(c("Contig", "Domain", "Family", "Species") %in% colnames(df)))
  # ctg005 has taxa_id = 0, meaning no rank was resolved at all.
  expect_false("ctg005" %in% df$Contig)
  expect_equal(nrow(df), 3)
})

test_that("read_taxonomy2 turns an unresolved rank into NA, not an empty string", {
  # merge_taxonomy.py leaves unresolved ranks empty. An empty string passes
  # every is.na() test, which makes an unassigned rank look assigned.
  df <- read_taxonomy2(fx("taxonomy_tse.tsv"))
  expect_true(is.na(df$Species[df$Contig == "ctg001"]))
  expect_true(is.na(df$Family[df$Contig == "ctg002"]))
  expect_equal(df$Domain[df$Contig == "ctg004"], "Viruses")
})

test_that("read_taxonomy2 errors on missing file and missing columns", {
  expect_error(read_taxonomy2("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_taxonomy2(bad), "missing columns")
})


# --- read_iphop ----------------------------------------------------------
# Header from iPHoP v1.3.x Host_prediction_to_genus_m90.csv, which is comma
# separated and spells the second column "AAI to closest RaFAH reference".
test_that("read_iphop parses the comma-separated genus prediction file", {
  df <- suppressWarnings(read_iphop(fx("iphop_Host_prediction_to_genus_m90.csv")))
  expect_true(all(c("Contig", "iphop_genus", "iphop_score", "iphop_aai2ref",
                    "iphop_methods") %in% colnames(df)))
  expect_true(is.numeric(df$iphop_score))
})

test_that("read_iphop keeps the highest-scoring host of a multi-row virus", {
  # iPHoP reports one row per candidate host genus.
  df <- suppressWarnings(read_iphop(fx("iphop_Host_prediction_to_genus_m90.csv")))
  expect_equal(sum(df$Contig == "ctg001"), 1)
  expect_equal(df$iphop_genus[df$Contig == "ctg001"], "Escherichia")
  expect_equal(df$iphop_n_predictions[df$Contig == "ctg001"], 2L)
})

test_that("read_iphop errors on missing file", {
  expect_error(read_iphop("/nonexistent/file.csv"), "not found")
})


# --- read_replicyc -------------------------------------------------------
# BACPHLIP leaves the first header field empty; fread names it V1. Copied from
# bacphlip/putative_vcontigs_pref1.fasta.bacphlip of the reference run.
test_that("read_replicyc parses the BACPHLIP table", {
  df <- read_replicyc(fx("bacphlip.bacphlip"), version = "v2")
  expect_true(all(c("Contig", "bacphlip_replicyc") %in% colnames(df)))
  expect_equal(nrow(df), 4)
  expect_equal(df$bacphlip_replicyc[df$Contig == "ctg001"], "virulent")
  expect_equal(df$bacphlip_replicyc[df$Contig == "ctg002"], "temperate")
})

test_that("read_replicyc auto-detects v2 identifiers", {
  df <- read_replicyc(fx("bacphlip.bacphlip"), version = "auto")
  expect_setequal(df$Contig, c("ctg001", "ctg002", "ctg004", "ctg005"))
})

test_that("read_replicyc rejects an unsupported tool instead of returning a table without Contig", {
  expect_error(read_replicyc(fx("bacphlip.bacphlip"), tool = "phatyp"),
               "bacphlip")
})

test_that("read_replicyc errors on missing file", {
  expect_error(read_replicyc("/nonexistent/file.tsv"), "not found")
})


# --- read_genomad --------------------------------------------------------
# Header from genomad/virus_genomad_summary.tsv of the reference run.
test_that("read_genomad parses the virus summary", {
  df <- read_genomad(fx("genomad_summary.tsv"))
  expect_true(all(c("Contig", "genomad_score", "genomad_taxonomy", "genomad_fdr",
                    "genomad_n_hallmarks", "genomad_seq_name") %in% colnames(df)))
  expect_equal(nrow(df), 5)
  expect_equal(df$genomad_score[df$Contig == "ctg004"], 0.99)
})

test_that("read_genomad keeps fdr numeric when geNomad wrote NA in every row", {
  # Without score calibration every fdr is the literal NA, which fread reads as
  # a logical column, so the rowData type would change between runs.
  df <- read_genomad(fx("genomad_summary.tsv"))
  expect_true(is.numeric(df$genomad_fdr))
  expect_true(all(is.na(df$genomad_fdr)))
})

test_that("read_genomad strips the |provirus_ decoration so the row joins to its contig", {
  df <- read_genomad(fx("genomad_summary.tsv"))
  expect_true("ctg003" %in% df$Contig)
  expect_equal(df$genomad_seq_name[df$Contig == "ctg003"], "ctg003|provirus_1200_7400")
})

test_that("read_genomad errors on missing file and missing columns", {
  expect_error(read_genomad("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_genomad(bad), "missing columns")
})


# --- read_dramv ----------------------------------------------------------
# DRAM-v writes an unnamed first column of gene IDs; `fasta` is the input FASTA
# stem, identical on every row.
test_that("read_dramv takes the gene ID from the index column, not from `fasta`", {
  df <- read_dramv(fx("dramv_annotations.tsv"))
  expect_true(all(c("Contig", "gene_id", "dramv_category") %in% colnames(df)))
  expect_equal(nrow(df), 8)
  # Using `fasta` would give every gene the same identifier.
  expect_equal(length(unique(df$gene_id)), 8)
  expect_false(any(df$gene_id == "final-viral-combined-for-dramv"))
})

test_that("read_dramv strips the -cat_N suffix that --prep-for-dramv appends", {
  df <- read_dramv(fx("dramv_annotations.tsv"))
  expect_setequal(unique(df$Contig), c("ctg001", "ctg002", "ctg004"))
})

test_that("read_dramv strips the __full form of the VirSorter2 suffix", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("\tfasta\tscaffold\trank",
               "g1\tstem\tctg001__full-cat_1\tA",
               "g2\tstem\tctg002__0_partial-cat_4\tV"), f)
  df <- read_dramv(f)
  expect_setequal(df$Contig, c("ctg001", "ctg002"))
})

test_that("read_dramv survives a run without the optional database columns", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("\tfasta\tscaffold\trank", "g1\tstem\tctg001\tA"), f)
  df <- read_dramv(f)
  expect_equal(nrow(df), 1)
  expect_true(is.na(df$dramv_ko))
})

test_that("read_dramv errors on missing file", {
  expect_error(read_dramv("/nonexistent/file.tsv"), "not found")
})

test_that(".summarize_dramv aggregates per contig", {
  df <- read_dramv(fx("dramv_annotations.tsv"))
  s <- vpfkit:::.summarize_dramv(df)
  expect_true("dramv_amg_count" %in% colnames(s))
  row <- s[s$Contig == "ctg001", ]
  expect_equal(row$dramv_amg_count, 2)
  expect_equal(row$dramv_viral_gene_count, 1)
  expect_equal(row$dramv_total_gene_count, 3)
})


# --- read_checkamg -------------------------------------------------------
# Header from checkamg/checkamg_results/final_results.tsv of the reference run.
# The column names contain spaces.
test_that("read_checkamg parses the CheckAMG final results", {
  df <- read_checkamg(fx("checkamg_final_results.tsv"))
  expect_true(all(c("Contig", "gene_id", "checkamg_class", "checkamg_confidence",
                    "checkamg_function", "checkamg_kegg_ko") %in% colnames(df)))
  expect_equal(nrow(df), 6)
  expect_setequal(unique(df$checkamg_class),
                  c("metabolic", "unclassified", "regulatory", "physiological"))
  expect_setequal(unique(df$checkamg_confidence), c("low", "high", "medium"))
})

test_that("read_checkamg accepts the results directory the pipeline publishes", {
  d <- withr::local_tempdir()
  file.copy(fx("checkamg_final_results.tsv"), file.path(d, "final_results.tsv"))
  df <- read_checkamg(d)
  expect_equal(nrow(df), 6)
})

test_that("read_checkamg errors on a directory without final_results.tsv", {
  d <- withr::local_tempdir()
  expect_error(read_checkamg(d), "final_results.tsv")
})

test_that(".summarize_checkamg counts AMG, APG and AReG per contig", {
  df <- read_checkamg(fx("checkamg_final_results.tsv"))
  s <- vpfkit:::.summarize_checkamg(df)
  row <- s[s$Contig == "ctg001", ]
  expect_equal(row$checkamg_gene_count, 3L)
  expect_equal(row$checkamg_amg_count, 1L)
  expect_equal(row$checkamg_areg_count, 1L)
  expect_equal(s$checkamg_apg_count[s$Contig == "ctg002"], 1L)
})


# --- read_vrhyme ---------------------------------------------------------
# vRhyme v1.1.0 writes `scaffold` then `bin`, with plain integer bin labels.
test_that("read_vrhyme parses the membership file", {
  df <- read_vrhyme(fx("vrhyme_best_bins_membership.tsv"))
  expect_true(all(c("Contig", "vrhyme_bin") %in% colnames(df)))
  expect_equal(nrow(df), 3)
  expect_equal(df$vrhyme_bin[df$Contig == "ctg001"], df$vrhyme_bin[df$Contig == "ctg004"])
})

test_that("read_vrhyme also accepts a generic contig/bin header", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("bin\tcontig", "1\tctg001"), f)
  df <- read_vrhyme(f)
  expect_equal(df$Contig, "ctg001")
})

test_that("read_vrhyme errors on missing file and on an unrecognizable header", {
  expect_error(read_vrhyme("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_vrhyme(bad), "missing columns")
})


# --- read_phist / read_pharokka / read_catbat ----------------------------
# No process in the current ViroProfiler pipeline produces these; the fixtures
# are synthetic and the readers are kept for externally produced tables.
test_that("read_phist parses predictions", {
  df <- read_phist(fx("phist_predictions.tsv"))
  expect_true(all(c("Contig", "phist_host", "phist_score") %in% colnames(df)))
  expect_equal(df$phist_host[df$Contig == "ctg001"], "Escherichia_coli")
})

test_that("read_phist errors on missing file and missing columns", {
  expect_error(read_phist("/nonexistent/file.tsv"), "not found")
  bad <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("col_a\tcol_b", "1\t2"), bad)
  expect_error(read_phist(bad), "missing columns")
})

test_that("read_pharokka parses the CDS table", {
  df <- read_pharokka(fx("pharokka_cds.tsv"))
  expect_true(all(c("Contig", "gene_id", "pharokka_function", "pharokka_card")
                  %in% colnames(df)))
  expect_equal(nrow(df), 7)
})

test_that(".summarize_pharokka aggregates per contig", {
  df <- read_pharokka(fx("pharokka_cds.tsv"))
  s <- vpfkit:::.summarize_pharokka(df)
  row <- s[s$Contig == "ctg002", ]
  expect_equal(row$pharokka_card_count, 1)
  expect_equal(row$pharokka_vfdb_count, 0)
})

test_that("read_catbat returns a Contig column like every other reader", {
  df <- read_catbat(fx("catbat_taxonomy.tsv"))
  expect_true("Contig" %in% colnames(df))
  expect_true(any(grepl("^CATBAT_", colnames(df))))
  expect_equal(nrow(df), 2)
})

test_that("read_catbat errors on missing file", {
  expect_error(read_catbat("/nonexistent/file.tsv"), "not found")
})


# --- read_coverm_log -----------------------------------------------------
test_that("read_coverm_log extracts library size and mapping rate", {
  f <- withr::local_tempfile(fileext = ".txt")
  writeLines(c(
    "[2026-08-13T14:40:20Z INFO  bird_tool_utils::clap_utils] CoverM version 0.8.0",
    "[2026-08-13T14:40:20Z INFO  coverm::contig] In sample 'HT02', found 4703 reads mapped out of 4706 total (99.94%)",
    "[2026-08-13T14:40:20Z INFO  coverm::contig] In sample 'UC20', found 12668 reads mapped out of 12668 total (100.00%)"), f)
  df <- read_coverm_log(f)
  expect_equal(nrow(df), 2)
  expect_equal(df$n_reads_total[df$sample_id == "HT02"], 4706)
  expect_equal(df$mapping_rate[df$sample_id == "UC20"], 1)
})

test_that("read_coverm_log warns when the log has no mapping lines", {
  f <- withr::local_tempfile(fileext = ".txt")
  writeLines("nothing useful here", f)
  expect_warning(read_coverm_log(f), "no 'In sample")
})


# --- .detect_vp_version --------------------------------------------------
test_that(".detect_vp_version detects the v1 -cat_ suffix", {
  detect <- vpfkit:::.detect_vp_version
  expect_equal(detect(c("ctg001-cat_1", "ctg002-cat_3")), "v1")
  expect_equal(detect(c("ctg001", "ctg002")), "v2")
  expect_equal(detect(c("ctg001", "ctg002-cat_5")), "v1")
  expect_equal(detect(character(0)), "v2")
})


# --- .is_absent_path -----------------------------------------------------
test_that("optional file arguments treat NA, empty and 'null' as not supplied", {
  # argparser fills unsupplied arguments with NA, Nextflow can interpolate an
  # empty string, and a shell can pass the literal "null".
  absent <- vpfkit:::.is_absent_path
  expect_true(absent(NULL))
  expect_true(absent(NA))
  expect_true(absent(NA_character_))
  expect_true(absent(""))
  expect_true(absent("  "))
  expect_true(absent("null"))
  expect_true(absent("NULL"))
  expect_true(absent("none"))
  expect_false(absent("/some/path.tsv"))
})

test_that("create_vpftse drops an optional path that does not exist, with a warning", {
  expect_warning(
    tse <- build_fixture_tse(fin_genomad = "/nonexistent/genomad.tsv"),
    "does not exist")
  expect_false("genomad_score" %in% colnames(SummarizedExperiment::rowData(tse)))
})

test_that("create_vpftse accepts NA for every optional path", {
  tse <- build_fixture_tse(fin_genomad = NA, fin_vrhyme = NA, fin_phist = NA,
                           fin_dramv = NA, fin_pharokka = NA, fin_checkamg = NA,
                           fin_iphop = NA, fin_metadata = NA, fin_coverm_log = NA,
                           fin_vircontigs = NA)
  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(nrow(tse), 5)
})


# --- create_vpftse -------------------------------------------------------
test_that("create_vpftse assembles a valid TSE", {
  tse <- build_fixture_tse()

  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(nrow(tse), 5)
  expect_equal(ncol(tse), 3)
  rd <- SummarizedExperiment::rowData(tse)
  expect_true("checkv_quality" %in% colnames(rd))
  expect_true("virsorter2_max_score" %in% colnames(rd))
})

test_that("create_vpftse names the coverage-depth assay trimmed_mean, not tmm", {
  # CoverM `--methods trimmed_mean` is the trimmed mean of per-base coverage
  # depth. TMM elsewhere means edgeR's trimmed mean of M-values, which this
  # pipeline never computes.
  tse <- build_fixture_tse()
  expect_setequal(SummarizedExperiment::assayNames(tse),
                  c("counts", "tpm", "trimmed_mean", "covfrac"))
  expect_false("tmm" %in% SummarizedExperiment::assayNames(tse))
})

test_that("create_vpftse records what each assay measures", {
  tse <- build_fixture_tse()
  info <- S4Vectors::metadata(tse)$viroprofiler
  expect_true(all(c("counts", "tpm", "trimmed_mean", "covfrac") %in% names(info$assays)))
  expect_match(info$assays$trimmed_mean, "not edgeR", ignore.case = TRUE)
  expect_true(!is.null(info$source_files$taxonomy))
})

test_that("create_vpftse refuses abundance files with different contig sets", {
  short <- withr::local_tempfile(fileext = ".tsv")
  lines <- readLines(fx("coverm_abundance_tpm.tsv"))
  writeLines(lines[1:4], short)
  expect_error(build_fixture_tse_alt <- create_vpftse(
    fin_abcount    = fx("coverm_abundance_count.tsv"),
    fin_abtpm      = short,
    fin_abtmm      = fx("coverm_abundance_tmm.tsv"),
    fin_abcov      = fx("coverm_abundance_covfrac.tsv"),
    fin_taxa       = fx("taxonomy_tse.tsv"),
    fin_checkv     = fx("checkv_quality_summary.tsv"),
    fin_virsorter2 = fx("virsorter2_final_viral_score.tsv"),
    fin_vibrant    = fx("vibrant_quality.tsv"),
    fin_dvf        = fx("dvf_virus.tsv"),
    fin_replicyc   = fx("bacphlip.bacphlip")), "different contig set")
})

test_that("create_vpftse reorders an abundance file whose contigs are shuffled", {
  shuffled <- withr::local_tempfile(fileext = ".tsv")
  lines <- readLines(fx("coverm_abundance_tpm.tsv"))
  writeLines(c(lines[1], rev(lines[-1])), shuffled)
  tse <- create_vpftse(
    fin_abcount    = fx("coverm_abundance_count.tsv"),
    fin_abtpm      = shuffled,
    fin_abtmm      = fx("coverm_abundance_tmm.tsv"),
    fin_abcov      = fx("coverm_abundance_covfrac.tsv"),
    fin_taxa       = fx("taxonomy_tse.tsv"),
    fin_checkv     = fx("checkv_quality_summary.tsv"),
    fin_virsorter2 = fx("virsorter2_final_viral_score.tsv"),
    fin_vibrant    = fx("vibrant_quality.tsv"),
    fin_dvf        = fx("dvf_virus.tsv"),
    fin_replicyc   = fx("bacphlip.bacphlip"))
  # ctg001 must carry its own TPM, not the value of the contig that happened to
  # sit in the same row of the shuffled file.
  expect_equal(SummarizedExperiment::assay(tse, "tpm")["ctg001", "HT02"], 15.5)
})

test_that("create_vpftse warns when a tool table matches no contig at all", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("seq_name\tvirus_score", "totally_different_id\t0.99"), f)
  expect_warning(build_fixture_tse(fin_genomad = f), "match the")
})

test_that("create_vpftse joins optional annotations", {
  tse <- build_fixture_tse(fin_genomad = fx("genomad_summary.tsv"),
                           fin_vrhyme = fx("vrhyme_best_bins_membership.tsv"),
                           fin_iphop = fx("iphop_Host_prediction_to_genus_m90.csv"))
  rd <- SummarizedExperiment::rowData(tse)
  expect_true(all(c("genomad_score", "vrhyme_bin", "iphop_genus") %in% colnames(rd)))
  expect_equal(rd[["genomad_score"]][rownames(rd) == "ctg004"], 0.99)
})


# --- colData -------------------------------------------------------------
test_that("create_vpftse builds colData from sample names when no metadata is given", {
  tse <- build_fixture_tse()
  cd <- SummarizedExperiment::colData(tse)
  expect_equal(rownames(cd), c("HT02", "UC20", "HT03"))
  expect_equal(as.character(cd$sample_name), c("HT02", "UC20", "HT03"))
})

test_that("create_vpftse joins a sample metadata file by name, not by position", {
  meta <- withr::local_tempfile(fileext = ".csv")
  # Deliberately in a different order than the assay columns.
  writeLines(c("sample_id,group,age",
               "HT03,control,45",
               "HT02,case,31",
               "UC20,control,52"), meta)
  tse <- build_fixture_tse(fin_metadata = meta)
  cd <- SummarizedExperiment::colData(tse)
  expect_equal(rownames(cd), c("HT02", "UC20", "HT03"))
  expect_equal(as.character(cd$group), c("case", "control", "control"))
  expect_equal(as.numeric(cd$age), c(31, 52, 45))
})

test_that("create_vpftse warns about metadata rows that match no sample", {
  meta <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("sample_id,group", "HT02,case", "NOT_A_SAMPLE,control"), meta)
  expect_warning(expect_warning(build_fixture_tse(fin_metadata = meta),
                                "no row for"),
                 "match no sample")
})

test_that("create_vpftse refuses both df_metadata and fin_metadata", {
  meta <- withr::local_tempfile(fileext = ".csv")
  writeLines(c("sample_id,group", "HT02,case"), meta)
  expect_error(build_fixture_tse(fin_metadata = meta,
                                 df_metadata = data.frame(sample_id = "HT02")),
               "not both")
})

test_that("create_vpftse adds sequencing depth from the CoverM log", {
  logf <- withr::local_tempfile(fileext = ".txt")
  writeLines(c(
    "[INFO coverm::contig] In sample 'HT02', found 4703 reads mapped out of 4706 total (99.94%)",
    "[INFO coverm::contig] In sample 'UC20', found 12668 reads mapped out of 12668 total (100.00%)",
    "[INFO coverm::contig] In sample 'HT03', found 900 reads mapped out of 1000 total (90.00%)"), logf)
  tse <- build_fixture_tse(fin_coverm_log = logf)
  cd <- SummarizedExperiment::colData(tse)
  expect_equal(as.numeric(cd$n_reads_total), c(4706, 12668, 1000))
  expect_equal(as.numeric(cd$mapping_rate)[3], 0.9)
})


# --- metadata ------------------------------------------------------------
test_that("create_vpftse stores gene-level annotations in metadata", {
  tse <- build_fixture_tse(fin_dramv = fx("dramv_annotations.tsv"),
                           fin_checkamg = fx("checkamg_final_results.tsv"))
  ga <- S4Vectors::metadata(tse)$gene_annotations
  expect_s3_class(ga, "data.frame")
  expect_true("source" %in% colnames(ga))
  expect_setequal(unique(ga$source), c("dramv", "checkamg"))
  expect_setequal(names(S4Vectors::metadata(tse)$gene_annotations_by_tool),
                  c("dramv", "checkamg"))
  rd <- SummarizedExperiment::rowData(tse)
  expect_true(all(c("dramv_amg_count", "checkamg_amg_count") %in% colnames(rd)))
})

test_that("create_vpftse keeps rowData aligned after a gene-level summary is merged", {
  tse <- build_fixture_tse(fin_dramv = fx("dramv_annotations.tsv"))
  rd <- SummarizedExperiment::rowData(tse)
  expect_equal(rownames(rd), rownames(SummarizedExperiment::assay(tse, "counts")))
  expect_equal(rd[["dramv_amg_count"]][rownames(rd) == "ctg001"], 2)
})

test_that("create_vpftse records pipeline provenance when it is supplied", {
  tse <- build_fixture_tse(pipeline_info = list(version = "1.2.3", profile = "test"))
  expect_equal(S4Vectors::metadata(tse)$viroprofiler$pipeline$version, "1.2.3")
})


# --- normalize_assay_names -----------------------------------------------
test_that("normalize_assay_names renames a legacy tmm assay and warns", {
  tse <- build_fixture_tse()
  SummarizedExperiment::assayNames(tse) <- c("counts", "tpm", "tmm", "covfrac")
  expect_warning(out <- normalize_assay_names(tse), "trimmed_mean")
  expect_true("trimmed_mean" %in% SummarizedExperiment::assayNames(out))
  expect_false("tmm" %in% SummarizedExperiment::assayNames(out))
})

test_that("normalize_assay_names leaves a current object untouched", {
  tse <- build_fixture_tse()
  expect_silent(out <- normalize_assay_names(tse))
  expect_setequal(SummarizedExperiment::assayNames(out),
                  c("counts", "tpm", "trimmed_mean", "covfrac"))
})


# --- annotate_viral_votes / create_vpftse_vir ----------------------------
test_that("create_vpftse_vir filters to viral contigs", {
  tse <- build_fixture_tse(fin_genomad = fx("genomad_summary.tsv"))
  tse_vir <- create_vpftse_vir(tse)
  expect_s4_class(tse_vir, "TreeSummarizedExperiment")
  expect_lte(nrow(tse_vir), nrow(tse))
  expect_gte(nrow(tse_vir), 1)
})

test_that("annotate_viral_votes records the evidence behind every contig", {
  tse <- annotate_viral_votes(build_fixture_tse(fin_genomad = fx("genomad_summary.tsv")))
  rd <- SummarizedExperiment::rowData(tse)
  expect_true(all(c("viral_vote_taxonomy", "viral_vote_checkv", "viral_vote_genomad",
                    "viral_vote_n", "viral_vote_evidence", "viral_selected")
                  %in% colnames(rd)))
  expect_true(is.logical(rd$viral_selected))
  # ctg001 is called by every source; the evidence string lists them.
  expect_true(grepl("genomad", rd[["viral_vote_evidence"]][rownames(rd) == "ctg001"]))
  sel <- S4Vectors::metadata(tse)$viral_selection
  expect_equal(sel$rule, "vote")
  expect_equal(sel$n_total, 5)
  expect_true(!is.null(sel$n_by_vote))
})

test_that("a missing vote column is skipped rather than collapsing the whole result", {
  # rowData$missing returns NULL, !is.na(NULL) is logical(0), and OR-ing a
  # zero-length vector yields logical(0), so tse[logical(0), ] used to return an
  # empty object with no error at all.
  tse <- build_fixture_tse()
  rd <- SummarizedExperiment::rowData(tse)
  SummarizedExperiment::rowData(tse) <- rd[, setdiff(colnames(rd), "Domain"), drop = FALSE]
  expect_warning(tse_vir <- create_vpftse_vir(tse), "votes skipped")
  expect_gt(nrow(tse_vir), 0)
})

test_that("create_vpftse_vir errors when no vote column is available at all", {
  tse <- build_fixture_tse()
  rd <- SummarizedExperiment::rowData(tse)
  vote_cols <- c("Domain", "checkv_quality", "virsorter2_max_score_group",
                 "vibrant_quality", "dvf_score", "genomad_score")
  SummarizedExperiment::rowData(tse) <- rd[, setdiff(colnames(rd), vote_cols), drop = FALSE]
  expect_error(create_vpftse_vir(tse), "None of the requested viral-identity votes")
})

test_that("an empty-string taxonomy rank does not count as a viral vote", {
  tse <- build_fixture_tse()
  rd <- SummarizedExperiment::rowData(tse)
  rd$Domain <- c("Duplodnaviria", "", "  ", NA, "Viruses")
  SummarizedExperiment::rowData(tse) <- rd
  tse <- annotate_viral_votes(tse, votes = "taxonomy")
  expect_equal(unname(SummarizedExperiment::rowData(tse)$viral_vote_taxonomy),
               c(TRUE, FALSE, FALSE, FALSE, TRUE))
})

test_that("the geNomad threshold is a parameter", {
  tse <- build_fixture_tse(fin_genomad = fx("genomad_summary.tsv"))
  strict <- annotate_viral_votes(tse, votes = "genomad", genomad_min_score = 0.99)
  loose <- annotate_viral_votes(tse, votes = "genomad", genomad_min_score = 0.1)
  expect_lt(sum(SummarizedExperiment::rowData(strict)$viral_selected),
            sum(SummarizedExperiment::rowData(loose)$viral_selected))
  expect_equal(S4Vectors::metadata(strict)$viral_selection$thresholds$genomad_min_score,
               0.99)
})

test_that("create_vpftse_vir can use the pipeline's own candidate list", {
  lst <- withr::local_tempfile(fileext = ".list")
  writeLines(c("ctg001", "ctg004"), lst)
  tse <- build_fixture_tse(fin_vircontigs = lst)
  expect_true("upstream_viral_candidate" %in%
                colnames(SummarizedExperiment::rowData(tse)))
  tse_vir <- create_vpftse_vir(tse, rule = "candidate")
  expect_equal(sort(rownames(tse_vir)), c("ctg001", "ctg004"))
  expect_equal(S4Vectors::metadata(tse_vir)$viral_selection$rule, "candidate")
})

test_that("create_vpftse_vir errors when the candidate rule is asked for without the column", {
  tse <- build_fixture_tse()
  expect_error(create_vpftse_vir(tse, rule = "candidate"), "upstream_viral_candidate")
})

test_that("create_vpftse_vir warns instead of silently returning nothing", {
  tse <- build_fixture_tse()
  rd <- SummarizedExperiment::rowData(tse)
  rd$Domain <- NA_character_
  rd$checkv_quality <- factor(NA, levels = levels(rd$checkv_quality))
  rd$virsorter2_max_score_group <- NA_character_
  rd$vibrant_quality <- factor(NA, levels = levels(rd$vibrant_quality))
  rd$dvf_score <- NA_real_
  SummarizedExperiment::rowData(tse) <- rd
  expect_warning(tse_vir <- create_vpftse_vir(tse), "zero rows")
  expect_equal(nrow(tse_vir), 0)
})


# --- abundance_adjust_by_covfrac -----------------------------------------
test_that("abundance_adjust_by_covfrac masks by breadth of coverage", {
  ab <- matrix(c(100, 200, 300, 400), 2, 2,
               dimnames = list(c("ctg1", "ctg2"), c("s1", "s2")))
  cf <- matrix(c(0.9, 0.2, 0.8, 0.76), 2, 2,
               dimnames = list(c("ctg1", "ctg2"), c("s1", "s2")))
  out <- abundance_adjust_by_covfrac(ab, cf)
  expect_equal(unname(out["ctg2", "s1"]), 0)
  expect_equal(unname(out["ctg1", "s1"]), 100)
  expect_equal(unname(out["ctg2", "s2"]), 400)
})

test_that("abundance_adjust_by_covfrac defaults to the 0.75 breadth cutoff", {
  # Roux et al. 2017 (PeerJ 5:e3817). The ViroProfiler paper used 0.5, which was
  # this function's previous default.
  expect_equal(formals(abundance_adjust_by_covfrac)$covfrac_threshold, 0.75)
  ab <- matrix(10, 1, 1, dimnames = list("ctg1", "s1"))
  cf <- matrix(0.6, 1, 1, dimnames = list("ctg1", "s1"))
  expect_equal(unname(abundance_adjust_by_covfrac(ab, cf)[1, 1]), 0)
  expect_equal(unname(abundance_adjust_by_covfrac(ab, cf, covfrac_threshold = 0.5)[1, 1]), 10)
})

test_that("abundance_adjust_by_covfrac refuses mismatched tables", {
  ab <- matrix(1:4, 2, 2, dimnames = list(c("ctg1", "ctg2"), c("s1", "s2")))
  cf <- matrix(1, 2, 2, dimnames = list(c("ctg1", "ctg9"), c("s1", "s2")))
  expect_error(abundance_adjust_by_covfrac(ab, cf), "different contigs or")
})

test_that("abundance_adjust_by_covfrac realigns a permuted coverage table", {
  ab <- matrix(c(10, 20), 2, 1, dimnames = list(c("ctg1", "ctg2"), "s1"))
  cf <- matrix(c(0.1, 0.9), 2, 1, dimnames = list(c("ctg2", "ctg1"), "s1"))
  out <- abundance_adjust_by_covfrac(ab, cf)
  expect_equal(unname(out[, 1]), c(10, 0))
})

test_that("abundance_adjust_by_covfrac treats a missing coverage fraction as absent", {
  ab <- matrix(c(10, 20), 2, 1, dimnames = list(c("ctg1", "ctg2"), "s1"))
  cf <- matrix(c(0.9, NA), 2, 1, dimnames = list(c("ctg1", "ctg2"), "s1"))
  out <- suppressMessages(abundance_adjust_by_covfrac(ab, cf))
  expect_equal(unname(out[, 1]), c(10, 0))
})

test_that("abundance_adjust_by_covfrac can restore column totals after masking", {
  ab <- matrix(c(90, 10), 2, 1, dimnames = list(c("ctg1", "ctg2"), "s1"))
  cf <- matrix(c(0.9, 0.1), 2, 1, dimnames = list(c("ctg1", "ctg2"), "s1"))
  out <- abundance_adjust_by_covfrac(ab, cf, renormalize = TRUE)
  expect_equal(sum(out), sum(ab))
})


# --- batch_create_vpftse -------------------------------------------------
test_that("batch_create_vpftse auto-discovers fixture files", {
  vpdir <- withr::local_tempdir()
  file.copy(fx("coverm_abundance_count.tsv"), file.path(vpdir, "abundance_contigs_count.tsv"))
  file.copy(fx("coverm_abundance_tpm.tsv"), file.path(vpdir, "abundance_contigs_tpm.tsv"))
  file.copy(fx("coverm_abundance_tmm.tsv"), file.path(vpdir, "abundance_contigs_trimmed_mean.tsv"))
  file.copy(fx("coverm_abundance_covfrac.tsv"), file.path(vpdir, "abundance_contigs_covered_fraction.tsv"))
  file.copy(fx("taxonomy_tse.tsv"), file.path(vpdir, "taxonomy_tse.tsv"))
  file.copy(fx("checkv_quality_summary.tsv"), file.path(vpdir, "quality_summary.tsv"))
  file.copy(fx("virsorter2_final_viral_score.tsv"), file.path(vpdir, "final-viral-score.tsv"))
  file.copy(fx("vibrant_quality.tsv"), file.path(vpdir, "VIBRANT_genome_quality_contigs.tsv"))
  file.copy(fx("dvf_virus.tsv"), file.path(vpdir, "dvf_virus.tsv"))
  file.copy(fx("bacphlip.bacphlip"), file.path(vpdir, "contigs.bacphlip"))

  tse <- suppressWarnings(batch_create_vpftse(vpdir))
  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(nrow(tse), 5)
  expect_equal(ncol(tse), 3)
  expect_true("checkv_quality" %in% colnames(SummarizedExperiment::rowData(tse)))
})

test_that("batch_create_vpftse errors on a missing directory or a missing required file", {
  expect_error(batch_create_vpftse("/nonexistent/dir"), "not found")
  expect_error(batch_create_vpftse(withr::local_tempdir()), "Required file not found")
})


####################################################################################
## Tests against the reference ViroProfiler run
##
## These are skipped when the reference data is absent, which is the normal case
## on CI. They are the only tests that prove the readers agree with real output
## rather than with a fixture written from the same assumption as the code.
####################################################################################

test_that("every reader parses the reference run's real output", {
  skip_without_real_data()

  cv <- read_checkv(file.path(REAL_VP, "checkv/quality_summary.tsv"))
  expect_gt(nrow(cv), 0)
  expect_true(is.factor(cv$checkv_quality))
  expect_false(anyNA(cv$checkv_quality))

  vs <- read_virsorter2(file.path(REAL_VP, "virsorter2/out_vs2/final-viral-score.tsv"))
  expect_gt(nrow(vs), 0)

  vb <- read_vibrant(file.path(
    REAL_VP, "vibrant/VIBRANT_contigs/VIBRANT_results_contigs/VIBRANT_genome_quality_contigs.tsv"))
  expect_gt(nrow(vb), 0)
  expect_false(anyNA(vb$vibrant_quality))

  gn <- read_genomad(file.path(REAL_VP, "genomad/virus_genomad_summary.tsv"))
  expect_gt(nrow(gn), 0)
  expect_true(is.numeric(gn$genomad_fdr))
  expect_true(all(gn$genomad_score >= 0 & gn$genomad_score <= 1))

  tx <- read_taxonomy2(file.path(REAL_VP, "taxonomy/taxonomy_tse.tsv"))
  expect_gt(nrow(tx), 0)
  expect_false(any(!is.na(tx$Species) & !nzchar(tx$Species)))

  bp <- read_replicyc(file.path(REAL_VP, "bacphlip/putative_vcontigs_pref1.fasta.bacphlip"))
  expect_gt(nrow(bp), 0)
  expect_true(all(bp$bacphlip_replicyc %in% c("virulent", "temperate")))

  ca <- read_checkamg(file.path(REAL_VP, "checkamg/checkamg_results"))
  expect_gt(nrow(ca), 0)
  expect_true(all(ca$checkamg_class %in%
                    c("metabolic", "physiological", "regulatory", "unclassified")))
  expect_true(all(ca$checkamg_confidence %in% c("high", "medium", "low")))
})

test_that("every reader's contig IDs join to the real abundance matrix", {
  skip_without_real_data()
  ab <- read_coverm(file.path(REAL_VP, "abundance/abundance_contigs_count.tsv.gz"))
  contigs <- rownames(ab)

  rate <- function(ids) mean(unique(ids) %in% contigs)
  # Every tool runs on a subset of the assembly, so not every contig is covered;
  # what must hold is that the tool reports no identifier the abundance matrix
  # does not have. A rate below 1 means the identifiers were mangled.
  #
  # CheckV is the one exception, and not because of identifiers: it runs on
  # contigs_cclib_long.fasta (23 sequences) while the abundance matrix is built
  # on contigs_nrclib.fasta (22), after ANI dereplication. Its surplus rows are
  # contigs that clustering removed, and they must be absent from the
  # dereplicated library rather than being mangled names.
  cv_ids <- unique(read_checkv(file.path(REAL_VP, "checkv/quality_summary.tsv"))$Contig)
  surplus <- setdiff(cv_ids, contigs)
  nrclib <- sub("^>", "", grep("^>", readLines(
    file.path(REAL_VP, "contiglib/contigs_nrclib.fasta")), value = TRUE))
  expect_length(intersect(surplus, nrclib), 0)
  expect_gte(rate(cv_ids), 0.95)
  expect_equal(rate(read_virsorter2(file.path(
    REAL_VP, "virsorter2/out_vs2/final-viral-score.tsv"))$Contig), 1)
  expect_equal(rate(read_vibrant(file.path(
    REAL_VP,
    "vibrant/VIBRANT_contigs/VIBRANT_results_contigs/VIBRANT_genome_quality_contigs.tsv"))$Contig), 1)
  expect_equal(rate(read_genomad(file.path(
    REAL_VP, "genomad/virus_genomad_summary.tsv"))$Contig), 1)
  expect_equal(rate(read_taxonomy2(file.path(REAL_VP, "taxonomy/taxonomy_tse.tsv"))$Contig), 1)
  expect_equal(rate(read_replicyc(file.path(
    REAL_VP, "bacphlip/putative_vcontigs_pref1.fasta.bacphlip"))$Contig), 1)
  expect_equal(rate(read_checkamg(file.path(REAL_VP, "checkamg/checkamg_results"))$Contig), 1)
})

test_that("create_vpftse rebuilds the reference run and reproduces its viral subset", {
  skip_without_real_data()
  # ViroProfiler feeds assets/no_dvf_scores.tsv into the DVF slot because
  # DeepVirFinder was removed from the pipeline; reconstruct that sentinel here
  # rather than depend on the pipeline checkout being present.
  dvf_sentinel <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("name\tlen\tscore\tpvalue\tqvalue",
               "__no_dvf_placeholder__\t0\t0\t1\t1"), dvf_sentinel)
  tse <- create_vpftse(
    fin_abcount    = file.path(REAL_VP, "abundance/abundance_contigs_count.tsv.gz"),
    fin_abtpm      = file.path(REAL_VP, "abundance/abundance_contigs_tpm.tsv.gz"),
    fin_abtmm      = file.path(REAL_VP, "abundance/abundance_contigs_trimmed_mean.tsv.gz"),
    fin_abcov      = file.path(REAL_VP, "abundance/abundance_contigs_covered_fraction.tsv.gz"),
    fin_taxa       = file.path(REAL_VP, "taxonomy/taxonomy_tse.tsv"),
    fin_checkv     = file.path(REAL_VP, "checkv/quality_summary.tsv"),
    fin_virsorter2 = file.path(REAL_VP, "virsorter2/out_vs2/final-viral-score.tsv"),
    fin_vibrant    = file.path(
      REAL_VP,
      "vibrant/VIBRANT_contigs/VIBRANT_results_contigs/VIBRANT_genome_quality_contigs.tsv"),
    fin_dvf        = dvf_sentinel,
    fin_replicyc   = file.path(REAL_VP, "bacphlip/putative_vcontigs_pref1.fasta.bacphlip"),
    fin_genomad    = file.path(REAL_VP, "genomad/virus_genomad_summary.tsv"),
    fin_checkamg   = file.path(REAL_VP, "checkamg/checkamg_results"),
    fin_coverm_log = file.path(REAL_VP, "abundance/log_contig_count.txt"),
    fin_vircontigs = file.path(REAL_VP, "vircontigs/putative_vcontigs_pref1.list"))

  expect_equal(ncol(tse), 2)
  expect_equal(nrow(tse), 22)
  expect_setequal(SummarizedExperiment::assayNames(tse),
                  c("counts", "tpm", "trimmed_mean", "covfrac"))

  cd <- SummarizedExperiment::colData(tse)
  expect_true(all(c("n_reads_total", "n_reads_mapped", "mapping_rate") %in% colnames(cd)))
  expect_equal(as.numeric(cd["HT02", "n_reads_total"]), 4706)

  # The reference run's published viral object has 18 contigs.
  tse_vir <- create_vpftse_vir(tse)
  expect_equal(nrow(tse_vir), 18)

  ref <- readRDS(file.path(REAL_VP, "results/viroprofiler_output.rds"))
  expect_setequal(rownames(tse_vir), rownames(ref))

  # The candidate rule reproduces the same set from the pipeline's own list.
  tse_cand <- create_vpftse_vir(tse, rule = "candidate", reannotate = TRUE)
  expect_equal(nrow(tse_cand), 18)
})

test_that("the reference run's DVF slot is a sentinel that matches no contig", {
  skip_without_real_data()
  sentinel <- file.path(REAL_VP, "..", "..", "..", "..")  # placeholder, see below
  # ViroProfiler feeds assets/no_dvf_scores.tsv, whose single row is
  # __no_dvf_placeholder__, purely because create_vpftse_vir() used to index the
  # dvf_score column unconditionally. Reconstruct it rather than depend on the
  # pipeline checkout being present.
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("name\tlen\tscore\tpvalue\tqvalue",
               "__no_dvf_placeholder__\t0\t0\t1\t1"), f)
  df <- read_dvf(f, thr_qvalue = 0.1)
  expect_equal(nrow(df), 0)
})

test_that("legacy objects saved with the tmm assay name still load", {
  skip_without_real_data()
  ref <- readRDS(file.path(REAL_VP, "results/viroprofiler_output.rds"))
  expect_true("tmm" %in% SummarizedExperiment::assayNames(ref))
  fixed <- suppressWarnings(normalize_assay_names(ref))
  expect_true("trimmed_mean" %in% SummarizedExperiment::assayNames(fixed))
  expect_equal(SummarizedExperiment::assay(fixed, "trimmed_mean"),
               SummarizedExperiment::assay(ref, "tmm"))
})
