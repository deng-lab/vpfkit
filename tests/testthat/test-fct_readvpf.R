# --- read_checkv ---
test_that("read_checkv parses fixture correctly", {
  fin <- test_path("fixtures", "checkv_quality_summary.tsv")
  df <- read_checkv(fin)

  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("checkv_quality" %in% colnames(df))
  expect_true("checkv_contig_length" %in% colnames(df))
  expect_true("checkv_completeness" %in% colnames(df))
  expect_equal(nrow(df), 5)
  expect_true(is.factor(df$checkv_quality))
  expect_equal(levels(df$checkv_quality),
               c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))
})

test_that("read_checkv errors on missing file", {
  expect_error(read_checkv("/nonexistent/file.tsv"), "not found")
})

test_that("read_checkv errors on missing columns", {
  bad_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines("col_a\tcol_b\n1\t2", bad_file)
  expect_error(read_checkv(bad_file), "missing columns")
})

# --- read_virsorter2 ---
test_that("read_virsorter2 parses fixture correctly", {
  fin <- test_path("fixtures", "virsorter2_final_viral_score.tsv")
  df <- read_virsorter2(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("virsorter2_max_score" %in% colnames(df))
  expect_true("virsorter2_max_score_group" %in% colnames(df))
  expect_equal(nrow(df), 5)
})

test_that("read_virsorter2 errors on missing file", {
  expect_error(read_virsorter2("/nonexistent/file.tsv"), "not found")
})

test_that("read_virsorter2 errors on missing columns", {
  bad_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines("col_a\tcol_b\n1\t2", bad_file)
  expect_error(read_virsorter2(bad_file), "missing columns")
})

# --- read_vibrant ---
test_that("read_vibrant parses fixture correctly", {
  fin <- test_path("fixtures", "vibrant_quality.tsv")
  df <- read_vibrant(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("vibrant_quality" %in% colnames(df))
  expect_true("vibrant_replicyc" %in% colnames(df))
  expect_true(is.factor(df$vibrant_quality))
  # ctg001 appears twice (with fragment), dedup keeps one
  expect_equal(sum(df$Contig == "ctg001"), 1)
})

test_that("read_vibrant errors on missing file", {
  expect_error(read_vibrant("/nonexistent/file.tsv"), "not found")
})

# --- read_dvf ---
test_that("read_dvf parses fixture correctly", {
  fin <- test_path("fixtures", "dvf_virus.tsv")
  df <- read_dvf(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("dvf_score" %in% colnames(df))
  # Default thresholds: score>0.9, pvalue<0.01, qvalue<0.01
  # Only ctg001 (0.98, 0.001, 0.002) and ctg004 (0.95, 0.0001, 0.0005) pass
  expect_equal(nrow(df), 2)
  expect_true(all(df$Contig %in% c("ctg001", "ctg004")))
})

test_that("read_dvf errors on missing file", {
  expect_error(read_dvf("/nonexistent/file.tsv"), "not found")
})

# --- read_coverm ---
test_that("read_coverm parses count fixture correctly", {
  fin <- test_path("fixtures", "coverm_abundance_count.tsv")
  df <- read_coverm(fin)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 5)
  expect_true(all(c("Sample_1", "Sample_2", "Sample_3") %in% colnames(df)))
  expect_true("ctg001" %in% rownames(df))
})

test_that("read_coverm errors on missing file", {
  expect_error(read_coverm("/nonexistent/file.tsv"), "not found")
})

# --- read_taxonomy2 ---
test_that("read_taxonomy2 parses fixture correctly", {
  fin <- test_path("fixtures", "mmseqs2_taxa.tsv")
  df <- read_taxonomy2(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("Family" %in% colnames(df))
  # ctg005 has taxa_id=0, should be filtered out
  expect_false("ctg005" %in% df$Contig)
  expect_equal(nrow(df), 3)
})

test_that("read_taxonomy2 errors on missing file", {
  expect_error(read_taxonomy2("/nonexistent/file.tsv"), "not found")
})

test_that("read_taxonomy2 errors on missing columns", {
  bad_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines("col_a\tcol_b\n1\t2", bad_file)
  expect_error(read_taxonomy2(bad_file), "missing columns")
})

# --- read_iphop ---
test_that("read_iphop parses fixture correctly", {
  fin <- test_path("fixtures", "iphop_genus.tsv")
  df <- read_iphop(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("iphop_genus" %in% colnames(df))
  expect_true("iphop_score" %in% colnames(df))
  expect_equal(nrow(df), 3)
})

test_that("read_iphop errors on missing file", {
  expect_error(read_iphop("/nonexistent/file.tsv"), "not found")
})

# --- read_replicyc ---
test_that("read_replicyc parses fixture correctly", {
  fin <- test_path("fixtures", "bacphlip.bacphlip")
  df <- read_replicyc(fin, version = "v2")
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("bacphlip_replicyc" %in% colnames(df))
  expect_equal(nrow(df), 4)
  expect_equal(df$bacphlip_replicyc[df$Contig == "ctg001"], "virulent")
  expect_equal(df$bacphlip_replicyc[df$Contig == "ctg002"], "temperate")
})

test_that("read_replicyc errors on missing file", {
  expect_error(read_replicyc("/nonexistent/file.tsv"), "not found")
})

# --- read_catbat ---
test_that("read_catbat parses fixture correctly", {
  fin <- test_path("fixtures", "catbat_taxonomy.tsv")
  df <- read_catbat(fin)
  expect_s3_class(df, "data.frame")
  expect_true("contig_id" %in% colnames(df))
  expect_true(all(grepl("CATBAT_|contig_id", colnames(df))))
  expect_equal(nrow(df), 2)
})

test_that("read_catbat errors on missing file", {
  expect_error(read_catbat("/nonexistent/file.tsv"), "not found")
})

# --- .detect_vp_version (internal) ---
test_that(".detect_vp_version detects v1 suffix", {
  detect <- vpfkit:::.detect_vp_version
  expect_equal(detect(c("ctg001-cat_1", "ctg002-cat_3")), "v1")
  expect_equal(detect(c("ctg001", "ctg002")), "v2")
  expect_equal(detect(c("ctg001", "ctg002-cat_5")), "v1")
})

# --- read_replicyc with auto version ---
test_that("read_replicyc auto-detects v2 (no -cat_ suffix)", {
  fin <- test_path("fixtures", "bacphlip.bacphlip")
  df <- read_replicyc(fin, version = "auto")
  expect_equal(nrow(df), 4)
  expect_true(all(df$Contig %in% c("ctg001", "ctg002", "ctg004", "ctg005")))
})

# --- create_vpftse ---
test_that("create_vpftse assembles a valid TSE", {
  fix <- test_path("fixtures")
  tse <- create_vpftse(
    fin_abcount   = file.path(fix, "coverm_abundance_count.tsv"),
    fin_abtpm     = file.path(fix, "coverm_abundance_tpm.tsv"),
    fin_abtmm     = file.path(fix, "coverm_abundance_tmm.tsv"),
    fin_abcov     = file.path(fix, "coverm_abundance_covfrac.tsv"),
    fin_taxa      = file.path(fix, "mmseqs2_taxa.tsv"),
    fin_checkv    = file.path(fix, "checkv_quality_summary.tsv"),
    fin_virsorter2 = file.path(fix, "virsorter2_final_viral_score.tsv"),
    fin_vibrant   = file.path(fix, "vibrant_quality.tsv"),
    fin_dvf       = file.path(fix, "dvf_virus.tsv"),
    fin_replicyc  = file.path(fix, "bacphlip.bacphlip")
  )

  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(length(SummarizedExperiment::assays(tse)), 4)
  expect_true(all(c("counts", "tpm", "tmm", "covfrac") %in% SummarizedExperiment::assayNames(tse)))
  expect_equal(nrow(tse), 5)
  expect_equal(ncol(tse), 3)
  rd <- SummarizedExperiment::rowData(tse)
  expect_true("checkv_quality" %in% colnames(rd))
  expect_true("virsorter2_max_score" %in% colnames(rd))
})

# --- create_vpftse_vir ---
test_that("create_vpftse_vir filters to viral contigs", {
  fix <- test_path("fixtures")
  tse <- create_vpftse(
    fin_abcount   = file.path(fix, "coverm_abundance_count.tsv"),
    fin_abtpm     = file.path(fix, "coverm_abundance_tpm.tsv"),
    fin_abtmm     = file.path(fix, "coverm_abundance_tmm.tsv"),
    fin_abcov     = file.path(fix, "coverm_abundance_covfrac.tsv"),
    fin_taxa      = file.path(fix, "mmseqs2_taxa.tsv"),
    fin_checkv    = file.path(fix, "checkv_quality_summary.tsv"),
    fin_virsorter2 = file.path(fix, "virsorter2_final_viral_score.tsv"),
    fin_vibrant   = file.path(fix, "vibrant_quality.tsv"),
    fin_dvf       = file.path(fix, "dvf_virus.tsv"),
    fin_replicyc  = file.path(fix, "bacphlip.bacphlip")
  )

  tse_vir <- create_vpftse_vir(tse)

  expect_s4_class(tse_vir, "TreeSummarizedExperiment")
  expect_lte(nrow(tse_vir), nrow(tse))
  expect_gte(nrow(tse_vir), 1)
})

# --- read_dramv ---
test_that("read_dramv parses fixture correctly", {
  fin <- test_path("fixtures", "dramv_annotations.tsv")
  df <- read_dramv(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("gene_id" %in% colnames(df))
  expect_true("dramv_category" %in% colnames(df))
  expect_equal(nrow(df), 8)
  # ctg001 has 2 AMG genes (rank A), 1 viral
  amg_ctg001 <- sum(df$dramv_category[df$Contig == "ctg001"] == "AMG")
  expect_equal(amg_ctg001, 2)
})

test_that("read_dramv errors on missing file", {
  expect_error(read_dramv("/nonexistent/file.tsv"), "not found")
})

test_that(".summarize_dramv aggregates per contig", {
  fin <- test_path("fixtures", "dramv_annotations.tsv")
  df <- read_dramv(fin)
  summary <- vpfkit:::.summarize_dramv(df)
  expect_true("dramv_amg_count" %in% colnames(summary))
  # ctg001: 2 AMG, 1 viral, 3 total
  row_ctg001 <- summary[summary$Contig == "ctg001", ]
  expect_equal(row_ctg001$dramv_amg_count, 2)
  expect_equal(row_ctg001$dramv_viral_gene_count, 1)
  expect_equal(row_ctg001$dramv_total_gene_count, 3)
})

# --- read_pharokka ---
test_that("read_pharokka parses fixture correctly", {
  fin <- test_path("fixtures", "pharokka_cds.tsv")
  df <- read_pharokka(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("gene_id" %in% colnames(df))
  expect_true("pharokka_function" %in% colnames(df))
  expect_true("pharokka_card" %in% colnames(df))
  expect_equal(nrow(df), 7)
})

test_that("read_pharokka errors on missing file", {
  expect_error(read_pharokka("/nonexistent/file.tsv"), "not found")
})

test_that(".summarize_pharokka aggregates per contig", {
  fin <- test_path("fixtures", "pharokka_cds.tsv")
  df <- read_pharokka(fin)
  summary <- vpfkit:::.summarize_pharokka(df)
  expect_true("pharokka_card_count" %in% colnames(summary))
  # ctg002: 1 CARD hit (beta-lactamase), 0 VFDB
  row_ctg002 <- summary[summary$Contig == "ctg002", ]
  expect_equal(row_ctg002$pharokka_card_count, 1)
  expect_equal(row_ctg002$pharokka_vfdb_count, 0)
})

# --- read_genomad ---
test_that("read_genomad parses fixture correctly", {
  fin <- test_path("fixtures", "genomad_summary.tsv")
  df <- read_genomad(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("genomad_score" %in% colnames(df))
  expect_true("genomad_taxonomy" %in% colnames(df))
  expect_equal(nrow(df), 5)
  expect_equal(df$genomad_score[df$Contig == "ctg004"], 0.99)
})

test_that("read_genomad errors on missing file", {
  expect_error(read_genomad("/nonexistent/file.tsv"), "not found")
})

test_that("read_genomad errors on missing columns", {
  bad_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines("col_a\tcol_b\n1\t2", bad_file)
  expect_error(read_genomad(bad_file), "missing columns")
})

# --- read_vrhyme ---
test_that("read_vrhyme parses fixture correctly", {
  fin <- test_path("fixtures", "vrhyme_bins.tsv")
  df <- read_vrhyme(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("vrhyme_bin" %in% colnames(df))
  expect_equal(nrow(df), 3)
  # ctg001 and ctg004 are in the same bin
  expect_equal(df$vrhyme_bin[df$Contig == "ctg001"], df$vrhyme_bin[df$Contig == "ctg004"])
})

test_that("read_vrhyme errors on missing file", {
  expect_error(read_vrhyme("/nonexistent/file.tsv"), "not found")
})

# --- read_phist ---
test_that("read_phist parses fixture correctly", {
  fin <- test_path("fixtures", "phist_predictions.tsv")
  df <- read_phist(fin)
  expect_s3_class(df, "data.frame")
  expect_true("Contig" %in% colnames(df))
  expect_true("phist_host" %in% colnames(df))
  expect_true("phist_score" %in% colnames(df))
  expect_equal(nrow(df), 3)
  expect_equal(df$phist_host[df$Contig == "ctg001"], "Escherichia_coli")
})

test_that("read_phist errors on missing file", {
  expect_error(read_phist("/nonexistent/file.tsv"), "not found")
})

test_that("read_phist errors on missing columns", {
  bad_file <- withr::local_tempfile(fileext = ".tsv")
  writeLines("col_a\tcol_b\n1\t2", bad_file)
  expect_error(read_phist(bad_file), "missing columns")
})

# --- batch_create_vpftse ---
test_that("batch_create_vpftse auto-discovers fixture files", {
  # Create a temporary directory mimicking ViroProfiler output
  vpdir <- withr::local_tempdir()
  fix <- test_path("fixtures")

  # Copy fixtures with ViroProfiler-standard names
  file.copy(file.path(fix, "coverm_abundance_count.tsv"), file.path(vpdir, "abundance_contigs_count.tsv"))
  file.copy(file.path(fix, "coverm_abundance_tpm.tsv"), file.path(vpdir, "abundance_contigs_tpm.tsv"))
  file.copy(file.path(fix, "coverm_abundance_tmm.tsv"), file.path(vpdir, "abundance_contigs_tmm.tsv"))
  file.copy(file.path(fix, "coverm_abundance_covfrac.tsv"), file.path(vpdir, "abundance_contigs_covfrac.tsv"))
  file.copy(file.path(fix, "mmseqs2_taxa.tsv"), file.path(vpdir, "taxa_mmseqs_formatted_all.tsv"))
  file.copy(file.path(fix, "checkv_quality_summary.tsv"), file.path(vpdir, "quality_summary.tsv"))
  file.copy(file.path(fix, "virsorter2_final_viral_score.tsv"), file.path(vpdir, "final-viral-score.tsv"))
  file.copy(file.path(fix, "vibrant_quality.tsv"), file.path(vpdir, "VIBRANT_genome_quality_contigs.tsv"))
  file.copy(file.path(fix, "dvf_virus.tsv"), file.path(vpdir, "dvf_virus.tsv"))
  file.copy(file.path(fix, "bacphlip.bacphlip"), file.path(vpdir, "contigs.bacphlip"))

  tse <- batch_create_vpftse(vpdir)

  expect_s4_class(tse, "TreeSummarizedExperiment")
  expect_equal(nrow(tse), 5)
  expect_equal(ncol(tse), 3)
  expect_true("checkv_quality" %in% colnames(SummarizedExperiment::rowData(tse)))
})

test_that("batch_create_vpftse errors on nonexistent directory", {
  expect_error(batch_create_vpftse("/nonexistent/dir"), "not found")
})

test_that("batch_create_vpftse errors when required file missing", {
  vpdir <- withr::local_tempdir()
  # Empty directory — should fail looking for required files
  expect_error(batch_create_vpftse(vpdir), "Required file not found")
})
