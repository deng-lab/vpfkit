vpf_gene_filter <- function(tse, assay = "counts", raw = tse) {
  list(
    tse = shiny::reactive(tse),
    raw = shiny::reactive(raw),
    assay = shiny::reactive(if (is.null(tse)) NULL else assay),
    audit = shiny::reactive(NULL),
    name = shiny::reactive("test")
  )
}

vpf_gene_demo <- function(pattern = "Gut virome") {
  demos <- vpf_demo_datasets()
  hit <- demos[grepl(pattern, names(demos), ignore.case = TRUE)]
  if (length(hit) == 0) NULL else readRDS(hit[[1]])
}

vpf_gene_real <- function() {
  p <- Sys.getenv("VPFKIT_REAL_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("module ui works", {
  ui <- mod_genes_ui(id = "test")
  golem::expect_shinytaglist(ui)
  expect_true("id" %in% names(formals(mod_genes_ui)))
})

test_that("no dataset renders placeholders", {
  testServer(mod_genes_server, args = list(r_filter = vpf_gene_filter(NULL)), {
    expect_null(gene_table())
    expect_false(is.null(output$status))
    expect_false(is.null(output$plt_phrog))
    expect_false(is.null(output$plt_amg))
    expect_false(is.null(output$plt_card))
    expect_false(is.null(output$plt_vfdb))
    expect_false(is.null(output$plt_per_contig))
  })
})

test_that("gene annotations are read from metadata and restricted to the selection", {
  tse <- vpf_gene_demo()
  skip_if(is.null(tse), "demo dataset not available")
  ga <- vpf_gene_annotations(tse)
  skip_if(is.null(ga), "demo dataset has no gene annotations")
  testServer(mod_genes_server, args = list(r_filter = vpf_gene_filter(tse)), {
    expect_equal(nrow(gene_table()), nrow(ga))
    expect_false(is.null(output$status))
    expect_false(is.null(output$plt_phrog))
    expect_false(is.null(output$plt_amg))
    expect_false(is.null(output$plt_card))
    expect_false(is.null(output$plt_vfdb))
    expect_false(is.null(output$plt_per_contig))
    expect_false(is.null(output$evidence_caveat))
  })

  # Filtering the contigs must shrink the gene table accordingly.
  subset_tse <- tse[1:5, ]
  testServer(mod_genes_server,
             args = list(r_filter = vpf_gene_filter(subset_tse, raw = tse)), {
    expect_lt(nrow(gene_table()), nrow(ga))
    expect_true(all(gene_table()$Contig %in% rownames(subset_tse)))
  })
})

test_that("a dataset with no gene table names the option that produces one", {
  tse <- vpf_gene_real()
  skip_if(is.null(tse), "real ViroProfiler output not available")
  expect_null(vpf_gene_annotations(tse))
  testServer(mod_genes_server, args = list(r_filter = vpf_gene_filter(tse)), {
    expect_null(gene_table())
    expect_false(is.null(output$status))
    expect_false(is.null(output$plt_phrog))
    expect_false(is.null(output$plt_amg))
    expect_false(is.null(output$plt_card))
    expect_false(is.null(output$plt_vfdb))
    expect_false(is.null(output$plt_per_contig))
    expect_null(output$evidence_caveat)
  })
})

test_that("a gene table with no CARD or VFDB hits still renders", {
  tse <- vpf_gene_demo()
  skip_if(is.null(tse), "demo dataset not available")
  ga <- vpf_gene_annotations(tse)
  skip_if(is.null(ga), "demo dataset has no gene annotations")
  ga$pharokka_card <- NA_character_
  ga$pharokka_vfdb <- NA_character_
  S4Vectors::metadata(tse)$gene_annotations <- ga
  testServer(mod_genes_server, args = list(r_filter = vpf_gene_filter(tse)), {
    expect_false(is.null(output$plt_card))
    expect_false(is.null(output$plt_vfdb))
  })
})

test_that("gene annotations survive a zero-contig selection", {
  tse <- vpf_gene_demo()
  skip_if(is.null(tse), "demo dataset not available")
  skip_if(is.null(vpf_gene_annotations(tse)), "demo dataset has no gene annotations")
  empty <- tse[integer(0), ]
  testServer(mod_genes_server,
             args = list(r_filter = vpf_gene_filter(empty, raw = tse)), {
    expect_equal(nrow(gene_table()), 0)
    expect_false(is.null(output$status))
    expect_false(is.null(output$plt_phrog))
    expect_false(is.null(output$plt_per_contig))
  })
})

vpf_gene_enriched <- function() {
  p <- Sys.getenv("VPFKIT_ENRICHED_TSE", unset = "")
  if (!nzchar(p) || !file.exists(p)) {
    p <- "/home/allen/data2/testdata/viroprofiler_real_full/results/viroprofiler_output_enriched.rds"
  }
  if (file.exists(p)) readRDS(p) else NULL
}

test_that("the panels adapt to CheckAMG columns when pharokka and DRAM-v are absent", {
  # A run may carry pharokka columns, DRAM-v columns, CheckAMG columns, or a
  # merge. Hard-coding one naming scheme leaves the tab blank for the others.
  tse <- vpf_gene_enriched()
  skip_if(is.null(tse), "enriched ViroProfiler output not available")
  ga <- vpf_gene_annotations(tse)
  expect_false(is.null(ga))
  expect_true("checkamg_class" %in% colnames(ga))
  expect_false(any(grepl("^pharokka_", colnames(ga))))
  expect_false(any(grepl("^dramv_", colnames(ga))))

  testServer(mod_genes_server, args = list(r_filter = vpf_gene_filter(tse)), {
    expect_equal(nrow(gene_table()), nrow(ga))
    expect_equal(first_column("pharokka_category", "checkamg_phrog"), "checkamg_phrog")
    expect_equal(first_column("dramv_ko", "checkamg_kegg_ko"), "checkamg_kegg_ko")
    expect_equal(first_column("checkamg_function", "checkamg_kegg_ko"), "checkamg_function")
    for (o in c("status", "plt_phrog", "plt_amg", "plt_card", "plt_vfdb",
                "plt_checkamg_class", "plt_checkamg_function", "plt_per_contig",
                "evidence_caveat")) {
      expect_false(is.null(output[[o]]), info = o)
    }
    # Proteins outside a strict viral region are flagged.
    expect_false(is.null(output$checkamg_status))
  })
})

test_that("gene annotations fall back to the per-tool tables", {
  tse <- vpf_gene_enriched()
  skip_if(is.null(tse), "enriched ViroProfiler output not available")
  by_tool <- S4Vectors::metadata(tse)$gene_annotations_by_tool
  skip_if(is.null(by_tool), "no per-tool gene tables")
  stripped <- tse
  S4Vectors::metadata(stripped)$gene_annotations <- NULL
  ga <- vpf_gene_annotations(stripped)
  expect_false(is.null(ga))
  expect_true("checkamg_class" %in% colnames(ga))
})
