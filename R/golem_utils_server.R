#' Inverted versions of in, is.null and is.na
#'
#' @noRd
#'
#' @examples
#' 1 %not_in% 1:10
#' not_null(NULL)
`%not_in%` <- Negate(`%in%`)

not_null <- Negate(is.null)

not_na <- Negate(is.na)

#' Removes the null from a vector
#'
#' @noRd
#'
#' @example
#' drop_nulls(list(1, NULL, 2))
drop_nulls <- function(x) {
  x[!sapply(x, is.null)]
}

#' If x is `NULL`, return y, otherwise return x
#'
#' @param x,y Two elements to test, one potentially `NULL`
#'
#' @noRd
#'
#' @examples
#' NULL %||% 1
"%||%" <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

#' If x is `NA`, return y, otherwise return x
#'
#' @param x,y Two elements to test, one potentially `NA`
#'
#' @noRd
#'
#' @examples
#' NA %|NA|% 1
"%|NA|%" <- function(x, y) {
  if (is.na(x)) {
    y
  } else {
    x
  }
}

#' Typing reactiveValues is too long
#'
#' @inheritParams reactiveValues
#' @inheritParams reactiveValuesToList
#'
#' @noRd
rv <- function(...) shiny::reactiveValues(...)
rvtl <- function(...) shiny::reactiveValuesToList(...)


# ---------------------------------------------------------------------------
# ViroProfiler-viewer domain helpers
#
# Everything below is used by the `mod_*` modules of the Shiny application.
# Functions are deliberately free of `shiny` reactive context so they can be
# unit tested directly.
# ---------------------------------------------------------------------------

#' Taxonomic ranks understood by the viewer, from broad to narrow
#' @noRd
VPF_TAXONOMY_RANKS <- c(
  "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
)


# -- Assay semantics --------------------------------------------------------

#' Dictionary of assay names produced by ViroProfiler
#'
#' Every entry describes what the assay actually measures. This matters because
#' the stored names are not self-explanatory: the assay historically written as
#' `tmm` holds CoverM's *trimmed mean coverage depth*, which is unrelated to
#' edgeR's TMM (trimmed mean of M-values) normalization. Using it as if it were
#' a normalized abundance silently produces wrong composition, diversity and
#' ordination results, so the viewer always shows the long description next to
#' the name.
#'
#' `role` is the contract the rest of the app relies on:
#' * `abundance` - may be used for composition, diversity and ordination.
#' * `detection` - a breadth/presence measure; must never be treated as an
#'   abundance.
#'
#' @return A named list of assay descriptors.
#' @noRd
vpf_assay_dictionary <- function() {
  list(
    counts = list(
      label = "Read counts",
      unit = "reads",
      role = "abundance",
      integer_counts = TRUE,
      description = paste(
        "CoverM read count: number of reads mapped to each contig.",
        "Raw integer counts; the only assay valid for count-based models."
      )
    ),
    tpm = list(
      label = "TPM",
      unit = "TPM",
      role = "abundance",
      integer_counts = FALSE,
      description = paste(
        "CoverM TPM: reads per kilobase scaled to a fixed library sum,",
        "so contig length and sequencing depth are both accounted for."
      )
    ),
    trimmed_mean = list(
      label = "Trimmed mean coverage depth",
      unit = "x-fold coverage",
      role = "abundance",
      integer_counts = FALSE,
      description = paste(
        "CoverM trimmed_mean: mean per-base coverage depth after discarding",
        "the 5% highest and 5% lowest covered positions. This is a coverage",
        "DEPTH in x-fold, NOT edgeR's TMM (trimmed mean of M-values)",
        "normalization."
      )
    ),
    tmm = list(
      label = "Trimmed mean coverage depth (stored as 'tmm')",
      unit = "x-fold coverage",
      role = "abundance",
      integer_counts = FALSE,
      description = paste(
        "CoverM trimmed_mean stored under the misleading name 'tmm'.",
        "It is a mean per-base coverage DEPTH in x-fold, NOT edgeR's TMM",
        "(trimmed mean of M-values) normalization. Newer ViroProfiler runs",
        "name this assay 'trimmed_mean'."
      )
    ),
    covfrac = list(
      label = "Covered fraction",
      unit = "fraction 0-1",
      role = "detection",
      integer_counts = FALSE,
      description = paste(
        "Fraction of the contig covered by at least one read (0-1).",
        "A breadth-of-coverage measure used to decide whether a contig is",
        "present at all. It is not an abundance and is excluded from",
        "composition, diversity and ordination."
      )
    ),
    covered_fraction = list(
      label = "Covered fraction",
      unit = "fraction 0-1",
      role = "detection",
      integer_counts = FALSE,
      description = paste(
        "Fraction of the contig covered by at least one read (0-1).",
        "A breadth-of-coverage measure used to decide whether a contig is",
        "present at all. It is not an abundance and is excluded from",
        "composition, diversity and ordination."
      )
    ),
    relabundance = list(
      label = "Relative abundance",
      unit = "proportion",
      role = "abundance",
      integer_counts = FALSE,
      description = "Sample-wise proportions derived from another assay."
    )
  )
}

#' Describe the assays present in a TSE
#'
#' When the object was written by a recent `create_vpftse()` it carries its own
#' assay documentation in `metadata(tse)$viroprofiler$assays`. That description
#' comes from the code that wrote the matrix, so it takes precedence over this
#' package's static dictionary.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A data frame with one row per assay: `assay`, `label`, `unit`,
#'   `role`, `integer_counts`, `description`, `known`.
#' @noRd
vpf_assay_info <- function(tse) {
  empty <- data.frame(
    assay = character(0), label = character(0), unit = character(0),
    role = character(0), integer_counts = logical(0),
    description = character(0), known = logical(0),
    stringsAsFactors = FALSE
  )
  if (is.null(tse)) {
    return(empty)
  }
  nms <- SummarizedExperiment::assayNames(tse)
  if (is.null(nms) || length(nms) == 0) {
    return(empty)
  }
  dict <- vpf_assay_dictionary()
  written <- S4Vectors::metadata(tse)[["viroprofiler"]][["assays"]]
  rows <- lapply(nms, function(nm) {
    d <- dict[[nm]]
    known <- !is.null(d)
    if (!known) {
      d <- list(
        label = nm, unit = "unknown", role = "abundance",
        integer_counts = FALSE,
        description = paste0(
          "Assay '", nm, "' is not part of the ViroProfiler dictionary. ",
          "Its units are unknown; interpret derived statistics with care."
        )
      )
    }
    if (!is.null(written) && !is.null(written[[nm]]) && is.character(written[[nm]])) {
      d$description <- written[[nm]]
      known <- TRUE
    }
    data.frame(
      assay = nm, label = d$label, unit = d$unit, role = d$role,
      integer_counts = isTRUE(d$integer_counts),
      description = d$description, known = known,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Named choices for an abundance-assay `selectInput`
#'
#' Detection assays such as `covfrac` are deliberately excluded: they are not
#' abundances and using them for composition or diversity is a silent error.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A named character vector suitable for `choices =`.
#' @noRd
vpf_abundance_assay_choices <- function(tse) {
  info <- vpf_assay_info(tse)
  info <- info[info$role == "abundance", , drop = FALSE]
  if (nrow(info) == 0) {
    return(character(0))
  }
  stats::setNames(info$assay, paste0(info$label, " (", info$assay, ")"))
}

#' Preferred default abundance assay
#' @noRd
vpf_default_assay <- function(tse) {
  choices <- vpf_abundance_assay_choices(tse)
  if (length(choices) == 0) {
    return(NULL)
  }
  for (pref in c("counts", "tpm", "trimmed_mean", "tmm")) {
    if (pref %in% choices) {
      return(pref)
    }
  }
  unname(choices[[1]])
}

#' Name of the covered-fraction assay, if the object has one
#' @noRd
vpf_covfrac_assay <- function(tse) {
  if (is.null(tse)) {
    return(NULL)
  }
  nms <- SummarizedExperiment::assayNames(tse)
  hit <- intersect(c("covfrac", "covered_fraction"), nms)
  if (length(hit) == 0) NULL else hit[[1]]
}

#' Is this assay raw integer counts?
#'
#' Chao1, rarefaction and negative-binomial models are only defined on raw
#' counts, so the UI has to be able to ask.
#'
#' @noRd
vpf_is_count_assay <- function(tse, assay) {
  if (is.null(tse) || is.null(assay) || !assay %in% SummarizedExperiment::assayNames(tse)) {
    return(FALSE)
  }
  dict <- vpf_assay_dictionary()
  if (!is.null(dict[[assay]]) && !isTRUE(dict[[assay]]$integer_counts)) {
    return(FALSE)
  }
  mat <- SummarizedExperiment::assay(tse, assay)
  if (length(mat) == 0) {
    return(FALSE)
  }
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  length(vals) > 0 && all(vals >= 0) && all(abs(vals - round(vals)) < 1e-8)
}


# -- Taxonomy ---------------------------------------------------------------

#' Turn blank and placeholder strings into `NA`
#'
#' Real ViroProfiler taxonomy tables encode "no assignment at this rank" as an
#' empty string rather than `NA`, which would otherwise show up as an unnamed
#' taxon in every composition plot.
#'
#' @param x A vector.
#' @return A character vector with blanks converted to `NA`.
#' @noRd
vpf_blank_na <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "na", "none", "None", "unknown", "Unknown", "-")] <- NA_character_
  x
}

#' Taxonomic ranks that actually carry information in this object
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A character vector of rank names, in taxonomic order.
#' @noRd
vpf_available_ranks <- function(tse) {
  if (is.null(tse) || nrow(tse) == 0) {
    return(character(0))
  }
  rd <- SummarizedExperiment::rowData(tse)
  present <- intersect(VPF_TAXONOMY_RANKS, colnames(rd))
  keep <- vapply(present, function(r) any(!is.na(vpf_blank_na(rd[[r]]))), logical(1))
  present[keep]
}

#' Row-data table with blanks normalized to `NA` in taxonomic columns
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A `data.frame`, with a `Contig` column holding the row names.
#' @noRd
vpf_row_table <- function(tse) {
  if (is.null(tse)) {
    return(data.frame())
  }
  df <- as.data.frame(SummarizedExperiment::rowData(tse), optional = TRUE)
  if (nrow(df) > 0 || ncol(df) > 0) {
    for (r in intersect(VPF_TAXONOMY_RANKS, colnames(df))) {
      df[[r]] <- vpf_blank_na(df[[r]])
    }
  }
  ids <- rownames(tse)
  if (is.null(ids)) {
    ids <- as.character(seq_len(nrow(tse)))
  }
  cbind(data.frame(Contig = ids, stringsAsFactors = FALSE), df)
}

#' Sample identifiers, preferring a `sample_name` column over assay colnames
#' @noRd
vpf_sample_names <- function(tse) {
  if (is.null(tse) || ncol(tse) == 0) {
    return(character(0))
  }
  cd <- SummarizedExperiment::colData(tse)
  if ("sample_name" %in% colnames(cd)) {
    nm <- as.character(cd[["sample_name"]])
    if (!anyNA(nm) && !any(nm == "")) {
      return(nm)
    }
  }
  nm <- colnames(tse)
  if (is.null(nm)) as.character(seq_len(ncol(tse))) else nm
}


# -- Annotation availability -------------------------------------------------

#' Catalogue of annotation families the viewer can display
#'
#' Each entry records which ViroProfiler option produces the annotation, so the
#' UI can tell the user exactly what to switch on rather than showing an empty
#' panel.
#'
#' @return A list of descriptors.
#' @noRd
vpf_annotation_catalog <- function() {
  list(
    list(
      key = "taxonomy", family = "Taxonomy",
      columns = VPF_TAXONOMY_RANKS, match = "any",
      param = "--use_vitap (VITAP) plus vConTACT3; merged by bin/merge_taxonomy.py",
      note = "ICTV ranks assigned to each contig."
    ),
    list(
      key = "checkv", family = "CheckV quality", prefix = "checkv_",
      param = "always on (CHECKV process)",
      note = "Completeness, contamination, provirus status and quality tier."
    ),
    list(
      key = "genomad", family = "geNomad", prefix = "genomad_",
      param = "always on (GENOMAD process)",
      note = "Virus score, topology, hallmark gene count and geNomad taxonomy."
    ),
    list(
      key = "virsorter2", family = "VirSorter2", prefix = "virsorter2_",
      param = "runs on the selected candidates to feed DRAM-v",
      note = "Per-group viral scores and hallmark counts."
    ),
    list(
      key = "vibrant", family = "VIBRANT", prefix = "vibrant_",
      param = "--use_vibrant (default true)",
      note = "Independent viral quality call and replication cycle."
    ),
    list(
      key = "replicyc", family = "Replication cycle",
      columns = c("bacphlip_replicyc", "replidec_replicyc", "vibrant_replicyc"),
      match = "any",
      param = "--replicyc bacphlip | replidec",
      note = "Temperate vs virulent lifestyle prediction."
    ),
    list(
      key = "iphop", family = "Host prediction (iPHoP)", prefix = "iphop_",
      param = "--use_iphop (default true)",
      note = "Predicted host genus with score and contributing methods."
    ),
    list(
      key = "phist", family = "Host prediction (PHIST)", prefix = "phist_",
      param = "PHIST host prediction",
      note = "k-mer based host prediction against a bacterial genome set."
    ),
    list(
      key = "viral_votes", family = "Viral-evidence votes",
      columns = c("viral_vote_n", "viral_vote_evidence", "viral_selected"),
      match = "any",
      param = "recorded by annotate_viral_votes() when the object is built",
      note = paste("Which detector marked each contig as viral, and how many did.",
                   "The votes are combined with OR, so a contig can enter on one",
                   "piece of evidence alone.")
    ),
    list(
      key = "dramv", family = "DRAM-v (contig level)", prefix = "dramv_",
      param = "--use_dram (default true)",
      note = "Per-contig auxiliary metabolic gene counts."
    ),
    list(
      key = "checkamg", family = "CheckAMG (contig level)", prefix = "checkamg_",
      param = "--use_checkamg (default true)",
      note = paste("Per-contig counts of metabolic (AMG), physiological (APG) and",
                   "regulatory (AReG) auxiliary genes, with a high-confidence subset.")
    ),
    list(
      key = "pharokka", family = "Pharokka (contig level)", prefix = "pharokka_",
      param = "pharokka annotation",
      note = "Per-contig gene, CARD and VFDB hit counts."
    ),
    list(
      key = "binning", family = "Viral MAG binning",
      columns = c("vrhyme_bin", "phamb_bin", "bin_id"), match = "any",
      param = "--binning vrhyme | phamb",
      note = "Bin membership from vRhyme or PHAMB."
    ),
    list(
      key = "dvf", family = "DeepVirFinder score", prefix = "dvf_",
      param = "produced for PHAMB (bin/genomad_to_dvf.py)",
      note = "Per-contig virus score in DeepVirFinder layout."
    )
  )
}

#' Which annotation families are present in this object?
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A data frame with `key`, `family`, `available`, `n_columns`,
#'   `columns`, `param` and `note`.
#' @noRd
vpf_annotation_status <- function(tse) {
  cat_list <- vpf_annotation_catalog()
  cols <- if (is.null(tse)) character(0) else colnames(SummarizedExperiment::rowData(tse))
  rd <- if (is.null(tse)) NULL else SummarizedExperiment::rowData(tse)

  informative <- function(nm) {
    if (is.null(rd) || !nm %in% colnames(rd)) {
      return(FALSE)
    }
    v <- rd[[nm]]
    any(!is.na(vpf_blank_na(v)))
  }

  rows <- lapply(cat_list, function(entry) {
    found <- if (!is.null(entry$prefix)) {
      grep(paste0("^", entry$prefix), cols, value = TRUE)
    } else {
      intersect(entry$columns, cols)
    }
    found <- found[vapply(found, informative, logical(1))]
    data.frame(
      key = entry$key, family = entry$family,
      available = length(found) > 0, n_columns = length(found),
      columns = paste(found, collapse = ", "),
      param = entry$param, note = entry$note,
      stringsAsFactors = FALSE
    )
  })
  status <- do.call(rbind, rows)

  ga <- vpf_gene_annotations(tse)
  rbind(status, data.frame(
    key = "gene_annotations", family = "Gene-level annotations",
    available = !is.null(ga), n_columns = if (is.null(ga)) 0L else ncol(ga),
    columns = if (is.null(ga)) "" else paste(colnames(ga), collapse = ", "),
    param = "--use_dram (DRAM-v) and pharokka; stored in metadata(tse)$gene_annotations",
    note = "Per-gene DRAM-v, PHROG, CARD and VFDB annotations.",
    stringsAsFactors = FALSE
  ))
}

#' Gene-level annotation table, or `NULL` when the run did not produce one
#'
#' Prefers the merged `gene_annotations` table. When only the per-tool tables
#' are present the widest of them is returned, so a run that produced, say,
#' CheckAMG output but no merge step is still browsable.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return A `data.frame` or `NULL`.
#' @noRd
vpf_gene_annotations <- function(tse) {
  if (is.null(tse)) {
    return(NULL)
  }
  md <- S4Vectors::metadata(tse)
  if (is.null(md)) {
    return(NULL)
  }
  usable <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (nrow(x) == 0 || ncol(x) == 0) NULL else x
  }
  # `[[` rather than `$`: partial matching would make a missing
  # `gene_annotations` silently resolve to `gene_annotations_by_tool`, whose
  # value is a list of tables rather than a table.
  ga <- usable(md[["gene_annotations"]])
  if (!is.null(ga)) {
    return(ga)
  }
  by_tool <- md[["gene_annotations_by_tool"]]
  if (is.null(by_tool) || length(by_tool) == 0) {
    return(NULL)
  }
  tables <- Filter(Negate(is.null), lapply(by_tool, usable))
  if (length(tables) == 0) {
    return(NULL)
  }
  tables[[which.max(vapply(tables, nrow, numeric(1)))]]
}

#' Gene-level identifier column used by the annotation tables
#'
#' @param df A gene annotation table.
#' @return The column name holding contig identifiers, or `NULL`.
#' @noRd
vpf_gene_contig_column <- function(df) {
  if (is.null(df)) {
    return(NULL)
  }
  hit <- intersect(c("Contig", "contig", "contig_id", "scaffold"), colnames(df))
  if (length(hit) == 0) NULL else hit[[1]]
}

#' Per-contig count of viral-evidence votes, or `NULL`
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return An integer vector, or `NULL` when the object records no votes.
#' @noRd
vpf_viral_votes <- function(tse) {
  if (is.null(tse) || nrow(tse) == 0) {
    return(NULL)
  }
  rd <- SummarizedExperiment::rowData(tse)
  if (!"viral_vote_n" %in% colnames(rd)) {
    return(NULL)
  }
  v <- vpf_as_numeric_column(rd[["viral_vote_n"]])
  if (is.null(v)) NULL else as.integer(v)
}

#' Summary of how the viral subset was selected, or `NULL`
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @return The `metadata(tse)$viral_selection` list, or `NULL`.
#' @noRd
vpf_viral_selection <- function(tse) {
  if (is.null(tse)) {
    return(NULL)
  }
  S4Vectors::metadata(tse)[["viral_selection"]]
}


# -- Experimental design -----------------------------------------------------

#' Names that identify a sample rather than describe it
#' @noRd
VPF_ID_COLUMNS <- c("sample", "sample_name", "sample_id", "sampleid", "sample.id",
                    "id", "name", "barcode", "run", "library", "library_id",
                    "file", "filename", "fastq_1", "fastq_2")

#' Column-data variables that can serve as a grouping factor
#'
#' A sample-identifier column is never a grouping variable: it produces `n`
#' groups of size one, which makes every group test meaningless. Such columns
#' are recognized both by name and by having one distinct value per sample.
#'
#' The one-value-per-sample rule is not applied when there are only two
#' samples, because there any two-level variable necessarily has one sample per
#' level, and refusing it would make an uploaded metadata table useless on a
#' pilot run. Whether a *test* may be run on the resulting design is decided
#' separately by [vpf_design_check()], which requires at least two samples per
#' group regardless.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @param max_levels Highest number of levels still treated as a grouping.
#' @return A character vector of column names.
#' @noRd
vpf_group_candidates <- function(tse, max_levels = 12) {
  if (is.null(tse) || ncol(tse) == 0) {
    return(character(0))
  }
  cd <- as.data.frame(SummarizedExperiment::colData(tse), optional = TRUE)
  if (ncol(cd) == 0) {
    return(character(0))
  }
  n <- ncol(tse)
  keep <- vapply(colnames(cd), function(nm) {
    if (tolower(nm) %in% VPF_ID_COLUMNS) {
      return(FALSE)
    }
    v <- cd[[nm]]
    if (is.numeric(v) && length(unique(stats::na.omit(v))) > max_levels) {
      return(FALSE)
    }
    lv <- unique(stats::na.omit(vpf_blank_na(v)))
    if (length(lv) < 2 || length(lv) > max_levels) {
      return(FALSE)
    }
    if (n > 2 && length(lv) >= n) {
      return(FALSE)
    }
    TRUE
  }, logical(1))
  colnames(cd)[keep]
}

#' Number of distinguishable group-label permutations
#'
#' Used to decide whether a permutation test can reach a given p-value at all.
#' `P = n! / prod(n_j!)` counts labelled allocations; dividing by the
#' factorials of the counts of equally sized groups removes relabelings of the
#' same partition.
#'
#' @param sizes Integer vector of group sizes.
#' @return A numeric scalar.
#' @noRd
vpf_permutation_count <- function(sizes) {
  sizes <- sizes[sizes > 0]
  if (length(sizes) < 2) {
    return(1)
  }
  n <- sum(sizes)
  log_p <- lfactorial(n) - sum(lfactorial(sizes))
  same_size <- table(sizes)
  log_p <- log_p - sum(lfactorial(as.numeric(same_size)))
  exp(log_p)
}

#' Evaluate which statistics the current design can support
#'
#' The thresholds follow from the mathematics of each procedure rather than
#' from taste:
#' * classical MDS has at most `n - 1` axes, so a two-dimensional PCoA needs
#'   `n >= 3`;
#' * a two-dimensional NMDS configuration is saturated for `n <= 5`, where
#'   stress is trivially near zero;
#' * a permutation test cannot reach `p < 0.05` unless more than 20
#'   distinguishable permutations exist.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @param group Name of a `colData` column, or `NULL`.
#' @param assay Assay used for the abundance-derived checks, or `NULL`.
#' @return A list of design facts and boolean guards.
#' @noRd
vpf_design_check <- function(tse, group = NULL, assay = NULL) {
  out <- list(
    n_samples = 0L, n_features = 0L, group = NULL, n_groups = 0L,
    group_sizes = integer(0), min_group_size = 0L,
    has_group = FALSE, n_permutations = 1,
    can_alpha = FALSE, can_chao1 = FALSE, can_group_test = FALSE,
    kw_asymptotic_ok = FALSE, wilcoxon_can_reach_05 = FALSE,
    can_pairwise_beta = FALSE, can_pcoa = FALSE, pcoa_saturated = FALSE,
    can_nmds = FALSE, can_permanova = FALSE, permanova_resolution_ok = FALSE,
    can_differential = FALSE, empty_libraries = character(0)
  )
  if (is.null(tse)) {
    return(out)
  }
  out$n_samples <- ncol(tse)
  out$n_features <- nrow(tse)
  n <- out$n_samples

  if (!is.null(assay) && assay %in% SummarizedExperiment::assayNames(tse) && nrow(tse) > 0) {
    libs <- colSums(SummarizedExperiment::assay(tse, assay), na.rm = TRUE)
    out$library_sizes <- libs
    out$empty_libraries <- vpf_sample_names(tse)[!is.finite(libs) | libs <= 0]
    out$can_alpha <- nrow(tse) > 0 && length(out$empty_libraries) < n
    out$can_chao1 <- out$can_alpha && vpf_is_count_assay(tse, assay)
  } else {
    out$can_alpha <- nrow(tse) > 0 && n >= 1
  }

  if (!is.null(group) && nzchar(group) &&
      group %in% colnames(SummarizedExperiment::colData(tse))) {
    g <- vpf_blank_na(SummarizedExperiment::colData(tse)[[group]])
    sizes <- table(g[!is.na(g)])
    out$group <- group
    out$n_groups <- length(sizes)
    out$group_sizes <- stats::setNames(as.integer(sizes), names(sizes))
    out$min_group_size <- if (length(sizes)) min(as.integer(sizes)) else 0L
    out$has_group <- out$n_groups >= 2
    out$n_permutations <- vpf_permutation_count(as.integer(sizes))
  }

  g <- out$n_groups
  min_g <- out$min_group_size

  out$can_group_test <- g >= 2 && min_g >= 2
  out$kw_asymptotic_ok <- out$can_group_test && min_g >= 5
  if (g == 2) {
    sz <- as.integer(out$group_sizes)
    out$wilcoxon_can_reach_05 <- choose(sum(sz), sz[[1]]) > 40
  }
  out$can_pairwise_beta <- n >= 2 && out$n_features > 0
  out$can_pcoa <- n >= 3 && out$n_features > 0
  out$pcoa_saturated <- n == 3
  out$can_nmds <- n >= 6 && out$n_features > 0
  out$can_permanova <- n >= 4 && g >= 2 && min_g >= 2
  out$permanova_resolution_ok <- out$can_permanova && out$n_permutations > 20
  out$can_differential <- g >= 2 && min_g >= 2
  out
}


# -- Data operations ---------------------------------------------------------

#' Prevalence of each contig across samples
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @param assay Assay name.
#' @param detection Value above which a contig counts as present.
#' @return A numeric vector of proportions, one per contig.
#' @noRd
vpf_prevalence <- function(tse, assay, detection = 0) {
  if (is.null(tse) || nrow(tse) == 0 || ncol(tse) == 0) {
    return(numeric(0))
  }
  if (!assay %in% SummarizedExperiment::assayNames(tse)) {
    return(numeric(0))
  }
  mat <- SummarizedExperiment::assay(tse, assay)
  present <- !is.na(mat) & mat > detection
  rowSums(present) / ncol(mat)
}

#' Subset rows of a TSE with a logical vector that may contain `NA`
#'
#' `NA` in a logical subscript is an error for `SummarizedExperiment`, and
#' missing annotations are common, so missing values are resolved explicitly.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @param keep Logical vector, possibly containing `NA`.
#' @param na_keeps Value substituted for `NA`.
#' @return The subset object.
#' @noRd
vpf_subset_rows <- function(tse, keep, na_keeps = FALSE) {
  if (is.null(tse)) {
    return(NULL)
  }
  if (is.null(keep) || length(keep) == 0) {
    return(tse)
  }
  if (length(keep) != nrow(tse)) {
    return(tse)
  }
  keep[is.na(keep)] <- isTRUE(na_keeps)
  tse[keep, ]
}

#' Coerce a row-data column to numeric without lexicographic surprises
#'
#' A numeric threshold applied to a factor compares level codes, and applied to
#' a character column compares strings, both of which silently select the wrong
#' contigs.
#'
#' @param x A vector.
#' @return A numeric vector, or `NULL` if the column cannot be coerced.
#' @noRd
vpf_as_numeric_column <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  if (is.factor(x)) {
    x <- as.character(x)
  }
  if (is.logical(x) && all(is.na(x))) {
    return(NULL)
  }
  suppressWarnings(num <- as.numeric(as.character(x)))
  if (all(is.na(num)) && any(!is.na(x))) {
    return(NULL)
  }
  num
}

#' Read a TSE from disk with a helpful error message
#'
#' @param path File path.
#' @return A list with `tse` and `error`; exactly one is non-`NULL`.
#' @noRd
vpf_read_tse <- function(path) {
  if (is.null(path) || !nzchar(path)) {
    return(list(tse = NULL, error = "No file selected."))
  }
  if (!file.exists(path)) {
    return(list(tse = NULL, error = paste0("File not found: ", path)))
  }
  obj <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(obj, "error")) {
    return(list(tse = NULL, error = paste0("Could not read the file: ", conditionMessage(obj))))
  }
  if (!methods::is(obj, "SummarizedExperiment")) {
    return(list(tse = NULL, error = paste0(
      "The file contains a '", paste(class(obj), collapse = "/"),
      "' object. A TreeSummarizedExperiment produced by create_vpftse() is required."
    )))
  }
  if (length(SummarizedExperiment::assayNames(obj)) == 0) {
    return(list(tse = NULL, error = "The object contains no assays."))
  }
  # Objects written before the assay was renamed carry CoverM's trimmed mean
  # under the name `tmm`, which reads as edgeR's unrelated TMM normalization.
  # Canonicalizing here means every module downstream sees one spelling.
  before <- SummarizedExperiment::assayNames(obj)
  obj <- tryCatch(normalize_assay_names(obj, quiet = TRUE),
                  error = function(e) obj)
  after <- SummarizedExperiment::assayNames(obj)
  renamed <- if (length(before) == length(after)) {
    stats::setNames(after[before != after], before[before != after])
  } else {
    character(0)
  }
  list(tse = obj, error = NULL, renamed_assays = renamed)
}

#' Directories that may hold the bundled demo datasets
#'
#' `inst/extdata` is where they live once the package is installed;
#' `data-raw` is where they live in a source checkout, which is what
#' `pkgload::load_all()` sees.
#'
#' @return A character vector of existing directories.
#' @noRd
vpf_demo_dirs <- function() {
  roots <- c(
    app_sys("extdata"),
    # Source checkout: app_sys("app/www") resolves to <pkg>/inst/app/www.
    {
      www <- app_sys("app/www")
      if (nzchar(www)) {
        c(file.path(dirname(dirname(dirname(www))), "data-raw"),
          file.path(dirname(dirname(www)), "extdata"))
      } else {
        character(0)
      }
    },
    file.path(getwd(), "data-raw"),
    file.path(getwd(), "inst", "extdata")
  )
  roots <- roots[nzchar(roots)]
  unique(roots[dir.exists(roots)])
}

#' Built-in demo datasets shipped with the package
#'
#' @return A named list of file paths that exist on this installation.
#' @noRd
vpf_demo_datasets <- function() {
  candidates <- list(
    "Gut virome, 60 contigs x 8 samples (Healthy vs Disease)" = "test_viroprofiler.rds",
    "Environmental virome, 40 contigs x 6 samples (Ocean vs Soil)" = "test_viroprofiler_2.rds"
  )
  dirs <- vpf_demo_dirs()
  out <- list()
  for (nm in names(candidates)) {
    hits <- file.path(dirs, candidates[[nm]])
    hits <- hits[file.exists(hits)]
    if (length(hits) > 0) {
      out[[nm]] <- normalizePath(hits[[1]])
    }
  }
  out
}

#' Attach an external sample-metadata table to `colData`
#'
#' ViroProfiler's samplesheet is restricted to `sample`, `fastq_1` and
#' `fastq_2`, so study metadata can never arrive with the pipeline output.
#' Joining it here is what makes every grouped analysis possible.
#'
#' @param tse A `TreeSummarizedExperiment`.
#' @param meta A data frame whose first column, or a column named
#'   `sample`/`sample_name`/`sample_id`, holds the sample identifiers.
#' @return A list with `tse`, `matched`, `unmatched_samples`, `unused_rows`
#'   and `error`.
#' @noRd
vpf_join_sample_metadata <- function(tse, meta) {
  if (is.null(tse)) {
    return(list(tse = NULL, matched = 0L, unmatched_samples = character(0),
                unused_rows = character(0), error = "No dataset loaded."))
  }
  if (is.null(meta) || nrow(meta) == 0) {
    return(list(tse = tse, matched = 0L, unmatched_samples = character(0),
                unused_rows = character(0), error = "The metadata table is empty."))
  }
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  id_col <- intersect(c("sample", "sample_name", "sample_id", "Sample", "SampleID", "sample.id"),
                      colnames(meta))
  id_col <- if (length(id_col) > 0) id_col[[1]] else colnames(meta)[[1]]
  keys <- trimws(as.character(meta[[id_col]]))

  targets <- vpf_sample_names(tse)
  alt <- colnames(tse)
  idx <- match(targets, keys)
  if (anyNA(idx) && !is.null(alt)) {
    idx2 <- match(alt, keys)
    if (sum(!is.na(idx2)) > sum(!is.na(idx))) {
      idx <- idx2
    }
  }
  if (all(is.na(idx))) {
    return(list(
      tse = tse, matched = 0L, unmatched_samples = targets,
      unused_rows = keys,
      error = paste0(
        "None of the sample identifiers matched. The dataset uses: ",
        paste(utils::head(targets, 5), collapse = ", "),
        if (length(targets) > 5) ", ..." else "",
        "; the uploaded table's '", id_col, "' column holds: ",
        paste(utils::head(keys, 5), collapse = ", "),
        if (length(keys) > 5) ", ..." else "", "."
      )
    ))
  }

  add <- meta[idx, setdiff(colnames(meta), id_col), drop = FALSE]
  rownames(add) <- colnames(tse)
  cd <- SummarizedExperiment::colData(tse)
  clash <- intersect(colnames(add), colnames(cd))
  if (length(clash) > 0) {
    colnames(add)[colnames(add) %in% clash] <- paste0(clash, "_meta")
  }
  SummarizedExperiment::colData(tse) <- cbind(cd, S4Vectors::DataFrame(add))
  list(
    tse = tse, matched = sum(!is.na(idx)),
    unmatched_samples = targets[is.na(idx)],
    unused_rows = setdiff(keys, targets),
    error = NULL
  )
}

#' Read a delimited or Excel metadata table
#'
#' @param path File path.
#' @param name Original file name, used to pick a parser.
#' @return A list with `data` and `error`.
#' @noRd
vpf_read_metadata_file <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  out <- tryCatch({
    if (ext %in% c("xlsx", "xls")) {
      openxlsx::read.xlsx(path)
    } else {
      sep <- if (ext == "csv") "," else if (ext %in% c("tsv", "txt")) "\t" else "auto"
      as.data.frame(data.table::fread(path, sep = sep, data.table = FALSE))
    }
  }, error = function(e) e)
  if (inherits(out, "error")) {
    return(list(data = NULL, error = paste0("Could not parse the table: ", conditionMessage(out))))
  }
  if (is.null(out) || ncol(out) < 2) {
    return(list(data = NULL, error = "The table needs a sample-identifier column and at least one metadata column."))
  }
  list(data = out, error = NULL)
}


# -- Plot helpers ------------------------------------------------------------

#' A plotly placeholder carrying an explanatory message
#'
#' Returning `NULL` from `renderPlotly()` leaves a blank rectangle with no
#' indication of what happened, which is how several degenerate cases used to
#' present themselves.
#'
#' @param message Text to display.
#' @return A plotly object.
#' @noRd
vpf_message_plot <- function(message) {
  plotly::layout(
    plotly::plot_ly(type = "scatter", mode = "markers", x = 0, y = 0,
                    marker = list(opacity = 0), hoverinfo = "none"),
    xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
    annotations = list(list(
      text = paste(strwrap(message, width = 60), collapse = "<br>"),
      showarrow = FALSE, x = 0.5, y = 0.5, xref = "paper", yref = "paper",
      font = list(size = 14, color = "#555555"), align = "center"
    ))
  )
}

#' Convert a ggplot to plotly, degrading to a message on failure
#' @noRd
vpf_ggplotly <- function(p, ...) {
  out <- tryCatch(plotly::ggplotly(p, ...), error = function(e) e)
  if (inherits(out, "error")) {
    return(vpf_message_plot(paste("Could not render the plot:", conditionMessage(out))))
  }
  out
}

#' Prepare a data frame for display in a table
#'
#' Rounds numeric columns and drops row names. Row names matter here because
#' `reactable()` renders any non-automatic row names as an extra unlabelled
#' first column, and subsetting a data frame is enough to make them
#' non-automatic, so a filtered table would silently grow a duplicate of its
#' key column.
#'
#' @param df A data frame.
#' @param digits Digits kept in numeric columns.
#' @return The data frame, ready to hand to `reactable()`.
#' @noRd
vpf_round_df <- function(df, digits = 4) {
  if (is.null(df) || ncol(df) == 0) {
    return(df)
  }
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], function(x) round(x, digits))
  rownames(df) <- NULL
  df
}

#' A discrete colour palette that stays readable for many taxa
#' @noRd
vpf_palette <- function(n) {
  base <- c(
    "#4E79A7", "#F28E2B", "#59A14F", "#E15759", "#B07AA1", "#76B7B2",
    "#EDC948", "#FF9DA7", "#9C755F", "#BAB0AC", "#86BCB6", "#D37295",
    "#A0CBE8", "#FFBE7D", "#8CD17D", "#FF9D9A", "#D4A6C8", "#B6992D"
  )
  if (n <= length(base)) {
    return(base[seq_len(n)])
  }
  grDevices::colorRampPalette(base)(n)
}
