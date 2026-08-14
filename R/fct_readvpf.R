####################################################################################
## Internal helpers
##
## Every reader in this file joins its result onto the contig identifiers of the
## abundance matrices. The helpers below exist because the ways that join can fail
## are all silent: a renamed column, a decorated sequence name, an empty-string
## rank read as a value, a tool table with two rows for one contig. None of those
## raise an error on their own; they produce an object that is the wrong size or
## carries the wrong values.
####################################################################################


#' Check that an input file exists and is not zero bytes
#'
#' A zero-byte file is a failed upstream process, not an empty result: a tool
#' that ran and found nothing still writes its header line. Distinguishing the
#' two here keeps the error message at the file that is actually broken.
#'
#' @param fin File path
#' @param what Tool name used in the error message
#' @return `fin`, invisibly
#' @noRd
.check_input_file <- function(fin, what) {
  if (is.null(fin) || length(fin) != 1L || is.na(fin) || !nzchar(fin)) {
    stop(what, " file path is empty or NA", call. = FALSE)
  }
  if (!file.exists(fin)) stop(what, " file not found: ", fin, call. = FALSE)
  if (dir.exists(fin)) stop(what, " path is a directory, not a file: ", fin, call. = FALSE)
  if (isTRUE(file.size(fin) == 0)) {
    stop(what, " file is empty (0 bytes), so it has no header line: ", fin,
         ". A tool that ran and found nothing still writes a header; treat this as a failed run.",
         call. = FALSE)
  }
  invisible(fin)
}


#' Decide whether an optional file argument means "absent"
#'
#' `argparser` fills unsupplied arguments with `NA`, Nextflow can interpolate an
#' empty string, and a shell can pass the literal `"null"`. All three reach
#' `create_vpftse()` as a value that is not `NULL`, and the old code only guarded
#' `fin_genomad` against one of them.
#'
#' @param x Candidate file path
#' @return TRUE when the argument should be treated as not supplied
#' @noRd
.is_absent_path <- function(x) {
  if (is.null(x)) return(TRUE)
  if (length(x) == 0L) return(TRUE)
  if (length(x) > 1L) return(FALSE)
  if (is.na(x)) return(TRUE)
  x <- trimws(as.character(x))
  if (!nzchar(x)) return(TRUE)
  tolower(x) %in% c("null", "na", "none", "false")
}


#' Normalize an optional file argument
#'
#' @param x Candidate file path
#' @param what Tool name used in messages
#' @param required_exists Warn (and drop) when a non-blank path does not exist
#' @return A single file path, or NULL
#' @noRd
.optional_path <- function(x, what, required_exists = TRUE) {
  if (.is_absent_path(x)) return(NULL)
  x <- trimws(as.character(x))
  if (required_exists && !file.exists(x)) {
    warning(what, " file was requested but does not exist, so it is skipped: ", x,
            call. = FALSE)
    return(NULL)
  }
  x
}


#' Read a delimited table with the missing-value conventions these tools use
#'
#' Every tool here writes unassigned fields either as the literal `NA` or as an
#' empty field. `data.table::fread()` maps the first to `NA` but keeps the second
#' as `""`, which then survives every downstream `is.na()` test. geNomad's
#' all-`NA` `fdr` column is read as `logical` for the same reason, so a column
#' that is numeric in one run is logical in the next.
#'
#' @param fin File path
#' @param what Tool name used in messages
#' @param ... Passed to [data.table::fread()]
#' @return data.frame
#' @noRd
.read_table <- function(fin, what, ...) {
  .check_input_file(fin, what)
  df <- tryCatch(
    data.table::fread(fin, na.strings = c("NA", "", "NaN"), data.table = FALSE, ...),
    error = function(e) stop("Could not parse ", what, " file ", fin, ": ",
                             conditionMessage(e), call. = FALSE)
  )
  if (ncol(df) == 0L) {
    stop(what, " file has no columns: ", fin, call. = FALSE)
  }
  as.data.frame(df, stringsAsFactors = FALSE)
}


#' Fail with a readable message when a table is missing required columns
#'
#' @param df data.frame
#' @param required Character vector of required column names
#' @param what Tool name used in the error message
#' @return `df`, invisibly
#' @noRd
.require_cols <- function(df, required, what) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(what, " file missing columns: ", paste(missing, collapse = ", "),
         ". Found: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  invisible(df)
}


#' Return a column if it exists, otherwise a vector of NA of the right length
#'
#' Several readers reached for columns such as `fdr` or `card_match` with
#' `.data$name`, which raises an error when the column is absent. Tool output
#' varies with the flags it was run under, so absence is normal.
#'
#' @param df data.frame
#' @param name Column name
#' @param fill Value used when the column is absent
#' @return A vector of length `nrow(df)`
#' @noRd
.col_or_na <- function(df, name, fill = NA) {
  if (name %in% colnames(df)) return(df[[name]])
  rep(fill, nrow(df))
}


#' Strip the decorations tools add to sequence names
#'
#' Each flag corresponds to one tool's renaming convention, and each reader turns
#' on only the ones its own tool applies. Turning them all on everywhere would
#' silently truncate a legitimate contig name that happens to contain the
#' pattern.
#'
#' @param x Character vector of sequence identifiers
#' @param vs2_suffix Strip VirSorter2's `||full`, `||<i>_partial`, `||lt2gene`
#' @param genomad_provirus Strip geNomad's `|provirus_<start>_<end>` decoration
#' @param vibrant_fragment Strip VIBRANT's `_fragment_N` prophage suffix
#' @param cat_suffix Strip the `-cat_N` VirSorter2 category suffix that
#'   `--prep-for-dramv` appends, and that ViroProfiler v1 carried throughout
#' @param dramv_suffix Strip the `__full` / `__<i>_partial` / `__lt2gene` form
#'   that `--prep-for-dramv` writes in place of VirSorter2's `||` decoration
#' @param what Tool name used in the message emitted when identifiers change
#' @return Character vector
#' @noRd
.normalize_contig_ids <- function(x, vs2_suffix = FALSE, genomad_provirus = FALSE,
                                  vibrant_fragment = FALSE, cat_suffix = FALSE,
                                  dramv_suffix = FALSE, what = "input") {
  original <- as.character(x)
  out <- trimws(original)
  if (vs2_suffix) out <- sub("\\|\\|.*$", "", out)
  if (genomad_provirus) out <- sub("\\|provirus_[0-9]+.*$", "", out)
  if (vibrant_fragment) out <- sub("_fragment_[0-9]+$", "", out)
  if (cat_suffix) out <- sub("-cat_[0-6]$", "", out)
  # `--prep-for-dramv` replaces the `||` in a VirSorter2 name with `__` and
  # appends `-cat_N`, so the -cat_N suffix has to come off first. Only the three
  # literal tokens are matched, because contig names in this pipeline already
  # contain `__` (for example `HT02__NODE_1_length_21363_cov_7.4`).
  if (dramv_suffix) out <- sub("__(full|[0-9]+_partial|lt2gene)$", "", out)
  n_changed <- sum(!is.na(out) & !is.na(original) & out != original)
  if (n_changed > 0) {
    message(what, ": normalized ", n_changed,
            " sequence identifier(s) by stripping tool-specific decorations.")
  }
  out
}


#' Replace empty strings with NA in character columns
#'
#' `taxonomy_tse.tsv` leaves unresolved ranks empty rather than writing `NA`, so
#' 15 of 16 contigs in the reference run carry `Species == ""`. An empty string
#' passes every `is.na()` test, which makes an unassigned rank look assigned.
#'
#' @param df data.frame
#' @param cols Columns to clean; default all character columns
#' @return data.frame
#' @noRd
.blank_to_na <- function(df, cols = NULL) {
  if (is.null(cols)) cols <- colnames(df)[vapply(df, is.character, logical(1))]
  for (cn in intersect(cols, colnames(df))) {
    v <- df[[cn]]
    v[!is.na(v) & !nzchar(trimws(v))] <- NA_character_
    df[[cn]] <- v
  }
  df
}


#' Reduce a tool table to one row per contig
#'
#' VIBRANT emits one row per prophage fragment, iPHoP one row per candidate host
#' genus, CheckV one row per provirus. Joining any of those onto the abundance
#' matrix multiplies rows; the multiplication then either aborts in
#' `column_to_rownames()` or, if the duplicates were made unique, produces a
#' rowData longer than the assay.
#'
#' @param df data.frame with a `Contig` column
#' @param what Tool name used in the message
#' @param order_by Optional column; rows are sorted by it before the first row of
#'   each contig is kept, so the retained row is the best one rather than an
#'   arbitrary one
#' @param decreasing Sort `order_by` in decreasing order
#' @param severity `"message"` when the tool is documented to report several rows
#'   per contig, so the collapse is expected; `"warning"` when it is not
#' @return data.frame with unique `Contig`
#' @noRd
.dedup_by_contig <- function(df, what, order_by = NULL, decreasing = FALSE,
                             severity = c("warning", "message")) {
  severity <- match.arg(severity)
  if (!"Contig" %in% colnames(df)) return(df)
  if (!anyDuplicated(df$Contig)) return(df)
  n_before <- nrow(df)
  if (!is.null(order_by) && order_by %in% colnames(df)) {
    key <- df[[order_by]]
    if (is.factor(key)) key <- as.integer(key)
    df <- df[order(key, decreasing = decreasing, na.last = TRUE), , drop = FALSE]
  }
  df <- df[!duplicated(df$Contig), , drop = FALSE]
  txt <- paste0(what, ": collapsed ", n_before, " rows to ", nrow(df),
                " unique contigs. Multiple rows per contig were reduced to one",
                if (!is.null(order_by)) paste0(", keeping the best '", order_by, "'") else
                  ", keeping the first", ".")
  if (severity == "warning") warning(txt, call. = FALSE) else message(txt)
  rownames(df) <- NULL
  df
}


#' Left-join a tool table onto the feature annotation, refusing to do it silently
#'
#' `dplyr::left_join()` repairs an overlapping column name by appending `.x` and
#' `.y`, which changes the rowData schema without any signal. It also multiplies
#' rows when the right-hand table has duplicate keys. Both are checked here.
#'
#' @param x Feature annotation data.frame with a `Contig` column
#' @param y Tool table with a `Contig` column
#' @param what Tool name used in messages
#' @return data.frame with the same number of rows as `x`
#' @noRd
.join_tool_table <- function(x, y, what) {
  if (is.null(y)) return(x)
  if (!"Contig" %in% colnames(y)) {
    stop(what, " table has no 'Contig' column, so it cannot be joined", call. = FALSE)
  }
  if (nrow(y) == 0L) {
    # A tool that ran and reported nothing must still contribute its columns, or
    # the rowData schema changes between runs and every consumer that indexes
    # those columns has to guess whether absence means "not run" or "found
    # nothing". Downstream code can tell the two apart from
    # metadata(tse)$viroprofiler$join_match.
    message(what, ": table has no rows; its columns are added as all-NA.")
    for (cn in setdiff(colnames(y), "Contig")) {
      x[[cn]] <- y[[cn]][NA_integer_][rep(1L, nrow(x))]
    }
    attr(x, "vpf_match") <- c(n_tool = 0, n_matched = 0, n_features = nrow(x))
    return(x)
  }
  y <- .dedup_by_contig(y, what)
  collisions <- setdiff(intersect(colnames(x), colnames(y)), "Contig")
  if (length(collisions) > 0) {
    stop(what, ": column name collision with columns already present: ",
         paste(collisions, collapse = ", "),
         ". Prefix the tool's columns instead of letting dplyr append .x/.y.",
         call. = FALSE)
  }
  matched <- sum(x$Contig %in% y$Contig)
  if (matched == 0L) {
    warning(what, ": none of its ", nrow(y), " contig identifiers match the ",
            nrow(x), " contigs in the abundance matrix. Every ", what,
            " column will be NA. This is almost always an identifier mismatch, ",
            "not a biological result. Example tool ID: '", y$Contig[1],
            "'; example abundance ID: '", x$Contig[1], "'.", call. = FALSE)
  }
  out <- dplyr::left_join(x, y, by = "Contig")
  if (nrow(out) != nrow(x)) {
    stop(what, ": join changed the row count from ", nrow(x), " to ", nrow(out),
         call. = FALSE)
  }
  attr(out, "vpf_match") <- c(
    n_tool = nrow(y), n_matched = matched, n_features = nrow(x)
  )
  out
}


#' Validate one abundance matrix
#'
#' @param m matrix
#' @param label Assay name used in messages
#' @return matrix
#' @noRd
.validate_assay <- function(m, label) {
  if (!is.matrix(m)) stop("Assay '", label, "' is not a matrix", call. = FALSE)
  if (!is.numeric(m)) {
    stop("Assay '", label, "' is not numeric. A single non-numeric column in the ",
         "CoverM table coerces the whole matrix to character, and a character ",
         "assay is stored without complaint.", call. = FALSE)
  }
  if (is.null(rownames(m)) || anyNA(rownames(m)) || !all(nzchar(rownames(m)))) {
    stop("Assay '", label, "' has missing or empty contig names", call. = FALSE)
  }
  if (is.null(colnames(m)) || anyNA(colnames(m)) || !all(nzchar(colnames(m)))) {
    stop("Assay '", label, "' has missing or empty sample names", call. = FALSE)
  }
  if (anyDuplicated(rownames(m))) {
    stop("Assay '", label, "' has duplicate contig names: ",
         paste(unique(rownames(m)[duplicated(rownames(m))]), collapse = ", "),
         call. = FALSE)
  }
  if (anyDuplicated(colnames(m))) {
    stop("Assay '", label, "' has duplicate sample names", call. = FALSE)
  }
  m
}


#' Align a set of assay matrices to a common contig and sample order
#'
#' `TreeSummarizedExperiment()` accepts matrices positionally whenever their
#' dimnames can be made to look compatible. Four CoverM invocations write four
#' files, and nothing guarantees that they agree on row order.
#'
#' @param assay_list Named list of matrices
#' @return The list, reordered to the first matrix
#' @noRd
.align_assays <- function(assay_list) {
  labels <- names(assay_list)
  for (i in seq_along(assay_list)) {
    assay_list[[i]] <- .validate_assay(assay_list[[i]], labels[i])
  }
  ref_rows <- rownames(assay_list[[1]])
  ref_cols <- colnames(assay_list[[1]])
  for (i in seq_along(assay_list)[-1]) {
    m <- assay_list[[i]]
    if (!setequal(rownames(m), ref_rows)) {
      stop("Assay '", labels[i], "' covers a different contig set than '", labels[1],
           "'. Missing: ", paste(utils::head(setdiff(ref_rows, rownames(m)), 5), collapse = ", "),
           "; extra: ", paste(utils::head(setdiff(rownames(m), ref_rows), 5), collapse = ", "),
           call. = FALSE)
    }
    if (!setequal(colnames(m), ref_cols)) {
      stop("Assay '", labels[i], "' covers a different sample set than '", labels[1],
           "'. Missing: ", paste(setdiff(ref_cols, colnames(m)), collapse = ", "),
           "; extra: ", paste(setdiff(colnames(m), ref_cols), collapse = ", "),
           call. = FALSE)
    }
    if (!identical(rownames(m), ref_rows) || !identical(colnames(m), ref_cols)) {
      message("Assay '", labels[i], "' was reordered to match '", labels[1], "'.")
      assay_list[[i]] <- m[ref_rows, ref_cols, drop = FALSE]
    }
  }
  assay_list
}


#' Canonical assay names and what each one measures
#'


#' Rename legacy assays in an existing TSE
#'
#' Objects saved by earlier versions of vpfkit carry an assay called `tmm` that
#' holds CoverM trimmed-mean coverage depth. Reading those objects must keep
#' working; writing that name must not.
#'
#' @param tse A `SummarizedExperiment` or `TreeSummarizedExperiment`
#' @param quiet Suppress the warning describing what was renamed
#' @return The object with canonical assay names
#' @export
#' @examples
#' \dontrun{
#' tse <- readRDS("viroprofiler_output.rds")
#' tse <- normalize_assay_names(tse)
#' SummarizedExperiment::assayNames(tse)
#' }
normalize_assay_names <- function(tse, quiet = FALSE) {
  nms <- SummarizedExperiment::assayNames(tse)
  if (is.null(nms)) return(tse)
  renamed <- character(0)
  aliases <- vpf_legacy_assay_renames()
  for (old in names(aliases)) {
    new <- aliases[[old]]
    if (old %in% nms && !new %in% nms) {
      nms[nms == old] <- new
      renamed <- c(renamed, paste0(old, " -> ", new))
    }
  }
  if (length(renamed) > 0) {
    SummarizedExperiment::assayNames(tse) <- nms
    if (!quiet) {
      warning("Renamed legacy assay(s): ", paste(renamed, collapse = "; "),
              ". CoverM 'trimmed_mean' is the trimmed mean of per-base coverage depth; ",
              "it is not edgeR's TMM normalization, and the old name implied that it was.",
              call. = FALSE)
    }
  }
  tse
}


#' Read a sample metadata table from CSV, TSV or XLSX
#'
#' @param fin File path
#' @return data.frame
#' @noRd
.read_metadata_table <- function(fin) {
  .check_input_file(fin, "Sample metadata")
  if (grepl("\\.xlsx?$", fin, ignore.case = TRUE)) {
    df <- openxlsx::read.xlsx(fin)
  } else {
    df <- data.table::fread(fin, na.strings = c("NA", "", "NaN"), data.table = FALSE)
  }
  as.data.frame(df, stringsAsFactors = FALSE)
}


#' Find the column of a metadata table that holds sample identifiers
#'
#' @param df data.frame
#' @return Column name
#' @noRd
.detect_sample_id_col <- function(df) {
  candidates <- c("sample_id", "sample", "sample_name", "sampleid", "sampleID",
                  "SampleID", "Sample", "Sample_ID", "sample_ID", "id", "ID", "name")
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  message("Sample metadata has no recognizable sample identifier column; ",
          "using the first column '", colnames(df)[1], "'.")
  colnames(df)[1]
}


#' Build colData, aligned to the sample order of the assays
#'
#' A metadata table in a different order than the assay columns is the one
#' failure that produces a plausible object with every phenotype attached to the
#' wrong sample, so samples are matched by name and never by position.
#'
#' @param samples Sample names, in assay column order
#' @param df_metadata Optional metadata data.frame
#' @return data.frame with rownames equal to `samples`
#' @noRd
.build_coldata <- function(samples, df_metadata = NULL) {
  base <- data.frame(sample_name = samples, row.names = samples,
                     stringsAsFactors = FALSE)
  if (is.null(df_metadata)) return(base)

  df_metadata <- as.data.frame(df_metadata, stringsAsFactors = FALSE)
  if (nrow(df_metadata) == 0L) {
    warning("Sample metadata has no rows; colData falls back to sample names only.",
            call. = FALSE)
    return(base)
  }
  id_col <- .detect_sample_id_col(df_metadata)
  ids <- trimws(as.character(df_metadata[[id_col]]))
  if (anyDuplicated(ids)) {
    stop("Sample metadata has duplicate identifiers in column '", id_col, "': ",
         paste(unique(ids[duplicated(ids)]), collapse = ", "), call. = FALSE)
  }

  missing_meta <- setdiff(samples, ids)
  extra_meta <- setdiff(ids, samples)
  if (length(missing_meta) > 0) {
    warning("Sample metadata has no row for ", length(missing_meta), " of ",
            length(samples), " samples: ", paste(missing_meta, collapse = ", "),
            ". Their metadata columns are NA.", call. = FALSE)
  }
  if (length(extra_meta) > 0) {
    warning("Sample metadata has ", length(extra_meta),
            " row(s) that match no sample in the abundance matrix and are dropped: ",
            paste(extra_meta, collapse = ", "),
            ". Check that the identifier column '", id_col,
            "' uses the same sample names as the abundance table.", call. = FALSE)
  }
  if (length(missing_meta) == length(samples)) {
    warning("Not one metadata identifier matches a sample name; colData falls back ",
            "to sample names only.", call. = FALSE)
    return(base)
  }

  idx <- match(samples, ids)
  aligned <- df_metadata[idx, setdiff(colnames(df_metadata), id_col), drop = FALSE]
  rownames(aligned) <- samples
  out <- cbind(base, aligned)
  clash <- intersect(colnames(base), colnames(aligned))
  if (length(clash) > 0) {
    out <- out[, !duplicated(colnames(out)), drop = FALSE]
  }
  out
}


#' Read the per-sample mapping statistics CoverM writes to stderr
#'
#' CoverM reports how many reads it mapped, and out of how many, for every sample.
#' ViroProfiler publishes those logs next to the abundance tables, and nothing else
#' in the output records the per-sample read counts at all. Without them a low
#' abundance cannot be told from a shallowly sequenced sample.
#'
#' @section What `n_reads_total` counts:
#' The reads presented to the mapper, which is what survived fastp and, when it
#' runs, decontamination. It is **not** the raw sequencing depth, and normalizing
#' by it is not the same as normalizing by library size. The raw counts are in
#' `fastp/<sample>.fastp.json`.
#'
#' @section Which log to read:
#' ViroProfiler runs CoverM once per method and writes one log each, and the
#' numbers differ between them: on the two-sample test set `log_contig_count.txt`
#' reports 4706 of 4706 reads mapped for HT02 while
#' `log_contig_trimmed_mean.txt` reports 4703, because coverage trimming discards
#' a few. Read the log belonging to the assay being interpreted; for library size
#' that is `log_contig_count.txt`, the invocation that produced the count matrix.
#'
#' @param fin A CoverM log file, e.g. `abundance/log_contig_count.txt`
#' @return data.frame with `sample_id`, `n_reads_mapped`, `n_reads_total`,
#'   `mapping_rate`
#' @export
#' @examples
#' \dontrun{
#' read_coverm_log("output/abundance/log_contig_count.txt")
#' }
read_coverm_log <- function(fin) {
  .check_input_file(fin, "CoverM log")
  lines <- readLines(fin, warn = FALSE)
  pat <- "In sample '(.+)', found ([0-9]+) reads mapped out of ([0-9]+) total"
  hit <- grep(pat, lines, value = TRUE)
  if (length(hit) == 0L) {
    warning("CoverM log has no 'In sample ... reads mapped' lines: ", fin,
            ". No sequencing depth information was extracted.", call. = FALSE)
    return(data.frame(sample_id = character(0), n_reads_mapped = numeric(0),
                      n_reads_total = numeric(0), mapping_rate = numeric(0),
                      stringsAsFactors = FALSE))
  }
  m <- regmatches(hit, regexec(pat, hit))
  out <- data.frame(
    sample_id      = vapply(m, `[`, character(1), 2),
    n_reads_mapped = as.numeric(vapply(m, `[`, character(1), 3)),
    n_reads_total  = as.numeric(vapply(m, `[`, character(1), 4)),
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$sample_id), , drop = FALSE]
  out$mapping_rate <- out$n_reads_mapped / out$n_reads_total
  rownames(out) <- NULL
  out
}


####################################################################################
## Readers
####################################################################################


#' Detect ViroProfiler pipeline version from contig IDs
#'
#' ViroProfiler v1 contigs have `-cat_N` suffixes from VirSorter2.
#'
#' @param contig_ids Character vector of contig identifiers
#' @return "v1" or "v2"
#' @noRd
.detect_vp_version <- function(contig_ids) {
  if (length(contig_ids) == 0L) return("v2")
  if (any(grepl("-cat_[0-6]$", contig_ids))) "v1" else "v2"
}


#' Read CAT/BAT results
#'
#' CAT/BAT is not part of the current ViroProfiler pipeline; this reader is kept
#' for datasets produced by earlier versions. The contig identifier is returned
#' as `Contig`, the name every other reader uses, so that the table can be joined
#' by [create_vpftse()].
#'
#' @param fin CAT/BAT taxonomy annotation file
#'
#' @return data.frame with `Contig` and `CATBAT_*` columns
#' @export
#'
read_catbat <- function(fin) {
  df <- .read_table(fin, "CAT/BAT", fill = TRUE, sep = "\t")
  .require_cols(df, "# contig", "CAT/BAT")
  ids <- .normalize_contig_ids(df[["# contig"]], what = "CAT/BAT")
  df <- df[, setdiff(colnames(df), "# contig"), drop = FALSE]
  colnames(df)[colnames(df) == "lineage scores"] <- "lineage_score"
  colnames(df) <- paste0("CATBAT_", colnames(df))
  df <- cbind(data.frame(Contig = ids, stringsAsFactors = FALSE), df)
  .dedup_by_contig(df, "CAT/BAT")
}


#' Read CheckV results
#'
#' Reads `quality_summary.tsv` as written by `checkv end_to_end`. Column names
#' are prefixed with `checkv_`, except `checkv_quality`, which CheckV already
#' names that way.
#'
#' CheckV marks an excised prophage with `provirus == "Yes"` and keeps the host
#' contig identifier in `contig_id`, so identifiers join directly; ViroProfiler's
#' `bin/run_checkv.sh` splits the two classes into separate files but does not
#' rename anything.
#'
#' @param fin CheckV `quality_summary.tsv`
#'
#' @return data.frame with a `Contig` column
#' @export
#'
read_checkv <- function(fin) {
  df <- .read_table(fin, "CheckV")
  .require_cols(df, c("contig_id", "checkv_quality", "completeness"), "CheckV")

  levels_quality <- c("Complete", "High-quality", "Medium-quality", "Low-quality",
                      "Not-determined")
  unknown <- setdiff(stats::na.omit(unique(df$checkv_quality)), levels_quality)
  if (length(unknown) > 0) {
    warning("CheckV: unrecognized checkv_quality value(s) become NA: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  ids <- .normalize_contig_ids(df$contig_id, what = "CheckV")
  df <- df[, setdiff(colnames(df), "contig_id"), drop = FALSE]
  df$checkv_quality <- factor(df$checkv_quality, levels = levels_quality)
  colnames(df) <- paste0("checkv_", colnames(df))
  colnames(df)[colnames(df) == "checkv_checkv_quality"] <- "checkv_quality"
  df <- cbind(data.frame(Contig = ids, stringsAsFactors = FALSE), df)
  .dedup_by_contig(df, "CheckV", order_by = "checkv_quality",
                   severity = "message")
}


#' Read a CoverM abundance table
#'
#' Reads one of the `abundance_contigs_*.tsv.gz` tables. The first column holds
#' contig identifiers; every other column is one sample.
#'
#' @param fpath CoverM output file
#' @param fbin2contig Optional two-column bin-to-contig mapping (no header). Rows
#'   are summed per bin, or maximized when `cov` is set. Pass `NULL` or `0` to
#'   skip.
#' @param cov Set to a non-zero value when the table is a coverage fraction, so
#'   that bins take the maximum rather than the sum of their contigs
#' @param sample_rename Optional named character vector of
#'   `pattern = replacement` applied to sample column names. Renaming samples is
#'   off by default: earlier versions rewrote any column matching `ds10Ms` to
#'   `Sample_`, which silently renamed samples in unrelated datasets.
#'
#' @return data.frame with contigs as rownames
#' @export
#'
read_coverm <- function(fpath, fbin2contig = 0, cov = 0, sample_rename = NULL) {
  df_abundance <- .read_table(fpath, "CoverM")
  if (!"Contig" %in% colnames(df_abundance)) {
    stop("CoverM file missing 'Contig' column. Found: ",
         paste(colnames(df_abundance), collapse = ", "), call. = FALSE)
  }
  if (ncol(df_abundance) < 2L) {
    stop("CoverM file has a 'Contig' column but no sample columns: ", fpath,
         call. = FALSE)
  }

  if (!is.null(sample_rename)) {
    nms <- colnames(df_abundance)
    for (pat in names(sample_rename)) {
      nms <- stringr::str_replace_all(nms, pat, sample_rename[[pat]])
    }
    colnames(df_abundance) <- nms
  }
  colnames(df_abundance)[colnames(df_abundance) == "Contig"] <- "genome_id"
  df_abundance$genome_id <- .normalize_contig_ids(df_abundance$genome_id,
                                                  what = "CoverM")

  sample_cols <- setdiff(colnames(df_abundance), "genome_id")
  non_numeric <- sample_cols[!vapply(df_abundance[sample_cols], is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop("CoverM file has non-numeric sample column(s): ",
         paste(non_numeric, collapse = ", "),
         ". as.matrix() would coerce the whole abundance matrix to character.",
         call. = FALSE)
  }

  if (!.is_absent_path(fbin2contig) && !identical(fbin2contig, 0)) {
    df_bin2contig <- data.table::fread(fbin2contig, header = FALSE, data.table = FALSE,
                                       col.names = c("bin_id", "genome_id"))
    df_abundance <- df_abundance %>%
      dplyr::left_join(df_bin2contig, by = "genome_id") %>%
      dplyr::mutate(bin_id = ifelse(is.na(.data$bin_id), .data$genome_id, .data$bin_id)) %>%
      dplyr::mutate(genome_id = .data$bin_id) %>%
      dplyr::select(-"bin_id") %>%
      dplyr::group_by(.data$genome_id)

    if (!identical(cov, 0) && !isFALSE(cov)) {
      df_abundance <- df_abundance %>% dplyr::summarise(dplyr::across(dplyr::everything(), max))
    } else {
      df_abundance <- df_abundance %>% dplyr::summarise(dplyr::across(dplyr::everything(), sum))
    }
    df_abundance <- as.data.frame(df_abundance, stringsAsFactors = FALSE)
  }

  if (anyDuplicated(df_abundance$genome_id)) {
    stop("CoverM file has duplicate contig identifiers: ",
         paste(utils::head(unique(df_abundance$genome_id[duplicated(df_abundance$genome_id)]), 5),
               collapse = ", "), call. = FALSE)
  }
  df_abundance %>% tibble::column_to_rownames("genome_id")
}


#' Adjust an abundance table by a coverage fraction table
#'
#' Zeroes every entry whose contig was covered over less than `covfrac_threshold`
#' of its length in that sample. Presence/absence in a virome is decided by
#' breadth of coverage rather than by depth, because a handful of reads piled on
#' one conserved region will otherwise be read as a present genome.
#'
#' @section Choice of threshold:
#' The default is 0.75, following Roux et al. (2017, *PeerJ* 5:e3817,
#' \doi{10.7717/peerj.3817}), whose simulations show that requiring at least 75%
#' of the contig covered at 1x, together with a strict identity cutoff,
#' suppresses false positives at a small cost in sensitivity. ViroProfiler runs
#' CoverM with `--min-read-percent-identity 0.95`, which is stricter than the
#' 90% used there. MetaPop applies 70% (Gregory et al. 2022,
#' \doi{10.1186/s40168-022-01231-0}). The ViroProfiler paper (Ru et al. 2023,
#' *Gut Microbes* 15:2192522, \doi{10.1080/19490976.2023.2192522}) used 0.5,
#' which was this function's previous default; raising it changes results, so
#' pass `covfrac_threshold = 0.5` to reproduce earlier analyses. CoverM's own
#' default minimum covered fraction is 0, a software default rather than a
#' biological recommendation.
#'
#' @section Order of operations:
#' Mask first, normalize second. Zeroing cells of an already-normalized matrix
#' leaves each column summing to something less than its original total, and by
#' a different amount in every sample, so the resulting compositions are not
#' comparable across samples. Apply this function to raw counts and derive
#' TPM or relative abundance afterwards, or pass `renormalize = TRUE` to rescale
#' each column back to its original total after masking.
#'
#' @param df_abundance An abundance table, ideally raw counts, e.g.
#'   `assay(tse, "counts")`
#' @param df_covfrac The matching coverage fraction table, e.g.
#'   `assay(tse, "covfrac")`
#' @param covfrac_threshold Minimum covered fraction to retain a value
#'   (default 0.75; see the section above)
#' @param renormalize Rescale every column to its pre-masking total. Use this
#'   only when `df_abundance` is already normalized (TPM, relative abundance) and
#'   raw counts are unavailable.
#'
#' @return matrix
#' @export
#'
abundance_adjust_by_covfrac <- function(df_abundance, df_covfrac,
                                        covfrac_threshold = 0.75,
                                        renormalize = FALSE) {
  if (covfrac_threshold <= 0) {
    warning("covfrac_threshold <= 0 keeps every value, including contigs with no ",
            "coverage at all.", call. = FALSE)
  }
  ab <- as.matrix(df_abundance)
  cf <- as.matrix(df_covfrac)
  if (!is.numeric(ab) || !is.numeric(cf)) {
    stop("Both tables must be numeric", call. = FALSE)
  }
  has_names <- !is.null(rownames(ab)) && !is.null(rownames(cf)) &&
    !is.null(colnames(ab)) && !is.null(colnames(cf))
  if (has_names) {
    if (!setequal(rownames(ab), rownames(cf)) || !setequal(colnames(ab), colnames(cf))) {
      stop("Abundance and coverage fraction tables describe different contigs or ",
           "samples, so multiplying them element by element would mask each contig ",
           "with another contig's coverage.", call. = FALSE)
    }
    cf <- cf[rownames(ab), colnames(ab), drop = FALSE]
  } else if (!identical(dim(ab), dim(cf))) {
    stop("Abundance and coverage fraction tables have different dimensions and no ",
         "dimnames to align them by.", call. = FALSE)
  } else {
    warning("One of the tables has no dimnames, so they are aligned by position. ",
            "Verify that both are in the same contig and sample order.", call. = FALSE)
  }

  # A missing covered fraction means the contig was not observed, so it is
  # masked. Comparing NA would otherwise propagate NA into the result.
  n_na <- sum(is.na(cf))
  if (n_na > 0) {
    message("Coverage fraction has ", n_na,
            " missing value(s); those entries are treated as not covered.")
  }
  mask <- !is.na(cf) & cf >= covfrac_threshold
  out <- ab * mask

  if (renormalize) {
    before <- colSums(ab, na.rm = TRUE)
    after <- colSums(out, na.rm = TRUE)
    scale <- ifelse(after > 0, before / after, 1)
    out <- sweep(out, 2, scale, `*`)
  }
  out
}


#' Convert reads per base to per-base depth
#'
#' CoverM's `reads_per_base` counts reads mapped per contig base; multiplying by
#' the read length gives bases per base, that is, fold coverage.
#'
#' The conversion assumes every mapped read contributes exactly `reads_len`
#' aligned bases, which it does not: fastp trims to variable lengths and the
#' aligner soft-clips. Treat the result as an approximation, and pass the mean
#' clean read length reported by fastp (`fastp/<sample>.fastp.json`, field
#' `summary$after_filtering$read1_mean_length`) rather than relying on the
#' default.
#'
#' @param fin Path of a CoverM `reads_per_base` table
#' @param reads_len Mean length of the clean reads
#'
#' @return matrix
#' @export
#'
rpb2bpb <- function(fin, reads_len = 150) {
  df_rpb <- .read_table(fin, "CoverM reads_per_base")
  if (!"Contig" %in% colnames(df_rpb)) {
    stop("CoverM reads_per_base file missing 'Contig' column", call. = FALSE)
  }
  df_rpb <- df_rpb %>% tibble::column_to_rownames("Contig")
  as.matrix(df_rpb * reads_len)
}


#' Read DeepVirFinder results
#'
#' DeepVirFinder is no longer run by ViroProfiler; the pipeline passes a
#' placeholder table so that the argument stays filled. The reader is kept for
#' datasets that do contain real DeepVirFinder output.
#'
#' Rows failing the thresholds are dropped, because [create_vpftse_vir()] reads
#' a non-missing `dvf_score` as a positive call. The thresholds applied are
#' recorded in `metadata(tse)$viral_selection`.
#'
#' @param fin DeepVirFinder result file
#' @param thr_score Minimum score
#' @param thr_pvalue Maximum p-value
#' @param thr_qvalue Maximum q-value
#' @param filter Keep only rows passing all three thresholds (default TRUE)
#'
#' @return data.frame with `Contig`, `dvf_score`, `dvf_pvalue`, `dvf_qvalue`
#' @export
#'
read_dvf <- function(fin, thr_score = 0.9, thr_pvalue = 0.01, thr_qvalue = 0.01,
                     filter = TRUE) {
  df <- .read_table(fin, "DVF")
  .require_cols(df, c("name", "score", "pvalue", "qvalue"), "DVF")
  out <- data.frame(
    Contig      = .normalize_contig_ids(df$name, what = "DVF"),
    dvf_score   = as.numeric(df$score),
    dvf_pvalue  = as.numeric(df$pvalue),
    dvf_qvalue  = as.numeric(df$qvalue),
    stringsAsFactors = FALSE
  )
  if (filter) {
    keep <- !is.na(out$dvf_score) & out$dvf_score > thr_score &
      !is.na(out$dvf_pvalue) & out$dvf_pvalue < thr_pvalue &
      !is.na(out$dvf_qvalue) & out$dvf_qvalue < thr_qvalue
    out <- out[keep, , drop = FALSE]
    rownames(out) <- NULL
  }
  .dedup_by_contig(out, "DVF", order_by = "dvf_score", decreasing = TRUE)
}


#' Read iPHoP host predictions
#'
#' Reads `Host_prediction_to_genus_mXX.csv`. The file is comma separated and its
#' header uses spaces (`AAI to closest reference`, `Host genus`,
#' `Confidence score`, `List of methods`); `data.table::fread()` detects the
#' separator, and the columns are matched by name with a positional fallback.
#'
#' iPHoP reports one row per candidate host genus, so a virus can appear several
#' times. The row with the highest confidence score is retained and the number of
#' candidates is kept in `iphop_n_predictions`.
#'
#' @param fin iPHoP genus prediction file
#'
#' @return data.frame with a `Contig` column
#' @export
#'
read_iphop <- function(fin) {
  df <- .read_table(fin, "iPHoP")
  if (ncol(df) < 5) stop("iPHoP file must have at least 5 columns", call. = FALSE)

  pick <- function(candidates, fallback_index) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) > 0) return(df[[hit[1]]])
    df[[fallback_index]]
  }
  out <- data.frame(
    Contig = .normalize_contig_ids(pick(c("Virus", "Contig"), 1L), what = "iPHoP"),
    iphop_aai2ref = suppressWarnings(as.numeric(
      pick(c("AAI to closest reference", "AAI to closest RaFAH reference",
             "iphop_aai2ref"), 2L))),
    iphop_genus = as.character(pick(c("Host genus", "iphop_genus"), 3L)),
    iphop_score = suppressWarnings(as.numeric(
      pick(c("Confidence score", "iphop_score"), 4L))),
    iphop_methods = as.character(pick(c("List of methods", "iphop_methods"), 5L)),
    stringsAsFactors = FALSE
  )
  counts <- table(out$Contig)
  out$iphop_n_predictions <- as.integer(counts[out$Contig])
  .dedup_by_contig(out, "iPHoP", order_by = "iphop_score", decreasing = TRUE,
                   severity = "message")
}


#' Read replication cycle predictions from BACPHLIP
#'
#' BACPHLIP writes a table whose first column header is empty, so
#' `data.table::fread()` names it `V1`.
#'
#' @param fin BACPHLIP result file
#' @param tool Tool name; only `"bacphlip"` is supported
#' @param version ViroProfiler version: `"auto"`, `"v1"` or `"v2"`
#'
#' @return data.frame with `Contig` and `bacphlip_replicyc`
#' @export
#'
read_replicyc <- function(fin, tool = "bacphlip", version = "auto") {
  if (!identical(tool, "bacphlip")) {
    stop("read_replicyc() supports tool = 'bacphlip' only; got '", tool, "'",
         call. = FALSE)
  }
  # BACPHLIP leaves the first header field empty, so fread has to be told that
  # the first line is a header.
  df <- .read_table(fin, "Replication cycle", header = TRUE)
  id_col <- if ("V1" %in% colnames(df)) "V1" else colnames(df)[1]
  .require_cols(df, c("Virulent", "Temperate"), "BACPHLIP")

  out <- data.frame(
    Contig = .normalize_contig_ids(df[[id_col]], what = "BACPHLIP"),
    bacphlip_replicyc = ifelse(df$Virulent > df$Temperate, "virulent", "temperate"),
    bacphlip_virulent_prob = as.numeric(df$Virulent),
    stringsAsFactors = FALSE
  )
  if (version == "auto") version <- .detect_vp_version(out$Contig)
  if (version == "v1") {
    out$Contig <- .normalize_contig_ids(out$Contig, cat_suffix = TRUE, what = "BACPHLIP")
  }
  .dedup_by_contig(out, "BACPHLIP")
}


#' Read an MMseqs2 taxonomy table
#'
#' MMseqs2 taxonomy was retired from ViroProfiler; taxonomy now comes from
#' VITAP and vConTACT3, merged by `bin/merge_taxonomy.py`, and is read by
#' [read_taxonomy2()]. This reader is kept for datasets produced earlier.
#'
#' @param fin Taxonomy annotation file
#' @param tool `"mmseqs"` or `"mmseq_ictv"`
#' @param version ViroProfiler version: `"auto"`, `"v1"` or `"v2"`
#'
#' @return data.frame with a `Contig` column
#' @export
#'
read_taxonomy <- function(fin, tool = "mmseqs", version = "auto") {
  if (tool == "mmseqs") {
    df <- .read_table(fin, "mmseqs taxonomy")
    .require_cols(df, c("genome_id", "superkingdom", "taxid"), "mmseqs taxonomy")
    df <- df[!is.na(df$taxid) & df$taxid != 0, , drop = FALSE]
    colnames(df)[colnames(df) == "genome_id"] <- "Contig"
    colnames(df)[colnames(df) == "superkingdom"] <- "Kingdom"
    colnames(df) <- stringr::str_to_title(colnames(df))
    df <- df[, setdiff(colnames(df), "Taxid"), drop = FALSE]
  } else if (tool == "mmseq_ictv") {
    df <- .read_table(fin, "mmseqs ICTV taxonomy", header = FALSE)
    df <- df[, c("V1", "V9"), drop = FALSE]
    colnames(df) <- c("Contig", "taxonomy")
  } else {
    stop("read_taxonomy() supports tool = 'mmseqs' or 'mmseq_ictv'; got '", tool, "'",
         call. = FALSE)
  }
  df <- .blank_to_na(df)
  if (version == "auto") version <- .detect_vp_version(df$Contig)
  if (version == "v1") {
    df$Contig <- .normalize_contig_ids(df$Contig, cat_suffix = TRUE,
                                       what = "mmseqs taxonomy")
  }
  .dedup_by_contig(df, "mmseqs taxonomy")
}


#' Read the merged ViroProfiler taxonomy table
#'
#' Reads `taxonomy_tse.tsv` written by `bin/merge_taxonomy.py`, which resolves
#' each rank independently from VITAP and vConTACT3. `Domain` holds the ICTV
#' realm, or the literal `Viruses` for a contig placed only at a lower rank;
#' `taxa_id` is 1 when a contig received any assignment and 0 when it received
#' none, and rows with 0 are dropped.
#'
#' Unresolved ranks are written as empty fields, not as `NA`, and are converted
#' to `NA` here so that an unassigned rank does not read as an assigned one.
#'
#' @param fin `taxonomy_tse.tsv`
#'
#' @return data.frame with `Contig` and one column per rank
#' @export
#'
read_taxonomy2 <- function(fin) {
  df <- .read_table(fin, "taxonomy")
  ranks <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus",
             "Species")
  .require_cols(df, c("contig_id", "taxa_id", ranks), "taxonomy")
  df <- df[!is.na(df$taxa_id) & df$taxa_id != 0, , drop = FALSE]
  out <- data.frame(
    Contig = .normalize_contig_ids(df$contig_id, what = "taxonomy"),
    stringsAsFactors = FALSE
  )
  for (r in ranks) out[[r]] <- as.character(df[[r]])
  out <- .blank_to_na(out, ranks)
  .dedup_by_contig(out, "taxonomy")
}


#' Read a tblastx table generated by Easyfig
#'
#' This is for adding links to gggenomes.
#'
#' @param fin `12.easyfig.out`
#' @param seq1 Query genome name
#' @param seq2 Target genome name
#' @param max_evalue E-value threshold
#' @param min_bitscore Bitscore threshold
#'
#' @return data.frame
#' @export
#'
read_tblastx <- function(fin, seq1 = 0, seq2 = 0, max_evalue = 0.001, min_bitscore = 50) {
  .check_input_file(fin, "tblastx")
  df <- data.table::fread(
    fin, header = FALSE, data.table = FALSE,
    col.names = c("seq_id", "seq_id2", "pident", "length", "mismatch", "gapopen",
                  "start", "end", "start2", "end2", "evalue", "bitscore"))

  if (!identical(seq1, 0)) df$seq_id <- seq1
  if (!identical(seq2, 0)) df$seq_id2 <- seq2
  df <- df[, c("seq_id", "start", "end", "seq_id2", "start2", "end2", "evalue",
               "bitscore", "pident"), drop = FALSE]
  df <- df[df$evalue < max_evalue & df$bitscore > min_bitscore, , drop = FALSE]
  rownames(df) <- NULL
  df
}


#' Import vConTACT2 results
#'
#' vConTACT2 is not part of the current ViroProfiler pipeline, which uses
#' vConTACT3 through `bin/merge_taxonomy.py`. This reader is kept for datasets
#' produced by earlier versions and has not been re-tested against real output.
#'
#' @param fin File `genome_by_genome_overview.csv` created by vConTACT2
#' @param assembler_label `"_NODE_"` for SPAdes
#' @param version ViroProfiler version: `"auto"`, `"v1"` or `"v2"`
#'
#' @return list of results
#' @export
#'
read_vcontact2 <- function(fin, assembler_label = "_NODE_", version = "auto") {
  . <- NULL
  seq_status <- NULL
  num_seqs <- NULL
  vConTACT_source <- NULL
  vConTACT_cluster_status <- NULL
  .check_input_file(fin, "vConTACT2")
  df_vcontact <- fread(fin) %>%
    mutate(source = ifelse(str_detect(.data$Genome, assembler_label), "queryseq", "refseq")) %>%
    mutate(cluster_status = ifelse(str_detect(.data[["VC Status"]], "Overlap"), "Overlap", .data[["VC Status"]])) %>%
    mutate(cluster_status = factor(.data$cluster_status)) %>%
    setnames(colnames(.), paste0("vConTACT_", str_replace_all(colnames(.), " ", "_")))

  vc2_refs <- df_vcontact %>%
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label, negate = TRUE)) %>%
    dplyr::filter(.data$vConTACT_VC != "")

  # Get VC stats
  vcontact_stats <- df_vcontact %>%
    dplyr::group_by(.data$vConTACT_source, .data$vConTACT_cluster_status) %>%
    dplyr::summarise(seqs_with_VC = sum(.data$vConTACT_VC != ""),
                     seqs_without_VC = sum(.data$vConTACT_VC == "")) %>%
    tidyr::gather(seq_status, num_seqs, -vConTACT_source, -vConTACT_cluster_status)

  plot_vc_stats <- ggplot(vcontact_stats, aes(y = .data$num_seqs, x = .data$seq_status, fill = .data$vConTACT_source, label = .data$num_seqs)) +
    geom_point(aes(color = .data$vConTACT_source), alpha = 0.5, size = 3) +
    facet_grid(.data$vConTACT_cluster_status ~ .) +
    geom_text_repel()

  # Annotate clusters using reference genomes
  vc2_contigs_vclst_anno <- df_vcontact %>%
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label)) %>%
    dplyr::filter(.data$vConTACT_VC != "") %>%
    dplyr::select("vConTACT_VC") %>%
    distinct() %>%
    inner_join(vc2_refs, by = "vConTACT_VC")

  # Whether contig clusters were annotated by reference (1) or not (0)
  df_vcontact <- df_vcontact %>%
    mutate(vConTACT_classified = ifelse(.data$vConTACT_VC %in% vc2_contigs_vclst_anno$vConTACT_VC, 1, 0)) %>%
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label)) %>%
    mutate(vConTACT_VC2 = ifelse(.data$vConTACT_VC_Status %in% c("Outlier", "Singleton"), .data$vConTACT_Genome, .data$vConTACT_VC)) %>%
    mutate(vConTACT_VC2 = ifelse(str_detect(.data$vConTACT_VC_Status, "Overlap"), .data$vConTACT_VC_Status, .data$vConTACT_VC2)) %>%
    setnames("vConTACT_Genome", "virsorter2_contig_id") %>%
    mutate(vConTACT_VC2 = str_replace_all(.data$vConTACT_VC2, "[\\s/]", "_")) %>%
    mutate(vConTACT_VC2 = str_replace_all(.data$vConTACT_VC2, "[\\(\\)]", ""))

  if (version == "auto") version <- .detect_vp_version(df_vcontact$virsorter2_contig_id)
  if (version == "v1") {
    df_vcontact <- df_vcontact %>%
      dplyr::mutate(Contig = str_replace(.data$virsorter2_contig_id, "-cat_[0-6]", "")) %>%
      dplyr::select(-"virsorter2_contig_id")
  }

  list("vc_tbl" = df_vcontact,
       "vc_stats" = vcontact_stats,
       "vc_plot" = plot_vc_stats,
       "vc_annotated" = vc2_contigs_vclst_anno,
       "vc_refs" = vc2_refs)
}


#' Read vConTACT2 results and return one data.frame
#'
#' @param fin File `genome_by_genome_overview.csv` created by vConTACT2
#' @param assembler_label `"_NODE_"` for SPAdes
#'
#' @return data.frame
#' @export
#'
read_vcontact2_simple <- function(fin, assembler_label = "_NODE_") {
  . <- NULL
  .check_input_file(fin, "vConTACT2")
  fread(fin) %>%
    dplyr::select(-"V1") %>%
    dplyr::filter(str_detect(.data$Genome, assembler_label)) %>%
    column_to_rownames("Genome") %>%
    setnames(colnames(.), paste0("vc_", str_replace_all(colnames(.), " ", "_"))) %>%
    rownames_to_column("contig_id")
}


#' Read VIBRANT genome quality results
#'
#' Reads `VIBRANT_genome_quality_<name>.tsv`, whose columns are `scaffold`,
#' `type` and `Quality`.
#'
#' VIBRANT names an excised prophage `<contig>_fragment_N`, so one input contig
#' can produce several rows. The suffix is stripped and the best quality is kept.
#' The quality levels are ordered best to worst — `complete circular`,
#' `high quality draft`, `medium quality draft`, `low quality draft` — which is
#' not their alphabetical order.
#'
#' @param fin VIBRANT genome quality file
#'
#' @return data.frame with `Contig`, `vibrant_replicyc` and `vibrant_quality`
#' @export
#'
read_vibrant <- function(fin) {
  df <- .read_table(fin, "VIBRANT")
  if (ncol(df) < 3) {
    stop("VIBRANT file must have at least 3 columns (scaffold, type, quality)",
         call. = FALSE)
  }
  pick <- function(candidates, fallback_index) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) > 0) return(df[[hit[1]]])
    df[[fallback_index]]
  }
  levels_quality <- c("complete circular", "high quality draft",
                      "medium quality draft", "low quality draft")
  quality <- as.character(pick(c("Quality", "quality"), 3L))
  unknown <- setdiff(stats::na.omit(unique(quality)), levels_quality)
  if (length(unknown) > 0) {
    warning("VIBRANT: unrecognized quality value(s) become NA: ",
            paste(unknown, collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    Contig = .normalize_contig_ids(pick(c("scaffold", "Scaffold"), 1L),
                                   vibrant_fragment = TRUE, what = "VIBRANT"),
    vibrant_replicyc = as.character(pick(c("type", "Type"), 2L)),
    vibrant_quality = factor(quality, levels = levels_quality),
    stringsAsFactors = FALSE
  )
  # Ordering by the factor keeps the best quality; the previous implementation
  # sorted by the string literal "vibrant_quality", which sorts nothing.
  .dedup_by_contig(out, "VIBRANT", order_by = "vibrant_quality",
                   severity = "message")
}


#' Read VirSorter2 results
#'
#' Reads `final-viral-score.tsv`. Column names are prefixed with `virsorter2_`.
#'
#' The per-group score columns vary with `--include-groups`, so only `seqname`,
#' `max_score` and `max_score_group` are required. ViroProfiler runs VirSorter2
#' with `--seqname-suffix-off`, so `seqname` carries no `||full` decoration; the
#' decoration is stripped anyway for output produced without that flag.
#'
#' In ViroProfiler, VirSorter2 runs on contigs that were already selected as
#' putative viruses, in order to produce the DRAM-v input. It is not an
#' independent detector.
#'
#' @param fin VirSorter2 `final-viral-score.tsv`
#'
#' @return data.frame with a `Contig` column
#' @export
#'
read_virsorter2 <- function(fin) {
  df <- .read_table(fin, "VirSorter2")
  .require_cols(df, c("seqname", "max_score", "max_score_group"), "VirSorter2")
  ids <- .normalize_contig_ids(df$seqname, vs2_suffix = TRUE, what = "VirSorter2")
  df <- df[, setdiff(colnames(df), "seqname"), drop = FALSE]
  colnames(df) <- paste0("virsorter2_", colnames(df))
  df <- cbind(data.frame(Contig = ids, stringsAsFactors = FALSE), df)
  .dedup_by_contig(df, "VirSorter2", order_by = "virsorter2_max_score",
                   decreasing = TRUE)
}


#' Read DRAM-v annotations
#'
#' Reads `dramv-annotate/annotations.tsv`. DRAM-v writes the gene identifier in
#' an unnamed first column, which `data.table::fread()` names `V1`; the separate
#' `fasta` column holds the stem of the input FASTA name — the same value on
#' every row, normally `final-viral-combined-for-dramv` — and `scaffold` holds
#' the sequence the gene was called on. Earlier versions of this reader used
#' `fasta` as the gene identifier, which gave every gene the same ID.
#'
#' The column block written for each database is conditional on how DRAM was
#' configured, so columns are located by name and the optional ones are filled
#' with `NA` when absent.
#'
#' DRAM-v is run on the VirSorter2 `--prep-for-dramv` FASTA, whose names carry a
#' `-cat_N` suffix and, unless VirSorter2 ran with `--seqname-suffix-off`, a
#' `__full` / `__<i>_partial` / `__lt2gene` token as well. Both are stripped so
#' that `Contig` joins to the assembled contig.
#'
#' @param fin DRAM-v `annotations.tsv`
#' @return data.frame with gene-level annotations
#' @export
read_dramv <- function(fin) {
  # DRAM-v leaves the first header field empty, so fread must be told there is
  # a header; otherwise it can decide the header line is data.
  df <- .read_table(fin, "DRAM-v", header = TRUE)
  .require_cols(df, c("scaffold", "rank"), "DRAM-v")
  gene_col <- if ("V1" %in% colnames(df)) "V1" else
    if ("gene" %in% colnames(df)) "gene" else
      if (identical(colnames(df)[1], "fasta")) "fasta" else colnames(df)[1]
  if (identical(gene_col, "fasta")) {
    warning("DRAM-v: no unnamed gene index column was found, so the 'fasta' column ",
            "is used as the gene identifier. In real DRAM-v output that column is ",
            "the input FASTA stem and is identical on every row.", call. = FALSE)
  }
  out <- data.frame(
    Contig = .normalize_contig_ids(df$scaffold, vs2_suffix = TRUE, cat_suffix = TRUE,
                                   dramv_suffix = TRUE, what = "DRAM-v"),
    gene_id = as.character(df[[gene_col]]),
    dramv_category = dplyr::case_when(
      df$rank == "A" ~ "AMG",
      df$rank == "V" ~ "viral",
      df$rank == "H" ~ "host",
      TRUE ~ "other"
    ),
    dramv_ko = as.character(.col_or_na(df, "ko_id")),
    dramv_pfam = as.character(.col_or_na(df, "pfam_hits")),
    dramv_vog = as.character(.col_or_na(df, "vogdb")),
    dramv_amg_flags = as.character(.col_or_na(df, "amg_flags")),
    stringsAsFactors = FALSE
  )
  if ("auxiliary_score" %in% colnames(df)) {
    out$auxiliary_score <- suppressWarnings(as.numeric(df$auxiliary_score))
  }
  out
}


#' Summarize DRAM-v annotations per contig
#'
#' @param df_dramv Output from [read_dramv()]
#' @return data.frame with per-contig summary columns
#' @noRd
.summarize_dramv <- function(df_dramv) {
  df_dramv %>%
    dplyr::group_by(.data$Contig) %>%
    dplyr::summarise(
      dramv_amg_count = sum(.data$dramv_category == "AMG", na.rm = TRUE),
      dramv_viral_gene_count = sum(.data$dramv_category == "viral", na.rm = TRUE),
      dramv_total_gene_count = dplyr::n(),
      .groups = "drop"
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
}


#' Read pharokka functional annotations
#'
#' pharokka is not part of the current ViroProfiler pipeline; this reader is kept
#' for externally produced annotations and has not been re-tested against real
#' output.
#'
#' @param fin pharokka CDS output TSV file
#' @return data.frame with gene-level annotations
#' @export
read_pharokka <- function(fin) {
  df <- .read_table(fin, "pharokka")
  .require_cols(df, c("gene", "contig"), "pharokka")
  data.frame(
    Contig = .normalize_contig_ids(df$contig, what = "pharokka"),
    gene_id = as.character(df$gene),
    pharokka_function = as.character(.col_or_na(df, "function")),
    pharokka_phrog = as.character(.col_or_na(df, "phrog")),
    pharokka_category = as.character(.col_or_na(df, "phrog_category")),
    pharokka_card = as.character(.col_or_na(df, "card_match")),
    pharokka_vfdb = as.character(.col_or_na(df, "vfdb_match")),
    stringsAsFactors = FALSE
  )
}


#' Summarize pharokka annotations per contig
#'
#' @param df_pharokka Output from [read_pharokka()]
#' @return data.frame with per-contig summary columns
#' @noRd
.summarize_pharokka <- function(df_pharokka) {
  df_pharokka %>%
    dplyr::group_by(.data$Contig) %>%
    dplyr::summarise(
      pharokka_gene_count = dplyr::n(),
      pharokka_card_count = sum(!is.na(.data$pharokka_card), na.rm = TRUE),
      pharokka_vfdb_count = sum(!is.na(.data$pharokka_vfdb), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
}


#' Read CheckAMG results
#'
#' Reads `final_results.tsv` from a `checkamg annotate` run, the table
#' ViroProfiler's `CHECKAMG` process checks for content before publishing. Pass
#' either the file or the `checkamg_results` directory that the process emits.
#'
#' CheckAMG classifies each protein as `metabolic` (an auxiliary metabolic gene,
#' AMG), `physiological` (AVG/APG), `regulatory` (AReG) or `unclassified`, and
#' grades the confidence that the protein really is of viral origin as `high`,
#' `medium` or `low`. Its column names contain spaces.
#'
#' @param fin `final_results.tsv`, or the `checkamg_results` directory
#' @return data.frame with gene-level annotations
#' @export
read_checkamg <- function(fin) {
  if (!.is_absent_path(fin) && dir.exists(fin)) {
    candidate <- file.path(fin, "final_results.tsv")
    if (!file.exists(candidate)) {
      candidate <- file.path(fin, "checkamg_results", "final_results.tsv")
    }
    if (!file.exists(candidate)) {
      stop("CheckAMG directory has no final_results.tsv: ", fin, call. = FALSE)
    }
    fin <- candidate
  }
  df <- .read_table(fin, "CheckAMG", check.names = FALSE)
  .require_cols(df, c("Protein", "Contig"), "CheckAMG")

  cls <- as.character(.col_or_na(df, "Protein Classification"))
  data.frame(
    Contig = .normalize_contig_ids(df$Contig, what = "CheckAMG"),
    gene_id = as.character(df$Protein),
    checkamg_class = cls,
    checkamg_confidence = as.character(.col_or_na(df, "Protein Viral Origin Confidence")),
    checkamg_in_viral_region = as.character(.col_or_na(df, "Protein in Strict Viral Region")),
    checkamg_function = as.character(.col_or_na(df, "Function")),
    checkamg_kegg_ko = as.character(.col_or_na(df, "KEGG KO")),
    checkamg_pfam = as.character(.col_or_na(df, "Pfam Accession")),
    checkamg_cazy = as.character(.col_or_na(df, "CAZy Family")),
    checkamg_phrog = as.character(.col_or_na(df, "PHROG Number")),
    stringsAsFactors = FALSE
  )
}


#' Summarize CheckAMG annotations per contig
#'
#' @param df_checkamg Output from [read_checkamg()]
#' @return data.frame with per-contig summary columns
#' @noRd
.summarize_checkamg <- function(df_checkamg) {
  df_checkamg %>%
    dplyr::group_by(.data$Contig) %>%
    dplyr::summarise(
      checkamg_gene_count = dplyr::n(),
      checkamg_amg_count = sum(.data$checkamg_class == "metabolic", na.rm = TRUE),
      checkamg_apg_count = sum(.data$checkamg_class == "physiological", na.rm = TRUE),
      checkamg_areg_count = sum(.data$checkamg_class == "regulatory", na.rm = TRUE),
      checkamg_high_confidence_count = sum(.data$checkamg_confidence == "high", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    as.data.frame(stringsAsFactors = FALSE)
}


#' Read a geNomad virus summary
#'
#' Reads `*_virus_summary.tsv`. `fdr` is only computed when geNomad runs with
#' score calibration enabled and is otherwise written as `NA` in every row, which
#' makes the column read as logical; it is coerced to numeric so that the rowData
#' schema does not change between runs.
#'
#' geNomad names an excised provirus `<contig>|provirus_N`; the decoration is
#' stripped so that the row joins to its host contig, and the untouched name is
#' kept in `genomad_seq_name`.
#'
#' @param fin geNomad virus summary TSV file
#' @return data.frame with a `Contig` column
#' @export
read_genomad <- function(fin) {
  df <- .read_table(fin, "geNomad")
  .require_cols(df, c("seq_name", "virus_score"), "geNomad")
  out <- data.frame(
    Contig = .normalize_contig_ids(df$seq_name, genomad_provirus = TRUE,
                                   what = "geNomad"),
    genomad_seq_name = as.character(df$seq_name),
    genomad_score = suppressWarnings(as.numeric(df$virus_score)),
    genomad_fdr = suppressWarnings(as.numeric(.col_or_na(df, "fdr"))),
    genomad_topology = as.character(.col_or_na(df, "topology")),
    genomad_taxonomy = as.character(.col_or_na(df, "taxonomy")),
    genomad_n_hallmarks = suppressWarnings(as.numeric(.col_or_na(df, "n_hallmarks"))),
    stringsAsFactors = FALSE
  )
  .dedup_by_contig(out, "geNomad", order_by = "genomad_score", decreasing = TRUE,
                   severity = "message")
}


#' Read vRhyme viral bin assignments
#'
#' Reads `out_vrhyme/vRhyme_best_bins.<N>.membership.tsv`, whose header is
#' `scaffold`, `bin`, in that order. A generic `contig` column is also accepted.
#'
#' `bin` holds a plain integer, not a `vRhyme_bin_1` style label; the
#' `vRhyme_<n>__` prefix appears only in the per-bin FASTA files. Unbinned
#' scaffolds are absent from the file, so they end up with `NA` after the join.
#'
#' @param fin vRhyme bin membership file
#' @return data.frame with `Contig` and `vrhyme_bin`
#' @export
read_vrhyme <- function(fin) {
  df <- .read_table(fin, "vRhyme")
  id_col <- intersect(c("scaffold", "contig", "Contig", "sequence"), colnames(df))
  bin_col <- intersect(c("bin", "Bin", "membership"), colnames(df))
  if (length(id_col) == 0L || length(bin_col) == 0L) {
    stop("vRhyme file missing columns: expected a scaffold/contig column and a bin ",
         "column. Found: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    Contig = .normalize_contig_ids(df[[id_col[1]]], what = "vRhyme"),
    vrhyme_bin = as.character(df[[bin_col[1]]]),
    stringsAsFactors = FALSE
  )
  .dedup_by_contig(out, "vRhyme")
}


#' Read PHIST host predictions
#'
#' PHIST is not part of the current ViroProfiler pipeline; this reader is kept
#' for externally produced predictions and has not been re-tested against real
#' output. Both the capitalized (`Virus`/`Host`/`Score`) and the lowercase
#' (`phage`/`host`) header spellings are accepted.
#'
#' @param fin PHIST predictions file
#' @return data.frame with a `Contig` column
#' @export
read_phist <- function(fin) {
  df <- .read_table(fin, "PHIST")
  pick_name <- function(candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) > 0) hit[1] else NA_character_
  }
  virus_col <- pick_name(c("Virus", "virus", "phage", "Phage"))
  host_col <- pick_name(c("Host", "host"))
  score_col <- pick_name(c("Score", "score", "#common kmers", "adj-p-value"))
  if (anyNA(c(virus_col, host_col, score_col))) {
    stop("PHIST file missing columns: expected virus, host and score columns. ",
         "Found: ", paste(colnames(df), collapse = ", "), call. = FALSE)
  }
  out <- data.frame(
    Contig = .normalize_contig_ids(df[[virus_col]], what = "PHIST"),
    phist_host = as.character(df[[host_col]]),
    phist_score = suppressWarnings(as.numeric(df[[score_col]])),
    phist_host_taxonomy = as.character(.col_or_na(df, "Host_taxonomy")),
    stringsAsFactors = FALSE
  )
  .dedup_by_contig(out, "PHIST", order_by = "phist_score", decreasing = TRUE)
}


#' Read a list of contig identifiers, one per line
#'
#' @param fin A `.list` file, e.g. `vircontigs/putative_vcontigs_pref1.list`
#' @return Character vector
#' @noRd
.read_contig_list <- function(fin) {
  .check_input_file(fin, "Contig list")
  ids <- trimws(readLines(fin, warn = FALSE))
  unique(ids[nzchar(ids)])
}


#' Merge a per-contig summary into a rownames-keyed feature annotation
#'
#' @param feature_anno data.frame with Contig rownames
#' @param summary_df data.frame with a Contig column to left-join
#' @param what Tool name used in messages
#' @return data.frame with Contig rownames
#' @noRd
.merge_summary_into_rowdata <- function(feature_anno, summary_df, what = "summary") {
  keep <- rownames(feature_anno)
  out <- feature_anno %>%
    tibble::rownames_to_column("Contig") %>%
    .join_tool_table(summary_df, what) %>%
    tibble::column_to_rownames("Contig")
  out[keep, , drop = FALSE]
}


####################################################################################
## TSE construction
####################################################################################


#' Description of each assay written by create_vpftse()
#'
#' Kept next to the constructor because the reason `tmm` was wrong is that no
#' single place said what any of the four matrices measured.
#'
#' @noRd
.VPF_ASSAY_DOC <- list(
  counts = paste(
    "Number of reads mapped to the contig (CoverM `--methods count`).",
    "Unit: reads. Not normalized for contig length or library size."),
  tpm = paste(
    "Transcripts per million (CoverM `--methods tpm`).",
    "Unit: parts per million; columns sum to 1e6. Normalized for both contig",
    "length and library size."),
  trimmed_mean = paste(
    "Trimmed mean of per-base coverage depth (CoverM `--methods trimmed_mean`),",
    "computed after discarding the highest and lowest coverage positions.",
    "Unit: fold coverage. This is NOT edgeR's TMM (trimmed mean of M-values)",
    "normalization; earlier versions of vpfkit stored it under the name `tmm`."),
  covfrac = paste(
    "Fraction of the contig covered by at least one read",
    "(CoverM `--methods covered_fraction`). Unit: proportion in [0, 1].")
)


#' Create a ViroProfiler TSE object
#'
#' Assembles the abundance matrices and every per-contig annotation table into a
#' `TreeSummarizedExperiment` over all assembled contigs.
#'
#' @section Assays:
#' Four assays are written, all produced by CoverM from the same BAM files:
#' \describe{
#'   \item{`counts`}{Reads mapped per contig (`--methods count`).}
#'   \item{`tpm`}{Transcripts per million (`--methods tpm`).}
#'   \item{`trimmed_mean`}{Trimmed mean of per-base coverage depth
#'     (`--methods trimmed_mean`), in fold coverage. This is not edgeR's TMM
#'     normalization; objects written before vpfkit 0.6 stored it as `tmm`, and
#'     [normalize_assay_names()] renames it when such an object is read.}
#'   \item{`covfrac`}{Covered fraction of the contig (`--methods
#'     covered_fraction`), in [0, 1].}
#' }
#'
#' @section Metadata:
#' `metadata(tse)$gene_annotations` holds the gene-level tables from DRAM-v,
#' pharokka and CheckAMG stacked into one data.frame with a `source` column;
#' `metadata(tse)$gene_annotations_by_tool` keeps them separate.
#' `metadata(tse)$viroprofiler` records the vpfkit version, the time the object
#' was built, the source file of every input, and the per-tool join match rates.
#'
#' @param fin_abcount File `abundance_contigs_count.tsv.gz` created by CoverM
#' @param fin_abtpm File `abundance_contigs_tpm.tsv.gz` created by CoverM
#' @param fin_abtmm File `abundance_contigs_trimmed_mean.tsv.gz` created by
#'   CoverM. The argument keeps its historical name; the assay it fills is called
#'   `trimmed_mean`.
#' @param fin_abcov File `abundance_contigs_covered_fraction.tsv.gz` created by
#'   CoverM
#' @param fin_taxa File `taxonomy_tse.tsv` created by `bin/merge_taxonomy.py`
#' @param fin_checkv File `quality_summary.tsv` created by CheckV
#' @param fin_virsorter2 File `final-viral-score.tsv` created by VirSorter2
#' @param fin_vibrant File `VIBRANT_genome_quality_contigs.tsv` created by
#'   VIBRANT (optional; VIBRANT is switchable with `--use_vibrant`)
#' @param fin_dvf File `dvf_virus.tsv` created by DeepVirFinder (optional;
#'   DeepVirFinder is no longer part of the pipeline)
#' @param fin_replicyc File `*.bacphlip` created by BACPHLIP
#' @param df_metadata Sample metadata as a data.frame (optional)
#' @param fin_genomad geNomad virus summary TSV (optional)
#' @param fin_vrhyme vRhyme bin membership file (optional)
#' @param fin_phist PHIST host predictions (optional)
#' @param fin_dramv DRAM-v `annotations.tsv` (optional)
#' @param fin_pharokka pharokka CDS output TSV (optional)
#' @param fin_checkamg CheckAMG `final_results.tsv`, or the `checkamg_results`
#'   directory (optional)
#' @param fin_iphop iPHoP `Host_prediction_to_genus_mXX.csv` (optional)
#' @param fin_metadata Sample metadata file, CSV/TSV/XLSX (optional). Samples are
#'   matched by name against the abundance table columns; unmatched names are
#'   reported rather than dropped silently. Mutually exclusive with
#'   `df_metadata`.
#' @param fin_coverm_log A CoverM log such as `abundance/log_contig_count.txt`
#'   (optional). Adds `n_reads_total`, `n_reads_mapped` and `mapping_rate` to
#'   colData.
#' @param fin_vircontigs The putative viral contig list written by
#'   `VIRCONTIGS_PRE`, e.g. `vircontigs/putative_vcontigs_pref1.list` (optional).
#'   Adds the logical rowData column `upstream_viral_candidate`, which records
#'   exactly which contigs the pipeline itself selected as putative viruses.
#' @param pipeline_info Named list of pipeline provenance (versions, parameters)
#'   stored in `metadata(tse)$viroprofiler$pipeline` (optional)
#'
#' @return TreeSummarizedExperiment object
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
create_vpftse <- function(fin_abcount, fin_abtpm, fin_abtmm, fin_abcov, fin_taxa,
                          fin_checkv, fin_virsorter2, fin_vibrant = NULL, fin_dvf = NULL,
                          fin_replicyc, df_metadata = NULL, fin_genomad = NULL,
                          fin_vrhyme = NULL, fin_phist = NULL, fin_dramv = NULL,
                          fin_pharokka = NULL, fin_checkamg = NULL, fin_iphop = NULL,
                          fin_metadata = NULL, fin_coverm_log = NULL,
                          fin_vircontigs = NULL, pipeline_info = NULL) {

  assay_list <- list(
    counts       = as.matrix(read_coverm(fin_abcount)),
    tpm          = as.matrix(read_coverm(fin_abtpm)),
    trimmed_mean = as.matrix(read_coverm(fin_abtmm)),
    covfrac      = as.matrix(read_coverm(fin_abcov))
  )
  assay_list <- .align_assays(assay_list)
  contigs <- rownames(assay_list$counts)
  samples <- colnames(assay_list$counts)

  feature_anno <- data.frame(Contig = contigs, stringsAsFactors = FALSE,
                             check.names = FALSE)
  match_stats <- list()
  join_tool <- function(anno, df, what) {
    out <- .join_tool_table(anno, df, what)
    st <- attr(out, "vpf_match")
    if (!is.null(st)) match_stats[[what]] <<- st
    attr(out, "vpf_match") <- NULL
    out
  }

  feature_anno <- join_tool(feature_anno, read_taxonomy2(fin_taxa), "taxonomy")
  feature_anno <- join_tool(feature_anno, read_checkv(fin_checkv), "CheckV")
  feature_anno <- join_tool(feature_anno, read_virsorter2(fin_virsorter2), "VirSorter2")
  feature_anno <- join_tool(feature_anno, read_replicyc(fin_replicyc), "BACPHLIP")

  # VIBRANT is switchable with --use_vibrant and DeepVirFinder was removed from
  # the pipeline altogether, so neither is required. A vote whose column is
  # absent is reported by annotate_viral_votes() rather than silently counted as
  # a negative result.
  fin_vibrant <- .optional_path(fin_vibrant, "VIBRANT")
  fin_dvf <- .optional_path(fin_dvf, "DVF")
  if (!is.null(fin_vibrant)) {
    feature_anno <- join_tool(feature_anno, read_vibrant(fin_vibrant), "VIBRANT")
  }
  if (!is.null(fin_dvf)) {
    feature_anno <- join_tool(feature_anno, read_dvf(fin_dvf, thr_qvalue = 0.1), "DVF")
  }

  # Optional per-contig annotations. Every one of these arrives as NULL, NA, ""
  # or the literal "null" depending on how the caller was invoked.
  fin_genomad    <- .optional_path(fin_genomad, "geNomad")
  fin_vrhyme     <- .optional_path(fin_vrhyme, "vRhyme")
  fin_phist      <- .optional_path(fin_phist, "PHIST")
  fin_dramv      <- .optional_path(fin_dramv, "DRAM-v")
  fin_pharokka   <- .optional_path(fin_pharokka, "pharokka")
  fin_checkamg   <- .optional_path(fin_checkamg, "CheckAMG")
  fin_iphop      <- .optional_path(fin_iphop, "iPHoP")
  fin_metadata   <- .optional_path(fin_metadata, "Sample metadata")
  fin_coverm_log <- .optional_path(fin_coverm_log, "CoverM log")
  fin_vircontigs <- .optional_path(fin_vircontigs, "Putative viral contig list")

  if (!is.null(fin_genomad)) {
    feature_anno <- join_tool(feature_anno, read_genomad(fin_genomad), "geNomad")
  }
  if (!is.null(fin_vrhyme)) {
    feature_anno <- join_tool(feature_anno, read_vrhyme(fin_vrhyme), "vRhyme")
  }
  if (!is.null(fin_phist)) {
    feature_anno <- join_tool(feature_anno, read_phist(fin_phist), "PHIST")
  }
  if (!is.null(fin_iphop)) {
    feature_anno <- join_tool(feature_anno, read_iphop(fin_iphop), "iPHoP")
  }
  if (!is.null(fin_vircontigs)) {
    candidates <- .read_contig_list(fin_vircontigs)
    feature_anno$upstream_viral_candidate <- contigs %in% candidates
    match_stats[["putative viral contigs"]] <- c(
      n_tool = length(candidates),
      n_matched = sum(feature_anno$upstream_viral_candidate),
      n_features = length(contigs))
  }

  feature_anno <- feature_anno %>% tibble::column_to_rownames("Contig")

  # Gene-level annotations live in metadata(); only their per-contig summaries
  # belong in rowData, which has one row per contig.
  gene_tables <- list()
  if (!is.null(fin_dramv)) {
    gene_tables$dramv <- read_dramv(fin_dramv)
    feature_anno <- .merge_summary_into_rowdata(
      feature_anno, .summarize_dramv(gene_tables$dramv), "DRAM-v summary")
  }
  if (!is.null(fin_pharokka)) {
    gene_tables$pharokka <- read_pharokka(fin_pharokka)
    feature_anno <- .merge_summary_into_rowdata(
      feature_anno, .summarize_pharokka(gene_tables$pharokka), "pharokka summary")
  }
  if (!is.null(fin_checkamg)) {
    gene_tables$checkamg <- read_checkamg(fin_checkamg)
    feature_anno <- .merge_summary_into_rowdata(
      feature_anno, .summarize_checkamg(gene_tables$checkamg), "CheckAMG summary")
  }

  # Sample metadata
  if (!is.null(fin_metadata) && !is.null(df_metadata)) {
    stop("Pass either df_metadata or fin_metadata, not both", call. = FALSE)
  }
  if (!is.null(fin_metadata)) df_metadata <- .read_metadata_table(fin_metadata)
  col_data <- .build_coldata(samples, df_metadata)
  if (!is.null(fin_coverm_log)) {
    depth <- read_coverm_log(fin_coverm_log)
    if (nrow(depth) > 0) {
      idx <- match(samples, depth$sample_id)
      if (all(is.na(idx))) {
        warning("CoverM log sample names match no assay column; sequencing depth ",
                "was not added. Log names: ",
                paste(depth$sample_id, collapse = ", "), call. = FALSE)
      } else {
        col_data$n_reads_total <- depth$n_reads_total[idx]
        col_data$n_reads_mapped <- depth$n_reads_mapped[idx]
        col_data$mapping_rate <- depth$mapping_rate[idx]
      }
    }
  }

  stopifnot(identical(rownames(feature_anno), contigs),
            identical(rownames(col_data), samples))

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = S4Vectors::SimpleList(assay_list),
    colData = MultiAssayExperiment::DataFrame(col_data),
    rowData = feature_anno)

  if (length(gene_tables) > 0) {
    combined <- dplyr::bind_rows(lapply(names(gene_tables), function(nm) {
      x <- gene_tables[[nm]]
      x$source <- nm
      x
    }))
    S4Vectors::metadata(tse)$gene_annotations <- combined
    S4Vectors::metadata(tse)$gene_annotations_by_tool <- gene_tables
  }

  S4Vectors::metadata(tse)$viroprofiler <- list(
    vpfkit_version = .vpfkit_version(),
    created = Sys.time(),
    assays = .VPF_ASSAY_DOC,
    source_files = Filter(Negate(is.null), list(
      abundance_count = fin_abcount, abundance_tpm = fin_abtpm,
      abundance_trimmed_mean = fin_abtmm, abundance_covfrac = fin_abcov,
      taxonomy = fin_taxa, checkv = fin_checkv, virsorter2 = fin_virsorter2,
      vibrant = fin_vibrant, dvf = fin_dvf, replicyc = fin_replicyc,
      genomad = fin_genomad, vrhyme = fin_vrhyme, phist = fin_phist,
      dramv = fin_dramv, pharokka = fin_pharokka, checkamg = fin_checkamg,
      iphop = fin_iphop, sample_metadata = fin_metadata,
      coverm_log = fin_coverm_log, vircontigs = fin_vircontigs)),
    join_match = match_stats,
    pipeline = pipeline_info
  )

  tse
}


#' Installed version of vpfkit, if it can be determined
#'
#' @return character(1)
#' @noRd
.vpfkit_version <- function() {
  tryCatch(as.character(getNamespaceVersion("vpfkit")),
           error = function(e) NA_character_)
}


#' Create a TSE object from a ViroProfiler output directory
#'
#' Auto-discovers ViroProfiler output files from a standard directory structure
#' and assembles a `TreeSummarizedExperiment` without manual file path
#' specification.
#'
#' @param vpdir Path to a ViroProfiler output directory
#' @param df_metadata Optional sample metadata data.frame
#' @param version Kept for backward compatibility and not used: each reader
#'   detects the identifier convention from the identifiers themselves
#' @return TreeSummarizedExperiment object
#' @export
batch_create_vpftse <- function(vpdir, df_metadata = NULL, version = "auto") {
  if (!dir.exists(vpdir)) stop("Directory not found: ", vpdir, call. = FALSE)
  if (!identical(version, "auto")) {
    warning("batch_create_vpftse() ignores `version`; each reader detects the ",
            "identifier convention itself. Passing '", version,
            "' has no effect.", call. = FALSE)
  }

  .find_file <- function(pattern, required = FALSE) {
    matches <- list.files(vpdir, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(matches) == 0) {
      if (required) {
        stop("Required file not found matching pattern '", pattern, "' in ", vpdir,
             call. = FALSE)
      }
      return(NULL)
    }
    if (length(matches) > 1) {
      message("Several files match '", pattern, "'; using ", matches[1])
    }
    matches[1]
  }

  fin_abcount <- .find_file("abundance_contigs_count", required = TRUE)
  fin_abtpm <- .find_file("abundance_contigs_tpm", required = TRUE)
  fin_abtmm <- .find_file("abundance_contigs_(tmm|trimmed_mean)", required = TRUE)
  fin_abcov <- .find_file("abundance_contigs_(covfrac|covered_fraction)", required = TRUE)
  fin_taxa <- .find_file("(taxonomy_tse|taxa_mmseqs_formatted).*\\.tsv", required = TRUE)
  fin_checkv <- .find_file("quality_summary\\.tsv", required = TRUE)
  fin_virsorter2 <- .find_file("final-viral-score", required = TRUE)
  fin_vibrant <- .find_file("VIBRANT_genome_quality", required = TRUE)
  fin_dvf <- .find_file("dvf_virus", required = TRUE)
  fin_replicyc <- .find_file("\\.bacphlip$", required = TRUE)

  create_vpftse(
    fin_abcount = fin_abcount,
    fin_abtpm = fin_abtpm,
    fin_abtmm = fin_abtmm,
    fin_abcov = fin_abcov,
    fin_taxa = fin_taxa,
    fin_checkv = fin_checkv,
    fin_virsorter2 = fin_virsorter2,
    fin_vibrant = fin_vibrant,
    fin_dvf = fin_dvf,
    fin_replicyc = fin_replicyc,
    df_metadata = df_metadata,
    fin_genomad = .find_file("virus_genomad_summary\\.tsv|virus_summary\\.tsv"),
    fin_vrhyme = .find_file("membership\\.tsv|bin_to_contig"),
    fin_dramv = .find_file("annotations\\.tsv"),
    fin_pharokka = .find_file("cds_final_merged_output"),
    fin_checkamg = .find_file("final_results\\.tsv"),
    fin_iphop = .find_file("Host_prediction_to_genus"),
    fin_coverm_log = .find_file("log_contig_count\\.txt"),
    fin_vircontigs = .find_file("putative_vcontigs.*\\.list")
  )
}


#' Definition of each viral-identity vote
#'
#' `column` is the rowData column the vote reads and `kind` is how it is read.
#'
#' @noRd
.VPF_VOTES <- list(
  taxonomy   = list(column = "Domain",                     kind = "not_na"),
  checkv     = list(column = "checkv_quality",             kind = "in_levels"),
  virsorter2 = list(column = "virsorter2_max_score_group", kind = "in_levels"),
  vibrant    = list(column = "vibrant_quality",            kind = "not_na"),
  dvf        = list(column = "dvf_score",                  kind = "not_na"),
  genomad    = list(column = "genomad_score",              kind = "min_score")
)


#' Record which evidence marks each contig as viral
#'
#' Adds one logical rowData column per vote (`viral_vote_<name>`), the number of
#' votes in favor (`viral_vote_n`), a semicolon-separated list of the votes that
#' fired (`viral_vote_evidence`), and the decision itself (`viral_selected`).
#' The rule, the thresholds and the votes that could not be evaluated are stored
#' in `metadata(tse)$viral_selection`.
#'
#' Because votes are combined with OR, this is a permissive union chosen for
#' sensitivity, not a consensus. It is also not a set of independent tests: in
#' ViroProfiler, geNomad, CheckV and VIBRANT already determined which contigs
#' entered the candidate set, VirSorter2 ran only on that candidate set, and the
#' merged taxonomy was computed from it. `rule = "candidate"` uses the pipeline's
#' own candidate list instead and is the more defensible criterion when
#' `create_vpftse()` was given `fin_vircontigs`.
#'
#' @param tse A TSE built by [create_vpftse()]
#' @param rule `"vote"` (default) for the OR of the votes below, or `"candidate"`
#'   for the pipeline's own putative viral contig list
#' @param votes Names of the votes to evaluate. Defaults to all of `taxonomy`,
#'   `checkv`, `virsorter2`, `vibrant`, `dvf`, `genomad`.
#' @param genomad_min_score Minimum geNomad `virus_score`. The default of 0.7 is
#'   the conventional cutoff for geNomad's uncalibrated score; geNomad's own
#'   recommendation when score calibration is enabled is to filter on FDR.
#' @param checkv_levels CheckV quality tiers counted as viral. CheckV quality
#'   grades completeness of sequences already submitted as putative viruses; it
#'   is evidence of assembly quality rather than of viral identity.
#' @param virsorter2_groups VirSorter2 `max_score_group` values counted as viral
#' @param candidate_col rowData column holding the upstream candidate flag
#' @return The same TSE with the vote columns added
#' @export
#' @importFrom SummarizedExperiment rowData
annotate_viral_votes <- function(tse,
                                 rule = c("vote", "candidate"),
                                 votes = names(.VPF_VOTES),
                                 genomad_min_score = 0.7,
                                 checkv_levels = c("Complete", "High-quality",
                                                   "Medium-quality"),
                                 virsorter2_groups = c("dsDNAphage", "NCLDV", "RNA",
                                                       "ssDNA", "lavidaviridae"),
                                 candidate_col = "upstream_viral_candidate") {
  rule <- match.arg(rule)
  rd <- SummarizedExperiment::rowData(tse)
  n <- nrow(rd)
  unknown <- setdiff(votes, names(.VPF_VOTES))
  if (length(unknown) > 0) {
    stop("Unknown vote(s): ", paste(unknown, collapse = ", "),
         ". Available: ", paste(names(.VPF_VOTES), collapse = ", "), call. = FALSE)
  }

  if (rule == "candidate") {
    if (!candidate_col %in% colnames(rd)) {
      stop("rowData has no '", candidate_col, "' column. Build the object with ",
           "create_vpftse(fin_vircontigs = ...) to record which contigs the ",
           "pipeline selected as putative viruses.", call. = FALSE)
    }
    selected <- as.logical(rd[[candidate_col]])
    if (length(selected) != n || anyNA(selected)) {
      stop("'", candidate_col, "' must be a logical vector of length ", n,
           " with no NA", call. = FALSE)
    }
    rd$viral_selected <- selected
    rd$viral_vote_n <- as.integer(selected)
    rd$viral_vote_evidence <- ifelse(selected, "upstream_candidate", NA_character_)
    SummarizedExperiment::rowData(tse) <- rd
    S4Vectors::metadata(tse)$viral_selection <- list(
      rule = "candidate", candidate_col = candidate_col,
      votes_used = character(0), votes_missing = character(0),
      thresholds = list(), n_total = n, n_selected = sum(selected),
      decided_at = Sys.time())
    return(tse)
  }

  vote_mat <- matrix(FALSE, nrow = n, ncol = 0)
  used <- character(0)
  missing_votes <- character(0)
  for (nm in votes) {
    spec <- .VPF_VOTES[[nm]]
    if (!spec$column %in% colnames(rd)) {
      missing_votes <- c(missing_votes, paste0(nm, " (", spec$column, ")"))
      next
    }
    x <- rd[[spec$column]]
    # A zero-length or wrong-length column would make the OR collapse to
    # logical(0) and silently select nothing.
    if (length(x) != n) {
      stop("rowData column '", spec$column, "' has length ", length(x),
           " but the object has ", n, " rows", call. = FALSE)
    }
    v <- switch(
      spec$kind,
      not_na = if (is.character(x)) !is.na(x) & nzchar(trimws(x)) else !is.na(x),
      in_levels = as.character(x) %in% switch(nm, checkv = checkv_levels,
                                              virsorter2 = virsorter2_groups),
      min_score = {
        s <- suppressWarnings(as.numeric(x))
        !is.na(s) & s >= genomad_min_score
      }
    )
    v[is.na(v)] <- FALSE
    vote_mat <- cbind(vote_mat, stats::setNames(data.frame(v), nm))
    used <- c(used, nm)
  }

  if (length(used) == 0L) {
    stop("None of the requested viral-identity votes could be evaluated: ",
         paste(missing_votes, collapse = ", "),
         ". Refusing to return an empty or unfiltered object. Check that the ",
         "upstream tables were joined and that their column names are unchanged.",
         call. = FALSE)
  }
  if (length(missing_votes) > 0) {
    warning("Viral-identity votes skipped because their rowData columns are absent: ",
            paste(missing_votes, collapse = ", "),
            ". The remaining votes were used: ", paste(used, collapse = ", "),
            call. = FALSE)
  }

  vote_mat <- as.matrix(vote_mat)
  selected <- rowSums(vote_mat) > 0
  for (nm in used) rd[[paste0("viral_vote_", nm)]] <- unname(vote_mat[, nm])
  rd$viral_vote_n <- as.integer(rowSums(vote_mat))
  rd$viral_vote_evidence <- vapply(seq_len(n), function(i) {
    hit <- used[vote_mat[i, ]]
    if (length(hit) == 0L) NA_character_ else paste(hit, collapse = ";")
  }, character(1))
  rd$viral_selected <- unname(selected)
  SummarizedExperiment::rowData(tse) <- rd

  S4Vectors::metadata(tse)$viral_selection <- list(
    rule = "vote",
    combination = "OR (permissive union, not a consensus)",
    votes_used = used,
    votes_missing = missing_votes,
    thresholds = list(genomad_min_score = genomad_min_score,
                      checkv_levels = checkv_levels,
                      virsorter2_groups = virsorter2_groups),
    n_total = n,
    n_selected = sum(selected),
    n_by_vote = colSums(vote_mat),
    decided_at = Sys.time()
  )
  tse
}


#' Create a TSE virome object
#'
#' Subsets a TSE built by [create_vpftse()] down to the contigs judged viral.
#' The judgement itself, and the record of which evidence produced it, come from
#' [annotate_viral_votes()]; see there for the rule, its thresholds and their
#' limitations.
#'
#' A TSE that has already been annotated keeps its existing `viral_selected`
#' column unless `reannotate = TRUE`.
#'
#' @param tse TSE object
#' @param rule `"vote"` (default) or `"candidate"`
#' @param reannotate Recompute the votes even if the object already carries them
#' @param ... Passed to [annotate_viral_votes()], e.g. `genomad_min_score`
#' @return TreeSummarizedExperiment object
#' @export
#'
create_vpftse_vir <- function(tse, rule = c("vote", "candidate"),
                              reannotate = FALSE, ...) {
  rule <- match.arg(rule)
  rd <- SummarizedExperiment::rowData(tse)
  if (reannotate || !"viral_selected" %in% colnames(rd)) {
    tse <- annotate_viral_votes(tse, rule = rule, ...)
    rd <- SummarizedExperiment::rowData(tse)
  }
  selected <- as.logical(rd$viral_selected)
  selected[is.na(selected)] <- FALSE
  if (length(selected) != nrow(rd)) {
    stop("Internal error: viral_selected has the wrong length", call. = FALSE)
  }
  if (!any(selected)) {
    warning("No contig was judged viral, so the returned object has zero rows.",
            call. = FALSE)
  }
  tse[selected, ]
}
