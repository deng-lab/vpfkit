## Internal utilities shared across vpfkit.
##
## Nothing in this file is exported. Everything is namespace-qualified on
## purpose: a package must never call library() at load time.

# ---------------------------------------------------------------------------
# Assay naming
# ---------------------------------------------------------------------------

## ViroProfiler's ABUNDANCE process runs `coverm contig --methods ...` six
## times and vpfkit stores four of those tables as assays. The names below are
## canonical; the aliases are names the same quantity has been stored under by
## earlier versions of the package, so that objects saved with either spelling
## keep working.
##
## `tmm` is a legacy alias for `trimmed_mean` and is *not* edgeR's Trimmed Mean
## of M-values. It is CoverM's trimmed mean of the per-base coverage depth.
.vpf_assay_aliases <- list(
  counts        = c("counts", "count"),
  tpm           = c("tpm", "TPM"),
  trimmed_mean  = c("trimmed_mean", "tmm", "trimmed.mean"),
  covfrac       = c("covfrac", "covered_fraction", "coverage_fraction"),
  reads_per_base = c("reads_per_base", "rpb"),
  rpkm          = c("rpkm", "RPKM")
)

## One-line description of what each quantity actually is, used by the report
## so that no figure shows an abundance without saying which one.
.vpf_assay_units <- c(
  counts         = "CoverM read count (reads mapped per contig per sample)",
  tpm            = "CoverM TPM (transcripts per million; column sums to 1e6)",
  trimmed_mean   = "CoverM trimmed mean of per-base coverage depth (fold)",
  covfrac        = "CoverM covered fraction (breadth of coverage, 0-1)",
  reads_per_base = "CoverM reads per base (reads mapped / contig length)",
  rpkm           = "CoverM RPKM (reads per kilobase per million mapped reads)"
)

#' Legacy assay names and what each should be renamed to
#'
#' Derived from `.vpf_assay_aliases` so that a spelling added there is honoured
#' both when resolving a requested assay and when rewriting a stored one. A
#' top-level constant would not do: R collates `fct_readvpf.R`, where the only
#' caller lives, before this file.
#'
#' @return Named character vector, legacy name -> canonical name.
#' @noRd
vpf_legacy_assay_renames <- function() {
  out <- character(0)
  for (canonical in names(.vpf_assay_aliases)) {
    legacy <- setdiff(.vpf_assay_aliases[[canonical]], canonical)
    if (length(legacy)) out[legacy] <- canonical
  }
  out
}

#' Canonical name of an assay, following legacy aliases
#'
#' @param name Assay name as the caller spelled it.
#' @return The canonical name, or `name` unchanged when it is not a known
#'   ViroProfiler quantity.
#' @noRd
vpf_canonical_assay <- function(name) {
  if (!is.character(name) || length(name) != 1L || is.na(name)) return(name)
  for (canonical in names(.vpf_assay_aliases)) {
    if (name %in% .vpf_assay_aliases[[canonical]]) return(canonical)
  }
  name
}

#' Human-readable description of an assay
#'
#' @param name Assay name (canonical or alias).
#' @return A one-line description, or a generic fallback for unknown assays.
#' @noRd
vpf_assay_description <- function(name) {
  canonical <- vpf_canonical_assay(name)
  if (canonical %in% names(.vpf_assay_units)) {
    unname(.vpf_assay_units[canonical])
  } else {
    "user-supplied assay (units unknown to vpfkit)"
  }
}

#' Resolve an assay name against the assays an object actually has
#'
#' Accepts either the canonical name or any legacy alias, in either direction:
#' asking for `"trimmed_mean"` finds an assay stored as `"tmm"`, and asking for
#' `"tmm"` finds one stored as `"trimmed_mean"`.
#'
#' @param x A `SummarizedExperiment`, or a character vector of assay names.
#' @param name Assay name requested by the caller.
#' @param required When `TRUE` (default) an unresolvable name is an error;
#'   when `FALSE` the function returns `NA_character_`.
#' @return The name as stored in `x`, or `NA_character_`.
#' @noRd
vpf_match_assay <- function(x, name, required = TRUE) {
  available <- if (is.character(x)) x else SummarizedExperiment::assayNames(x)
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`assay.type` must be a single non-empty assay name.", call. = FALSE)
  }
  if (name %in% available) return(name)

  canonical <- vpf_canonical_assay(name)
  candidates <- .vpf_assay_aliases[[canonical]]
  hit <- if (is.null(candidates)) character(0) else intersect(candidates, available)
  if (length(hit) >= 1L) return(hit[[1L]])

  if (!required) return(NA_character_)
  stop(
    sprintf(
      "Assay '%s' not found. Available assays: %s.",
      name,
      if (length(available)) paste(sprintf("'%s'", available), collapse = ", ") else "none"
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Filesystem
# ---------------------------------------------------------------------------

#' Validate an output path and make sure its parent directory exists
#'
#' Writing to a directory that does not exist otherwise fails deep inside
#' `file()` with a bare "cannot open the connection", which does not say what
#' the caller did wrong.
#'
#' @param file Output path.
#' @param create_dir Create the parent directory when missing.
#' @return The path, expanded and made absolute.
#' @noRd
vpf_prepare_outfile <- function(file, create_dir = TRUE) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stop("`file` must be a single non-empty file path.", call. = FALSE)
  }
  file <- path.expand(file)
  if (dir.exists(file)) {
    stop(sprintf("`file` points at an existing directory: '%s'.", file), call. = FALSE)
  }
  parent <- dirname(file)
  if (!dir.exists(parent)) {
    if (!create_dir) {
      stop(sprintf("Output directory does not exist: '%s'.", parent), call. = FALSE)
    }
    if (!dir.create(parent, recursive = TRUE, showWarnings = FALSE)) {
      stop(sprintf("Could not create output directory: '%s'.", parent), call. = FALSE)
    }
  }
  if (file.access(parent, mode = 2L) != 0L) {
    stop(sprintf("Output directory is not writable: '%s'.", parent), call. = FALSE)
  }
  if (file.exists(file) && file.access(file, mode = 2L) != 0L) {
    stop(sprintf("Output file exists and is not writable: '%s'.", file), call. = FALSE)
  }
  ## Absolute paths survive any working-directory change made downstream.
  if (!isTRUE(startsWith(file, "/"))) {
    file <- file.path(normalizePath(parent, winslash = "/", mustWork = TRUE), basename(file))
  }
  file
}

#' Assert that an object is a SummarizedExperiment
#'
#' @param tse Object to check.
#' @param arg Argument name to quote in the error message.
#' @return `TRUE`, invisibly.
#' @noRd
vpf_assert_se <- function(tse, arg = "tse") {
  if (!methods::is(tse, "SummarizedExperiment")) {
    stop(
      sprintf("`%s` must be a SummarizedExperiment (got %s).", arg, class(tse)[[1L]]),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Abundance masking
# ---------------------------------------------------------------------------

#' Zero out abundances of contigs whose coverage breadth is too low
#'
#' A contig can collect reads over a short conserved stretch while the rest of
#' its length has no coverage at all. Requiring a minimum *breadth* of coverage
#' (the fraction of the contig covered at least once) before an abundance is
#' believed removes most of those false positives. The widely used convention
#' is 75% breadth combined with a strict identity threshold (Roux et al., 2017,
#' PeerJ 5:e3817, doi:10.7717/peerj.3817); the ViroProfiler paper uses 50%
#' (Ru et al., 2023, Gut Microbes, doi:10.1080/19490976.2023.2192522).
#'
#' Masking should be applied to raw counts *before* normalization. Masking an
#' already-normalized assay such as TPM leaves each sample with its own
#' residual scaling factor, so columns no longer sum to a common total
#' (Kuzub et al., 2025, Nat Commun, doi:10.1038/s41467-025-61478-7). This
#' function warns when it is asked to mask a normalized assay.
#'
#' @param tse A `SummarizedExperiment`.
#' @param raw_abundance_name Assay to mask. Legacy assay names are accepted.
#' @param covfrac_name Assay holding the covered fraction.
#' @param thr Minimum covered fraction, in `[0, 1]`. Contigs at exactly `thr`
#'   are kept. The default follows Roux et al. (2017); pass 0.5 to reproduce
#'   analyses made with the ViroProfiler paper's threshold.
#' @return `tse`, with the masked values written back into the same assay.
#' @noRd
refind_abundance <- function(tse, raw_abundance_name, covfrac_name = "covfrac", thr = 0.75) {
  vpf_assert_se(tse)
  if (!is.numeric(thr) || length(thr) != 1L || is.na(thr)) {
    stop("`thr` must be a single non-missing number.", call. = FALSE)
  }
  if (thr < 0 || thr > 1) {
    stop(sprintf("`thr` must be between 0 and 1 (got %s).", format(thr)), call. = FALSE)
  }

  abundance_name <- vpf_match_assay(tse, raw_abundance_name)
  covfrac_stored <- vpf_match_assay(tse, covfrac_name, required = FALSE)
  if (is.na(covfrac_stored)) {
    warning(
      sprintf("Assay '%s' not found; abundances left unmasked.", covfrac_name),
      call. = FALSE
    )
    return(tse)
  }

  abundance <- as.matrix(SummarizedExperiment::assay(tse, abundance_name))
  covfrac <- as.matrix(SummarizedExperiment::assay(tse, covfrac_stored))

  ## Both matrices come from the same object, so they are aligned by
  ## construction. Check anyway: element-wise multiplication in R aligns by
  ## position and never by name, so a mismatch would silently mask the wrong
  ## contigs rather than fail.
  if (!identical(dim(abundance), dim(covfrac))) {
    stop("Abundance and covered-fraction assays have different dimensions.", call. = FALSE)
  }

  if (vpf_canonical_assay(abundance_name) %in% c("tpm", "rpkm")) {
    warning(
      sprintf(
        paste0(
          "Masking '%s', which is already normalized. Column sums will no longer ",
          "be comparable across samples. Mask raw counts and re-normalize instead."
        ),
        abundance_name
      ),
      call. = FALSE
    )
  }

  ## A missing covered fraction means the contig was never observed to be
  ## covered, so it is treated as absent rather than propagated as NA.
  keep <- !is.na(covfrac) & covfrac >= thr
  abundance[!keep] <- 0

  SummarizedExperiment::assay(tse, abundance_name) <- abundance
  tse
}
