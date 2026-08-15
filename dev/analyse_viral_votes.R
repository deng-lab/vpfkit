# What does the viral vote actually select?
#
#   Rscript dev/analyse_viral_votes.R [path/to/object_all_contigs.rds]
#
# annotate_viral_votes() combines its votes with OR, so the viral set is a union
# of the detectors' sensitivities. Two things about that are worth measuring
# rather than assuming, and both need an object from *before* the viral subset
# was taken -- `*_all_contigs.rds` -- because the subset has already dropped
# every contig the union rejected.
#
#   1. How far the union is from `rule = "candidate"`, the pipeline's own list.
#      This is the cost of changing the default.
#   2. Whether the votes are independent. They are not: geNomad, CheckV and
#      VIBRANT decide the candidate set, VirSorter2 runs only on that set, and
#      the merged taxonomy is computed from it. Pairwise Jaccard puts a number
#      on it, which matters because `viral_vote_n` reads like a count of
#      independent confirmations and is not one.
suppressMessages(pkgload::load_all(rprojroot::find_root(rprojroot::is_r_package),
                                   quiet = TRUE))
suppressMessages(library(SummarizedExperiment))

args <- commandArgs(trailingOnly = TRUE)
path <- if (length(args) > 0) args[[1]] else Sys.getenv(
  "VPFKIT_VOTE_TSE",
  "/mnt/nas26/testdata/viroprofiler_16sample/run_final/results/viroprofiler_output_all_contigs.rds"
)
if (!file.exists(path)) stop("No object at ", path)

tse <- readRDS(path)
rd <- as.data.frame(rowData(tse))
cat("object   :", basename(path), "\ncontigs  :", nrow(rd), "\n\n")

vote_cols <- setdiff(grep("^viral_vote_", names(rd), value = TRUE),
                     c("viral_vote_n", "viral_vote_evidence"))
if (length(vote_cols) == 0) {
  stop("No viral_vote_* columns. Was this object built by annotate_viral_votes()?")
}

cat("=== votes in favour, per tool ===\n")
per <- vapply(vote_cols, function(c) sum(rd[[c]], na.rm = TRUE), integer(1))
names(per) <- sub("^viral_vote_", "", names(per))
print(sort(per, decreasing = TRUE))

cat("\n=== viral_vote_n (0 = not selected) ===\n")
print(table(rd$viral_vote_n))

cat("\n=== contigs resting on exactly one vote ===\n")
print(table(rd$viral_vote_evidence[rd$viral_vote_n == 1]))

if ("upstream_viral_candidate" %in% names(rd)) {
  uni <- rd$viral_vote_n >= 1
  cand <- as.logical(rd$upstream_viral_candidate)
  cat("\n=== union (rule = \"vote\") vs the pipeline candidate list ===\n")
  print(table(union = uni, candidate = cand))
  cat(sprintf("\nunion only    : %d\ncandidate only: %d\n",
              sum(uni & !cand), sum(cand & !uni)))
} else {
  cat("\n(no upstream_viral_candidate column: create_vpftse() had no fin_vircontigs)\n")
}

cat("\n=== pairwise Jaccard between votes ===\n")
m <- as.matrix(rd[, vote_cols, drop = FALSE])
colnames(m) <- sub("^viral_vote_", "", colnames(m))
m[is.na(m)] <- FALSE
jac <- outer(seq_len(ncol(m)), seq_len(ncol(m)), Vectorize(function(i, j) {
  u <- sum(m[, i] | m[, j])
  if (u == 0) NA_real_ else round(sum(m[, i] & m[, j]) / u, 3)
}))
dimnames(jac) <- list(colnames(m), colnames(m))
print(jac)
cat("\nHigh off-diagonal values are expected and are the point: these votes",
    "\nshare an upstream, so they are not independent evidence.\n")
