# Do the viral-vote findings hold on a second, independent dataset?
#
# The reference run measured this over 161 contigs, of which 100 were selected.
# p0075 (pig CRC virome, 73 samples) offers a candidate pool of ~17,000 contigs
# of which 5,150 were selected -- two orders of magnitude larger, a different
# host, a different study, and an older ViroProfiler whose detector tables use
# the older column names, mapped here onto the same five vote semantics.
#
# The pool, not the viral subset. Running this on p0075's tse_vir_contig.rds
# instead gives Jaccard 0.94-1.00 and no single-vote contigs, which is not a
# result: that object has already dropped everything the union rejected, so
# every detector necessarily agrees on what is left.
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
D <- if (length(args) > 0) args[[1]] else Sys.getenv(
  "VPFKIT_VOTE_TABLES",
  "/home/allen/github/rujinlong/p0075-pigCRCvirome/analyses/data/00-raw/d101-import_virome"
)
if (!dir.exists(D)) stop("No detector tables at ", D)

checkv  <- fread(cmd = paste("zcat", file.path(D, "quality_summary.tsv.gz")), showProgress = FALSE)
vs2     <- fread(file.path(D, "virsorter2-score.tsv"), showProgress = FALSE)
vibrant <- fread(file.path(D, "VIBRANT_genome_quality.tsv"), showProgress = FALSE)
genomad <- fread(file.path(D, "genomad_virus_summary.tsv"), showProgress = FALSE)

# The pool every detector was run against: CheckV sees every candidate contig.
pool <- unique(checkv$contig_id)
cat(sprintf("candidate pool: %d contigs\n\n", length(pool)))

CHECKV_LEVELS <- c("Complete", "High-quality", "Medium-quality")
VS2_GROUPS <- c("dsDNAphage", "NCLDV", "RNA", "ssDNA", "lavidaviridae")

hit <- list(
  checkv     = checkv[checkv_quality %in% CHECKV_LEVELS, contig_id],
  virsorter2 = vs2[max_score_group %in% VS2_GROUPS, seqname],
  vibrant    = unique(vibrant$scaffold),
  genomad    = genomad[virus_score >= 0.7, seq_name]
)
m <- sapply(hit, function(ids) pool %in% ids)
rownames(m) <- pool

cat("=== votes in favour, per tool (of the pool) ===\n")
print(sort(colSums(m), decreasing = TRUE))

n <- rowSums(m)
cat("\n=== vote_n distribution ===\n"); print(table(n))
cat(sprintf("\nunion (>=1 vote): %d of %d (%.1f%%)\n",
            sum(n >= 1), length(pool), 100 * sum(n >= 1) / length(pool)))

cat("\n=== contigs resting on exactly one vote ===\n")
print(table(colnames(m)[apply(m[n == 1, , drop = FALSE], 1, which)]))

cat("\n=== pairwise Jaccard ===\n")
k <- ncol(m)
jac <- outer(seq_len(k), seq_len(k), Vectorize(function(i, j) {
  u <- sum(m[, i] | m[, j]); if (u == 0) NA_real_ else round(sum(m[, i] & m[, j]) / u, 3)
}))
dimnames(jac) <- list(colnames(m), colnames(m))
print(jac)

core <- c("virsorter2", "vibrant", "genomad")
cc <- jac[core, core]
cat(sprintf("\ncore detectors (%s) pairwise: %.3f - %.3f\n", paste(core, collapse = "/"),
            min(cc[upper.tri(cc)]), max(cc[upper.tri(cc)])))
cat(sprintf("checkv against the other three   : %.3f - %.3f\n",
            min(jac["checkv", core]), max(jac["checkv", core])))
cat("\nreference run, same measurement: core 0.79-0.87, checkv 0.09-0.11\n")
