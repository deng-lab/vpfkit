# Packaging changes required by the reader/TSE audit

Everything below has to happen outside `R/fct_readvpf.R`, `R/data.R`,
`tests/testthat/test-fct_readvpf.R` and `tests/testthat/fixtures/`, which is why it is
listed here instead of being applied directly.

## 1. `devtools::document()`

Four functions are newly exported with `@export` roxygen tags and need to reach
`NAMESPACE` and `man/`:

| Function | Purpose |
|---|---|
| `read_checkamg()` | Reads `checkamg_results/final_results.tsv`, or the directory the `CHECKAMG` process publishes. There was no CheckAMG reader at all. |
| `read_coverm_log()` | Parses the per-sample library size and mapping rate out of a CoverM log, for `colData`. |
| `normalize_assay_names()` | Renames the legacy `tmm` assay to `trimmed_mean` when an object saved by an earlier version is read. |
| `annotate_viral_votes()` | Writes the per-contig viral-identity votes into `rowData` and the rule into `metadata()`. `create_vpftse_vir()` calls it, so a caller that wants the audit trail on the *all-contig* object has to call it directly. |

No function was removed and no existing export changed name.

## 2. `DESCRIPTION`: two packages to add to `Imports`

Both are base R distribution packages, but `R CMD check` fails a `::` call into a
package that is not declared.

```
Imports:
    stats,
    utils,
```

- `stats::na.omit()`, `stats::setNames()` — used when validating tool vocabularies and
  when building the vote matrix.
- `utils::head()` — used to truncate the identifier lists inside error messages.

No new third-party dependency was introduced. `openxlsx` (already in `Imports`, and
already in `NAMESPACE` via `importFrom(openxlsx, read.xlsx)`) is what reads an `.xlsx`
sample metadata file.

## 3. `R/mod_vpfilter.R:130` — required, or the Shiny app loses an assay

The assay that holds CoverM's trimmed mean of per-base coverage depth is now written as
`trimmed_mean`. `mod_vpfilter.R` still maps the dropdown label to the old name:

```r
"trimmed mean" = "tmm",
```

That has to become `"trimmed mean" = "trimmed_mean"`. The module's fallback
(`if (!assay_name %in% assayNames(tse_obj)) assay_name <- abdc_metric()`) does not help,
because `abdc_metric()` returns the same stale `"tmm"`.

For objects loaded from an existing `.rds`, the module should call
`normalize_assay_names(tse)` right after `readRDS()`; it renames `tmm` in place and warns
once. Writing only ever produces the new name.

The same module can now consume, if useful:

- `rowData(tse)$viral_vote_*`, `viral_vote_n`, `viral_vote_evidence`, `viral_selected`
- `metadata(tse)$viral_selection` — the rule, the thresholds, and the per-vote tallies
- `metadata(tse)$viroprofiler$join_match` — how many contigs each tool actually matched
- `metadata(tse)$gene_annotations_by_tool` — the gene tables kept separate, alongside the
  stacked `metadata(tse)$gene_annotations` the module already reads (which now carries a
  `source` column)
- `colData(tse)$n_reads_total`, `n_reads_mapped`, `mapping_rate`, plus whatever columns
  the user's metadata file contributed

## 4. `tests/testthat/test-fct_export.R` and `test-fct_report.R`

Both build their test TSE from `fixtures/mmseqs2_taxa.tsv`. MMseqs2 taxonomy was retired
from ViroProfiler; the file that actually feeds `fin_taxa` is `taxonomy_tse.tsv`, written
by `bin/merge_taxonomy.py`. The fixture has been rewritten to the real `taxonomy_tse.tsv`
schema and is now present under **both** names so that neither test breaks. Once those two
files are updated to `fixtures/taxonomy_tse.tsv`, `fixtures/mmseqs2_taxa.tsv` can be
deleted.

Renamed or added fixtures, for reference:

| Old | New | Why |
|---|---|---|
| `mmseqs2_taxa.tsv` | `taxonomy_tse.tsv` | the real producer is `merge_taxonomy.py`, not MMseqs2 |
| `iphop_genus.tsv` | `iphop_Host_prediction_to_genus_m90.csv` | iPHoP writes CSV, not TSV |
| `vrhyme_bins.tsv` | `vrhyme_best_bins_membership.tsv` | real header is `scaffold`, `bin` |
| — | `checkamg_final_results.tsv` | new reader |
| — | `virsorter2_final_viral_score_multigroup.tsv` | VirSorter2's score columns vary with `--include-groups` |

## 5. `data-raw/test_viroprofiler.rds` and `test_viroprofiler_2.rds`

Both synthetic objects still carry the `tmm` assay name. `dev/make_test_data.R` builds
them, so regenerating them is the clean fix; until then the Shiny app has to call
`normalize_assay_names()` after `readRDS()`, which is the same thing it needs for any
object a user produced with an earlier release.
