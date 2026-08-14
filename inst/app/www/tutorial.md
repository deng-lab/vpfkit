## ViroProfiler-viewer

An interactive viewer for the results of the
[ViroProfiler](https://github.com/deng-lab/viroprofiler) pipeline. It reads the
`TreeSummarizedExperiment` (TSE) object that `create_vpftse()` writes, usually
named `viroprofiler_output.rds`, and lets you filter contigs, inspect their
annotations, and compare samples without writing code.

Nothing is required to try it: the **Data** tab ships two bundled example
datasets that load with a single click.

## How to run it

### From R

```r
# install.packages("remotes")
remotes::install_github("deng-lab/vpfkit")
vpfkit::run_app()
```

`run_app()` accepts two arguments worth knowing about:

```r
vpfkit::run_app(
  max_upload_mb     = 500,  # Shiny's own default is only 5 MB
  allow_server_path = TRUE  # let the app read a path on this machine
)
```

`max_upload_mb` also reads the `VPFKIT_MAX_UPLOAD_MB` environment variable. If
the app sits behind a reverse proxy, that proxy may reject a large upload before
Shiny sees it, and no R setting can override it.

### With Docker

```sh
docker run --rm -p 3838:3838 denglab/viroprofiler-viewer
```

## The tabs

### Data

Load a dataset from one of three sources, then read what it actually contains.

Three panels matter here:

* **Annotations available** lists every annotation family the viewer knows
  about, whether this object carries it, and the ViroProfiler option that
  produces it. A run without `--use_iphop` has no host predictions, and this is
  where you see that rather than finding an empty plot later.
* **Assays** explains what each abundance matrix measures. Read it before
  choosing a metric; see the warning about `tmm` below.
* **Library sizes** shows total signal per sample. Strongly unequal libraries
  affect observed richness and every untransformed comparison.

**Adding sample metadata.** ViroProfiler's samplesheet is restricted to
`sample`, `fastq_1` and `fastq_2`, so study variables such as treatment group,
subject or timepoint can never arrive with the pipeline output. Upload them here
as a CSV, TSV or Excel table whose first column, or a column named `sample`,
`sample_name` or `sample_id`, holds the sample identifiers. Every other column
is added to `colData`. This is what enables grouped composition plots, group
tests and PERMANOVA; without it the Diversity tab correctly refuses to run them.

### Filter

Select the contigs every other tab will use. Controls appear only for
annotations this dataset actually carries, so a filter can never be applied
against a missing column.

Available where the annotation exists: contig length, CheckV quality tier,
CheckV completeness, geNomad virus score, VirSorter2 score, number of
viral-evidence votes, prevalence across samples, and total abundance.

Three behaviours are worth knowing:

* **Coverage masking.** Reads from one genome often map to a related contig over
  a short stretch, which looks like presence. Abundance is therefore set to zero
  in every sample where the contig's covered fraction is below the threshold.
  Missing coverage counts as absent. This is the single most important guard
  against reading "TPM > 0" as confident detection.
* **Viral-evidence votes.** When the object records which detector called each
  contig viral, the *Viral evidence* panel shows the distribution and the
  minimum-votes filter becomes available. The detectors are combined with OR, so
  the vote count is a sensitivity ranking rather than a consensus: a contig
  admitted on one low-quality VIBRANT call contributes to every diversity index
  exactly as much as one that all four detectors agreed on. Requiring at least
  two votes is the cheapest way to drop that noise.
* **The audit table** records how many contigs each step removed and lists the
  filters that were *not* applied because their annotation is absent. A step
  that silently drops everything is the failure mode this table exists to
  prevent.

### Taxonomy

Stacked composition at any rank that carries assignments, optionally faceted by
a metadata variable, in relative or absolute units.

Contigs with no assignment at the chosen rank are kept and labelled
`Unclassified`. In a virome this is usually the largest and most interesting
fraction; hiding it rescales every other bar and misrepresents the community.
You can hide it, but the app says so explicitly when you do.

The **Classification depth** panel shows how many contigs carry an assignment at
each rank. A steep fall toward Genus and Species is normal for viral taxonomy,
not a defect in the run.

### Diversity

Alpha diversity per sample, dissimilarities between samples, ordination, and
PERMANOVA.

Each procedure has a sample size below which its result is not merely
under-powered but undefined, and the app enforces those limits instead of
printing a number that cannot mean anything:

| Procedure | Enabled when | Why |
|---|---|---|
| Per-sample alpha diversity | any sample with signal | computed independently per sample |
| Chao1 | the active assay holds raw integer counts | it is estimated from singletons and doubletons |
| Pairwise dissimilarity | at least 2 samples | one distance is descriptive but not a pattern |
| PCoA, two axes | at least 3 samples | classical MDS has at most `n - 1` axes |
| NMDS, two axes | at least 6 samples | two dimensions are saturated below that, so stress is trivially near zero |
| Group test | at least 2 groups, at least 2 samples per group | one sample per group leaves no within-group variability |
| PERMANOVA | at least 4 samples, at least 2 groups, at least 2 per group | otherwise the design has no estimable residual |

Where a test is possible but its resolution is inadequate, the app says so. A
two-sided exact Wilcoxon test cannot reach `p < 0.05` until the allocation count
exceeds 40, first achieved by a balanced four-versus-four design; a permutation
test cannot reach `p < 0.05` unless more than 20 distinguishable permutations
exist.

Note also that these values describe a *selected subset*. The viral TSE is
already the subset of the assembly that passed the pipeline's viral-evidence
rules, and your filters narrow it further. Diversity values are comparable
between samples processed identically; they are not absolute community
diversities.

### Contigs

Every per-contig annotation column, searchable and downloadable, with a
drill-down for a single contig showing its annotations and its abundance in each
sample. The **Column overview** panel reports how much of each column is filled,
which is the quickest way to see which pipeline stages produced usable output.

Large tables switch to server-side mode automatically: the rows stay in R and
only the visible page travels to the browser, while search and sort still apply
to the whole table.

### Host & lifestyle

Predicted host taxa from iPHoP or PHIST, predicted replication cycle from
BACPHLIP, Replidec or VIBRANT, and provirus status from CheckV and geNomad.

A host prediction links a viral sequence to a candidate host by sequence
evidence. It is not a demonstration that the phage infects that host, and the
panel says so rather than leaving the impression that it is an observation.

### Genes

Gene-level annotations read from `metadata(tse)$gene_annotations`, restricted to
the contigs that pass the current filters. Which annotator produced them varies
by run, so each panel takes whichever column is present rather than assuming one
naming scheme:

| Panel | Prefers | Falls back to |
|---|---|---|
| Functional categories | `pharokka_category` | `checkamg_phrog` |
| Auxiliary genes | `dramv_ko`, restricted to DRAM-v's AMG calls | `checkamg_kegg_ko` |
| CARD, VFDB | `pharokka_card`, `pharokka_vfdb` | - |
| CheckAMG | `checkamg_class`, `checkamg_function` | `checkamg_kegg_ko`, `checkamg_pfam`, `checkamg_cazy` |

CheckAMG separates auxiliary genes into metabolic (AMG), physiological (APG) and
regulatory (AReG) classes and records whether each protein lies inside a strict
viral region. Proteins outside one are flagged, because an auxiliary-gene call
there may belong to residual host sequence rather than to the virus.

A CARD or VFDB match is a sequence similarity hit. It may lie in residual host
sequence on a contaminated contig, and it is not a demonstrated resistance or
virulence phenotype. Check the contig's CheckV contamination before drawing a
conclusion from one.

### Compare

Load a second dataset and compare shared and unique taxa, per-sample diversity,
and composition at a rank both objects carry.

No between-dataset significance test is offered. Two runs differ in assembly,
sequencing depth and filtering as much as in biology, and contig identifiers are
assigned per assembly, so two independent runs share none unless both were
dereplicated against a common vOTU catalogue.

### Export

Downloads of the filtered TSE, the abundance matrix, contig annotations, sample
metadata and gene annotations, all reflecting the current selection.

The **provenance** download is a plain-text record of the dataset, the active
assay, every filter step with its effect, the meaning of each assay, and the
package versions used. Keep it with the exported tables: it is what makes a
figure reproducible months later.

## Reading the abundance metrics correctly

ViroProfiler runs [CoverM](https://github.com/wwood/CoverM), which reports
several different quantities. They are not interchangeable, and the assay names
alone do not say what they measure.

| Assay | What it is | Use it for |
|---|---|---|
| `counts` | number of reads mapped to the contig | anything, and the only valid input for count-based models |
| `tpm` | reads per kilobase scaled to a fixed library sum | composition, diversity, ordination |
| `trimmed_mean` | mean per-base coverage **depth**, in x-fold, after discarding the 5% highest and 5% lowest covered positions | composition and ordination, with its units in mind |
| `covfrac` | fraction of the contig covered by at least one read, 0 to 1 | deciding whether a contig is present at all |

Two consequences the viewer enforces:

* **`covfrac` is never offered as an abundance.** It measures breadth of
  coverage, not amount. Using it for composition or diversity produces numbers
  that look fine and mean nothing.
* **An assay named `tmm` is not edgeR's TMM.** Older ViroProfiler runs stored
  CoverM's trimmed mean under the name `tmm`, which in bioinformatics normally
  denotes the trimmed mean of M-values normalization from edgeR. It is a
  coverage depth. Newer runs name it `trimmed_mean`. Reading an older object
  renames the assay in memory and says so on the Data tab; the file on disk is
  untouched.

Where the object carries its own assay documentation, written by
`create_vpftse()`, the Data tab shows that text rather than the description
above, because it comes from the code that produced the matrix.

## Getting help

Report problems at
[github.com/deng-lab/vpfkit/issues](https://github.com/deng-lab/vpfkit/issues).
Including the provenance file from the Export tab makes a report much faster to
act on.
