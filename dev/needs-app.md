# Shiny app: what the maintainer still has to do

Everything below is outside the files the app work owns. Nothing here is
optional if the package is to pass `R CMD check` after `devtools::document()`.

## 1. Run `devtools::document()`

New `@noRd` roxygen blocks were added throughout `R/mod_*.R`,
`R/golem_utils_ui.R` and `R/golem_utils_server.R`, and `run_app()` gained two
documented arguments. Only `run_app()` is exported, so `man/run_app.Rd` is the
one file that changes.

```r
devtools::document()
```

`run_app()`'s new arguments:

| Argument | Default | Purpose |
|---|---|---|
| `max_upload_mb` | `500`, or `VPFKIT_MAX_UPLOAD_MB` | Raises Shiny's 5 MB upload ceiling, which is smaller than a real result object. |
| `allow_server_path` | `TRUE` outside production | Whether the app may read a dataset from a path on the host. |

## 2. DESCRIPTION

### Move to `Suggests`

```
Suggests:
    callr,
    chromote,
    knitr,
    markdown,
    quarto,
    rmarkdown,
    shinytest2,
    testthat (>= 3.0.0),
    withr
```

| Package | Why | If absent |
|---|---|---|
| `markdown` | `shiny::includeMarkdown()` renders the About tab. | `vpf_render_markdown()` falls back to `commonmark`, then to a `<pre>` block. The tab still works. |
| `shinytest2`, `chromote` | The browser end-to-end test in `test-app.R`. | The test skips. |
| `callr` | The HTTP-launch test in `test-app.R`. | The test skips. |
| `withr` | Already in `Suggests`; used by the new module tests. | - |

### Already in `Imports` and now used more heavily

`grDevices` is used by `vpf_palette()` for `colorRampPalette()`, and `tools`
by `vpf_read_metadata_file()` and the export filenames. Both are base R
packages that ship with every installation, but `R CMD check` wants them
declared:

```
Imports:
    ...
    grDevices,
    methods,
    stats,
    tools,
    utils,
```

`methods` is needed by `vpf_read_tse()` (`methods::is()`), `stats` and `utils`
throughout.

### Not added, deliberately

- **`bslib`** - Codex confirmed that `shinythemes::shinytheme("flatly")` plus
  `navbarPage()` is not deprecated, and the brief said to keep it. `bslib`
  would bring Bootstrap 5, `card()`, `value_box()` and dark mode, but it is a
  visual rewrite with real migration traps (Bootstrap 3 selectors, fillable
  flex layouts, `navbar_options()` since bslib 0.9). Worth doing as its own
  change, not folded into a bug-fix pass.
- **`DT`** - CRAN `reactable` 0.4.5 has no server-side data support; only the
  unreleased 0.4.5.9000 does, and it needs V8. Rather than add `DT` as a second
  table dependency, `mod_bigtable` implements server-side paging directly:
  filtering, sorting and slicing all happen in R and only the visible page is
  serialized. If a future dataset makes that insufficient, `DT` with
  `server = TRUE` is the drop-in replacement.
- **`shinyWidgets`** - the app uses only base Shiny inputs. The single control
  that needed more than `selectInput` could give, the contig picker on the
  Contigs tab, is a `selectizeInput` updated with
  `updateSelectizeInput(server = TRUE)`, so a hundred-thousand-contig object
  never ships its identifiers into the page.
- **`mirai` / `future` / `promises`** - `shiny::ExtendedTask` is the current
  recommendation for genuinely slow work. On the datasets this viewer targets
  (thousands of contigs, tens of samples) ordination and PERMANOVA finish in
  under a second, so the app uses `withProgress()` instead. Revisit if a user
  reports a blocked session.

## 3. Ship the demo datasets inside the package

The Data tab offers two bundled example datasets so the app is usable with no
files at all. They are found by `vpf_demo_datasets()`, which searches
`inst/extdata` first and falls back to `data-raw/` for a source checkout. The
fallback is what makes them work under `pkgload::load_all()`; an installed
package needs them copied:

```r
dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
file.copy(
  c("data-raw/test_viroprofiler.rds", "data-raw/test_viroprofiler_2.rds"),
  "inst/extdata/", overwrite = TRUE
)
```

Roughly 32 kB together. Without this step the Data tab's demo option is empty
in an installed package and several tests skip.

## 4. A latent bug in `R/utils.R`, which the app no longer depends on

`refind_abundance()` and `read_bins()` call `assay()`, `assay<-`, `fread()`,
`column_to_rownames()` and `left_join()` without qualification. Those names
resolve today only because `R/utils.R` begins with

```r
library(tidyverse)
library(data.table)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
```

Top-level `library()` calls in a package's `R/` directory run at **install**
time, not at load time. Under `pkgload::load_all()` they run in the current
session and mask the problem; in an installed package they do not, and neither
`assay` nor `assay<-` is in vpfkit's imports environment:

```r
imp <- parent.env(asNamespace("vpfkit"))
exists("assay", envir = imp, inherits = FALSE)     # FALSE
exists("assay<-", envir = imp, inherits = FALSE)   # FALSE
```

So `refind_abundance()` fails at runtime in an installed package unless the
user happens to have `SummarizedExperiment` attached.

The app no longer calls it: `vpf_mask_by_covfrac()` in
`R/golem_utils_server.R` does the same masking with every call fully
qualified, and additionally treats a missing covered fraction as absent rather
than letting `NA` propagate into the abundance matrix. `R/utils.R` was not
touched because it belongs to another work stream, but the `library()` calls
should be replaced with `@importFrom` tags or `pkg::` prefixes.

## 5. `.Rbuildignore`

`dev/` and `data-raw/` are already ignored. If step 3 is taken, `inst/extdata`
must **not** be ignored.

## 6. What the app expects from the reader side

Consumed where present, degraded gracefully where not:

| Field | Used by | Behaviour when absent |
|---|---|---|
| `metadata(tse)$viroprofiler$assays` | Data tab assay descriptions | Falls back to the app's own dictionary. |
| `metadata(tse)$viral_selection` | Filter tab, Viral evidence | Panel explains that `annotate_viral_votes()` records it. |
| `rowData$viral_vote_n`, `viral_vote_*` | Minimum-votes filter, Diversity caveat | Filter control is not rendered; audit records "not applied". |
| `metadata(tse)$gene_annotations` | Genes tab | Names `--use_dram` and pharokka as the source. |
| `metadata(tse)$gene_annotations_by_tool` | Genes tab fallback | Used only when the merged table is absent. |
| `rowData$checkamg_*` | Genes tab, CheckAMG panel | Names `--use_checkamg`. |
| `rowData$iphop_*`, `phist_*` | Host tab | Names `--use_iphop`. |
| `rowData$bacphlip_replicyc` and friends | Lifestyle tab | Names `--replicyc`. |
| `colData` grouping columns | Diversity, Taxonomy | Every group test stays disabled and says why. |

`normalize_assay_names()` is called inside `vpf_read_tse()`, so the app sees
`trimmed_mean` whether the object on disk uses that name or the legacy `tmm`.

## 7. Files deleted from `inst/app/www/`

`custom.sass` (one space), `custom2.css` (an rmarkdown TOC rule this app never
loads), `handlers.js` (an empty custom message handler) and `script.js` (an
empty ready handler) were dead. `custom.css` was empty and now holds the app's
styling; `footer.html` is new and is what the About tab includes.
