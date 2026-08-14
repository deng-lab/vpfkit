#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom data.table fread
#' @importFrom data.table setnames
#' @importFrom dplyr across where distinct
#' @importFrom dplyr any_of
#' @importFrom dplyr case_when
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr inner_join
#' @importFrom dplyr mutate
#' @importFrom dplyr n
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr summarise
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggrepel geom_text_repel
#' @importFrom openxlsx read.xlsx
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom reactable colDef
#' @importFrom reactable reactable
#' @importFrom reactable reactableOutput
#' @importFrom reactable renderReactable
#' @importFrom rlang .data
#' @importFrom stringr str_detect
#' @importFrom stringr str_replace
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_to_title
#' @importFrom SummarizedExperiment rowData
#' @importFrom tibble column_to_rownames
#' @importFrom tibble rownames_to_column
## usethis namespace: end
NULL

#' Dependencies that are real but invisible to `R CMD check`
#'
#' Two packages in `Imports` are never named by any R source file here, and
#' would otherwise be reported as unused and eventually dropped:
#'
#' * `R.utils` is what `data.table::fread()` delegates gz decompression to.
#'   Without it a clean installation cannot read a single one of ViroProfiler's
#'   abundance tables, all of which are `.tsv.gz`.
#' * `here` is evaluated by `inst/golem-config.yml`, which golem reads at
#'   startup; the failure without it is `there is no package called 'here'`
#'   from inside `golem::get_golem_options()`.
#'
#' @return `NULL`, invisibly. Never called.
#' @noRd
vpf_declared_but_uncalled_imports <- function() {
  R.utils::gunzip
  here::here
  invisible(NULL)
}
