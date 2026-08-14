#' The application server-side
#'
#' The loaded object and the filtered object are each computed once, in
#' `mod_data_input_server()` and `mod_vpfilter_server()`, and handed to every
#' other module as reactives. Nothing downstream re-reads the file or re-applies
#' the filters, so all tabs necessarily describe the same contigs.
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  r_input <- mod_data_input_server("data")
  r_filter <- mod_vpfilter_server("filter", r_input)

  mod_composition_server("composition", r_filter)
  mod_diversity_server("diversity", r_filter)
  mod_features_server("features", r_filter)
  mod_host_server("host", r_filter)
  mod_genes_server("genes", r_filter)
  mod_compare_server("compare", r_filter)
  mod_export_server("export", r_filter)
  mod_about_server("about")

  output$status_bar <- renderText({
    raw <- r_filter$raw()
    if (is.null(raw)) {
      return("No dataset loaded. Open the Data tab to load a bundled example or your own result.")
    }
    keep <- r_filter$tse()
    n_keep <- if (is.null(keep)) 0L else nrow(keep)
    sprintf("%s  |  %s of %s contigs kept  |  %d samples  |  assay: %s",
            r_filter$name() %||% "dataset",
            format(n_keep, big.mark = ","),
            format(nrow(raw), big.mark = ","),
            ncol(raw),
            r_filter$assay() %||% "none")
  })
}
