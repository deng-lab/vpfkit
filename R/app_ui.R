#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    navbarPage(
      title = "ViroProfiler-viewer",
      id = "main_nav",
      theme = shinythemes::shinytheme("flatly"),
      header = div(
        class = "vpf-statusbar",
        textOutput("status_bar", inline = TRUE)
      ),
      tabPanel("Data", mod_data_input_ui("data")),
      tabPanel("Filter", mod_vpfilter_ui("filter")),
      tabPanel("Taxonomy", mod_composition_ui("composition")),
      tabPanel("Diversity", mod_diversity_ui("diversity")),
      tabPanel("Contigs", mod_features_ui("features")),
      tabPanel("Host & lifestyle", mod_host_ui("host")),
      tabPanel("Genes", mod_genes_ui("genes")),
      tabPanel("Compare", mod_compare_ui("compare")),
      tabPanel("Export", mod_export_ui("export")),
      tabPanel("About", mod_about_ui("about"))
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "ViroProfiler-viewer"
    )
  )
}
