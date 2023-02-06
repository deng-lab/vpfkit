#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # Your application UI logic
    navbarPage("ViroProfiler-viewer",
               tabPanel("Home", mod_vpfilter_ui("vpfilter_1")),
               # tabPanel("Genome Browser", fluidPage(mod_jb2_ui("jb"))),
               # tabPanel("Tutorial", fixedPage(withMathJax(includeMarkdown("inst/app/www/tutorial.md")), hr(), div(class="footer", includeHTML("inst/app/www/footer.html")))),
               theme = shinythemes::shinytheme("flatly"))
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
      app_title = "vpfkit"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
