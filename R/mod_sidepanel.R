#' sidepanel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_sidepanel_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarPanel(
      fileInput(inputId = "file1", label = "Upload viroprofiler_output.RData", multiple = FALSE, accept = c(".RData")),
      numericInput("min_ctglen", "Minimum contig length:", value = 10000, step = 1000),
      selectInput("abdc_metric", "Abundance metric:", choices = c("read counts", "trimmed mean"), selected = "read counts"),
      sliderInput("min_covfrac", "Min Coverage Fraction:", min = 0, max = 1, value = 0.5, step = 0.05, ticks = T),
      selectInput("taxa_rank", "Taxonomy rank", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), selected = "Family"),
      uiOutput("transform_on"),
      uiOutput("transform_method"),
      uiOutput("ctg_feature_col"),
      fluidRow(column(6, selectInput("checkv_quality", "CheckV quality", choices = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"), selected = c("Complete", "High-quality", "Medium-quality", "Low-quality"), multiple = TRUE)),
               column(3, radioButtons("filter_logic", "Logic", choices = c("AND", "OR"), selected = "AND")),
               column(5, selectInput("vs2_category", "VirSorter2 category", choices = c(1,2,3,4,5,6), selected = c(1,2,4,5), multiple = TRUE))),
      numericInput("abundance_min_threshold", "Minimum abundance to show", value = 0.001, min=0),
      selectInput("narm", "Remove NA in glom:", choices = c(TRUE, FALSE), selected = FALSE),
      span("More information on ViroProfiler is "),
      a(href="https://github.com/deng-lab/viroprofiler", "on GitHub", target="_blank"),
      br(),
      span("Source code of this app is available "),
      a(href="https://github.com/deng-lab/viroprofiler-viewer", "here", target="_blank"),
      width = 3
    )
  )
}

#' sidepanel Server Functions
#'
#' @noRd 
mod_sidepanel_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    

  })
}

## Copy in UI
# mod_sidepanel_ui("sidepanel_ui")

## Copy in server
# mod_sidepanel_server("sidepanel_ui")
