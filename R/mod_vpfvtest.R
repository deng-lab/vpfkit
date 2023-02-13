#' vpfvtest UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_vpfvtest_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      fileInput(inputId = ns("file1"), label = "Upload viroprofiler_output.rds", multiple = FALSE, accept = c(".rds")),
      numericInput(ns("min_ctglen"), "Minimum contig length:", value = 10000, step = 1000),
      selectInput(ns("abdc_metric"), "Abundance metric:", choices = c("read counts", "trimmed mean"), selected = "read counts"),
      sliderInput(ns("min_covfrac"), "Min Coverage Fraction:", min = 0, max = 1, value = 0.5, step = 0.05, ticks = T),
      selectInput(ns("taxa_rank"), "Taxonomy rank", choices = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), selected = "Family"),
      uiOutput(ns("transform_on")),
      uiOutput(ns("transform_method")),
      uiOutput(ns("ctg_feature_col")),
      fluidRow(column(6, selectInput(ns("checkv_quality"), "CheckV quality", choices = c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"), selected = c("Complete", "High-quality", "Medium-quality", "Low-quality"), multiple = TRUE)),
               column(3, radioButtons(ns("filter_logic"), "Logic", choices = c("AND", "OR"), selected = "AND")),
               column(5, selectInput(ns("vs2_category"), "VirSorter2 category", choices = c(1,2,3,4,5,6), selected = c(1,2,4,5), multiple = TRUE))),
      numericInput(ns("abundance_min_threshold"), "Minimum abundance to show", value = 0.001, min=0),
      selectInput(ns("narm"), "Remove NA in glom:", choices = c(TRUE, FALSE), selected = FALSE),
      span("More information on ViroProfiler is "),
      a(href="https://github.com/deng-lab/viroprofiler", "on GitHub", target="_blank"),
      br(),
      span("Source code of this app is available "),
      a(href="https://github.com/deng-lab/vpfkit", "here", target="_blank"),
      width = 3
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Data overview",
                 conditionalPanel("output.fileUploaded == false",
                                  fluidRow(column(12, span(textOutput(ns("text1")), style="font-size: 1.2em;color:red;")))),
                 fluidRow(column(12, reactableOutput(ns("tbl_rowdata")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_checkv_qc"))),
                          column(6, plotlyOutput(ns("plt_completeness"))),
                          column(6, plotlyOutput(ns("plt_adiversity"))),
                          column(6, plotlyOutput(ns("plt_bdiversity"))),
                 ),
        ),
        # == JBrowser UI from the JBrowser module ==
        # tabPanel("Genome Browser",
        #          fluidPage(JBrowserUI("jb")))
      ),
      width = 9
    )
  )
}

#' vpfvtest Server Functions
#'
#' @noRd
mod_vpfvtest_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns

    tse_raw <- reactive({
      if (is.null(input$file1)) {
        return(NULL)
      } else {
        tse <- readRDS(input$file1$datapath)
        # print(input$checkv_quality)
        # print(tse)
        # print("hello world")
        # feature_anno <- SummarizedExperiment::rowData(tse) %>% SummarizedExperiment::as.data.frame()
        # tse_sel <- tse[1:5, 1]
        # print(tse_sel)
        # # print(feature_anno)
        # if (input$filter_logic == "AND") {
        #   tse_sel <- tse[feature_anno[["checkv_quality"]] %in% input$checkv_quality, ]
        # } else {
        #   tse_sel <- tse[feature_anno[["checkv_quality"]] %in% input$checkv_quality, ]
        # }
      }})

    output$fileUploaded <- reactive ({
      return(!is.null(tse_raw()))
    })

    outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

    output$text1 <- renderText({
      paste("Please upload a 'viroprofiler_output.rds' file.")
    })

    abdc_metric <- reactive({
      if (input$abdc_metric == "read counts") {
        return("count")
      } else if (input$abdc_metric == "trimmed mean") {
        return("trimmed_mean")
      }
    })

  })
}

## To be copied in the UI
# mod_vpfvtest_ui("vpfvtest_1")

## To be copied in the server
# mod_vpfvtest_server("vpfvtest_1")
