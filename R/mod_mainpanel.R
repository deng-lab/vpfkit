#' mainpanel UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_mainpanel_ui <- function(id) {
  ns <- NS(id)
  tagList(
    mainPanel(
      tabsetPanel(
        tabPanel("Data overview",
                 conditionalPanel("output.fileUploaded == false",
                                  fluidRow(column(12, span(textOutput("text1"), style="font-size: 1.2em;color:red;")))),
                 fluidRow(column(12, reactableOutput("tbl_rowdata"))),
                 fluidRow(column(6, plotlyOutput("plt_checkv_qc")),
                          column(6, plotlyOutput("plt_completeness")),
                          column(6, plotlyOutput("plt_adiversity")),
                          column(6, plotlyOutput("plt_bdiversity")),
                 ),
        ),
        tabPanel("Assay",
                 fluidRow(column(12, reactableOutput("tbl_abundance"))),
                 fluidRow(column(12, uiOutput("dl_button"))),
                 fluidRow(column(12, plotlyOutput("plt_abundance_barplot"))),
        ),
        tabPanel("Metadata", reactableOutput("tbl_smeta"),
                 fluidRow(column(12, reactableOutput("tbl_smeta2")))
        ),
        tabPanel("Gene annotation",
                 fluidRow(column(12, reactableOutput("tbl_geneanno"))),
                 fluidRow(column(6, plotlyOutput("plt_AMG")),
                          column(6, plotlyOutput("plt_PFAM"))),
                 fluidRow(column(6, plotlyOutput("plt_CARD")),
                          column(6, plotlyOutput("plt_VF"))),
        ),
        # == JBrowser UI from the JBrowser module ==
        # tabPanel("Genome Browser",
        #          fluidPage(JBrowserUI("jb")))
      ),
      width = 9
    )
  )
}

#' mainpanel Server Functions
#'
#' @noRd 
mod_mainpanel_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  })
}

## Copy in UI
# mod_mainpanel_ui("mainpanel_ui")

## Copy in server
# mod_mainpanel_server("mainpanel_ui")
