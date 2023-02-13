#' vpfilter UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_vpfilter_ui <- function(id){
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
        tabPanel("Assay",
                 fluidRow(column(12, reactableOutput(ns("tbl_abundance")))),
                 fluidRow(column(12, uiOutput(ns("dl_button")))),
                 fluidRow(column(12, plotlyOutput(ns("plt_abundance_barplot")))),
        ),
        tabPanel("Metadata", reactableOutput(ns("tbl_smeta")),
                 fluidRow(column(12, reactableOutput(ns("tbl_smeta2"))))
        ),
        tabPanel("Gene annotation",
                 fluidRow(column(12, reactableOutput(ns("tbl_geneanno")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_AMG"))),
                          column(6, plotlyOutput(ns("plt_PFAM")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_CARD"))),
                          column(6, plotlyOutput(ns("plt_VF")))),
        ),
        # == JBrowser UI from the JBrowser module ==
        # tabPanel("Genome Browser",
        #          fluidPage(JBrowserUI("jb")))
      ),
      width = 9
    )
  )
}

#' vpfilter Server Functions
#'
#' @noRd
mod_vpfilter_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    tse_raw <- reactive({
        if (is.null(input$file1)) {
            return(NULL)
        } else {
            tse <- readRDS(input$file1$datapath)
            # print(tse)
            # feature_anno <- SummarizedExperiment::rowData(tse) %>% data.frame()
            # # print(feature_anno)
            # if (input$filter_logic == "AND") {
            #     tse_sel <- tse[feature_anno[["checkv_quality"]] %in% input$checkv_quality, ]
            # } else {
            #     tse_sel <- tse[feature_anno[["checkv_quality"]] %in% input$checkv_quality, ]
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

    tse <- reactive({
      if (is.null(tse_raw())) {
        return(NULL)
      } else {
        tse <- tse_raw()
        rdata <- SummarizedExperiment::rowData(tse) %>% data.frame()
        rsel <- rdata$checkv_contig_length>input$min_ctglen
        tse <- tse[rsel,]
        tse <- refind_abundance(tse, abdc_metric(), "covfrac", input$min_covfrac)
        # tse <- tse[rowData()]
      }
    })

    feature_meta <- reactive({
      SummarizedExperiment::rowData(tse()) %>% data.frame()
    })

    output$tbl_smeta <- renderReactable({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
          SummarizedExperiment::colData(tse()) %>%
            data.frame() %>%
            mutate(across(where(is.numeric), round, 4)) %>%
              reactable(
                fullWidth = FALSE,
                wrap = FALSE,
                resizable = TRUE,
                striped = TRUE,
                highlight = TRUE,
                bordered = TRUE,
              )
        }
    })

    # ========= Data overview ==========
    output$tbl_rowdata <- renderReactable({
      if (is.null(tse_raw())) {
        return(NULL)
      } else {
        feature_meta() %>%
          reactable(wrap = F, resizable = T)
      }
    })

    output$plt_checkv_qc <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            feature_meta() %>%
                # dplyr::filter(!checkv_quality %in% c("Not-determined", "Low-quality")) %>%
                ggplot(aes(x=checkv_contig_length, fill=checkv_quality)) +
                ggplot2::geom_histogram(alpha=0.8)
        }
    })

    output$plt_completeness <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            feature_meta() %>%
              ggplot(aes(x=checkv_completeness, fill=checkv_quality)) +
              ggplot2::geom_histogram(alpha=0.8)
        }
    })

    output$plt_adiversity <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            scater::plotColData(tse(), "richness_observed", "condition", colour_by = "condition")
        }
    })

    output$plt_bdiversity <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            plot_beta_diversity(tse(), name = "Beta-diversity", method = "bray", exprs_values = "relabundance", ntop=nrow(tse()), NMDS = T)
        }
    })


    # ========== Abundance table ==========
    output$tbl_abundance <- renderReactable({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
          SummarizedExperiment::assay(tse(), abdc_metric()) %>%
                as.data.frame() %>%
                rownames_to_column("feature") %>%
                mutate(across(where(is.numeric), round, 4)) %>%
                reactable(
                  fullWidth = FALSE,
                  wrap = FALSE,
                  resizable = TRUE,
                  striped = TRUE,
                  highlight = TRUE,
                  elementId = "table_abundance",
                  columns = list(feature = colDef(width = 350)))
        }
    })

    output$dl_button <- renderUI({
      if (is.null(tse_raw())) {
        return(NULL)
      } else {
      htmltools::browsable(
        tagList(
          tags$button(
            tagList(fontawesome::fa("download"), "Download as CSV"),
            onclick = "Reactable.downloadDataCSV('table_abundance', 'abundance.csv')"
          )))
      }
    })

    output$plt_abundance_barplot <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            miaViz::plotAbundance(tse(), abund_values=abdc_metric(), rank = input$taxa_rank, use_relative=F)
        }
    })

  })
}

## To be copied in the UI
# mod_vpfilter_ui("vpfilter_1")

## To be copied in the server
# mod_vpfilter_server("vpfilter_1")
