library(shiny)
library(here)
library(shinyWidgets)
library(tidyverse)
library(reactable)
library(plotly)
library(shinythemes)
library(miaViz)
library(mia)
library(scater)
library(markdown)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
source(here("R/utils.R"))


# Functions
plot_beta_diversity <- function(tse_obj, name, NMDS=FALSE, scale = F, ...) {
  if (NMDS == TRUE) {
    tse_obj <- runNMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    xlab <- "Axis 1"
    ylab <- "Axis 2"
  } else {
    tse_obj <- runMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    e <- attr(reducedDim(tse_obj, name), "eig")
    rel_eig <- 100 * e/sum(e[e>0])
    xlab <- paste("Axis 1 (", round(rel_eig[[1]], 2), "%)", sep = "")
    ylab <- paste("Axis 2 (", round(rel_eig[[2]], 2), "%)", sep = "")
  }
  p <- plotReducedDim(tse_obj, name, colour_by = "condition") +
    xlab(xlab) +
    ylab(ylab) +
    stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=colour_by), linetype=2) +
    ggtitle(name) +
    theme_bw()
  return(p)
}


# ================ JBrowser Module ========================================

library(JBrowseR)
data_server <- serve_data("./JBrowser_data/")

JBrowserUI <- function(id) {
  tagList(
    # this adds to the browser to the UI, and specifies the output ID in the server
    JBrowseROutput(NS(id,"browserOutput"))
  )
}

JBrowserServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # create the necessary JB2 assembly configuration
    assembly <- assembly(
      "http://127.0.0.1:5000/phanotate_nts.fasta.gz",
      bgzip = TRUE
    )
    # create configuration for a JB2 GFF Feature Track
    annotations_track <- track_feature(
      # "https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz",
      "http://127.0.0.1:5000/trnascan_out.gff.gz",
      assembly
    )
    annotations_track2 <- track_feature(
      # "https://jbrowse.org/genomes/sars-cov2/sars-cov2-annotations.sorted.gff.gz",
      "http://127.0.0.1:5000/vpf.genbank.sorted.gff.gz",
      assembly
    )
    # create the tracks array to pass to browser
    tracks <- tracks(annotations_track, annotations_track2)
    theme <- theme("#333", "#ff6200")
    default_session <- default_session(
      assembly,
      annotations_track2
    )
    # link the UI with the browser widget
    output$browserOutput <- renderJBrowseR(
      JBrowseR(
        "View",
        assembly = assembly,
        tracks = tracks,
        theme = theme
        # defaultSession = default_session
      )
    )
  })
}




# ========= sidebar ==========
sideP <- sidebarPanel(
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

# ========= main panel ==========
mainP <- mainPanel(
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

navp_raw <- tabPanel("Overview", sidebarLayout(sideP, mainP))
navp_tutorial <- tabPanel("Tutorial", fixedPage(withMathJax(includeMarkdown("www/tutorial.md")), hr(), div(class="footer", includeHTML("www/footer.html"))))
# == JBrowser UI from the JBrowser module ==
navp_JBrowser <- tabPanel("Genome Browser", fluidPage(JBrowserUI("jb")))
ui <- navbarPage("ViroProfiler-viewer", navp_raw, navp_JBrowser, navp_tutorial, theme = shinytheme("united"))

# load("viroprofiler_output.RData")
# ========= server =================
server <- function(input, output) {
    tse_raw <- reactive({
        if (is.null(input$file1)) {
            return(NULL)
        } else {
            load(input$file1$datapath)
            feature_anno <- rowData(tse)
            if (input$filter_logic == "AND") {
                tse_sel <- tse[feature_anno$checkv_quality %in% input$checkv_quality &
                               feature_anno$vs2_category %in% input$vs2_category, ]
            } else {
                tse_sel <- tse[feature_anno$checkv_quality %in% input$checkv_quality |
                               feature_anno$vs2_category %in% input$vs2_category, ]
            }
        }})

    output$fileUploaded <- reactive ({
        return(!is.null(tse_raw()))
    })

    outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

    output$text1 <- renderText({
        paste("Please upload a 'viroprofiler_output.RData' file.")
        })

    abdc_metric <- reactive({
        if (input$abdc_metric == "read counts") {
            return("counts")
        } else if (input$abdc_metric == "trimmed mean") {
            return("trimmed_mean")
        }
    })

    tse <- reactive({
        tse <- tse_raw()
        tse <- tse[rowData(tse)$checkv_contig_length>input$min_ctglen,]
        tse <- refind_abundance(tse, abdc_metric(), "covfrac", input$min_covfrac)
        # tse <- tse[rowData()]
    })


    sample_meta <- reactive({
        colData(tse()) %>% 
            as.data.frame() %>%
            mutate(across(where(is.numeric), round, 4))
    })

    feature_meta <- reactive({
        rowData(tse()) %>% as.data.frame()
    })

    output$tbl_smeta <- renderReactable({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            sample_meta() %>%
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
                geom_histogram(alpha=0.8)
        }
    })

    output$plt_completeness <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            feature_meta() %>%
                ggplot(aes(x=checkv_completeness, fill=checkv_quality)) +
                geom_histogram(alpha=0.8)
        }
    })

    output$plt_adiversity <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
            plotColData(tse(), "richness_observed", "condition", colour_by = "condition")
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
              assay(tse(), abdc_metric()) %>%
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
            plotAbundance(tse(), abund_values=abdc_metric(), rank = input$taxa_rank, use_relative=F)
        }
    })
    
    # == JBrowser server from the JBrowser module ==
    JBrowserServer("jb")
    
}


# Run the application 
shinyApp(ui = ui, server = server)
