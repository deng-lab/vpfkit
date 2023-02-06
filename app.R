library("shiny")
library("here")
library("shinyWidgets")
library("tidyverse")
library("reactable")
library("plotly")
library("shinythemes")
library("miaViz")
library("mia")
library("scater")
library("markdown")
library("JBrowseR")

# Load packages
packages <- c("shiny", "here", "shinyWidgets", "tidyverse", "reactable", "plotly", "shinythemes", "miaViz", "mia", "scater", "markdown", "JBrowseR")
for (package in packages) {
  if (!require(package, character.only = TRUE)) {
    renv::install(package)
  }
}
lapply(packages, library, character.only = TRUE)


Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)
source(here("R/utils.R"))
source(here("R/mod_jbrowse.R"))
source(here("R/mod_sidepanel.R"))
source(here("R/mod_mainpanel.R"))

navp_raw <- tabPanel("Overview", sidebarLayout(mod_sidepanel_ui("sidepanel_ui"), mod_mainpanel_ui("mainpanel_ui")))
navp_tutorial <- tabPanel("Tutorial", fixedPage(withMathJax(includeMarkdown("www/tutorial.md")), hr(), div(class="footer", includeHTML("www/footer.html"))))
navp_JBrowser <- tabPanel("Genome Browser", fluidPage(mod_jb2_ui("jb")))
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
    mod_jb2_server("jb")

}


# Run the application
shinyApp(ui = ui, server = server)
