library(shiny)
library(shinythemes)
library(here)
library(tidyverse)
library(data.table)
library(reactable)
library(plotly)
library(scales)
library(speedyseq)
library(crosstalk)

# -------- load data --------
load(here("data/vpf.RData"))
checkv <- ctganno$checkv
catbat <- ctganno$catbat
vcontact <- ctganno$vcontact
taxonomy <- ctganno$taxonomy
bins <- ctganno$bins
# ctgmeta <- data.frame(tax_table(pseq)) %>% 
#   mutate(contig_length=as.integer(contig_length))
ctgmeta <- ctganno_merged
smeta <- data.frame(sample_data(pseq)) %>% 
  mutate(samplename=rownames(.))
abundance <- data.frame(otu_table(pseq))

minlen <- min(ctganno_merged$contig_length)
maxlen <- max(ctganno_merged$contig_length)


# -------- UI -------------
ui <- navbarPage("ViroProfiler Viewer",
            theme = shinytheme("cerulean"),
             tabPanel("Select",
                      titlePanel("Filter contigs"),
                      sidebarLayout(
                        sidebarPanel(
                          uiOutput("contig_len"),
                          uiOutput("checkv_quality"),
                          uiOutput("ctg_feature_col"),
                          uiOutput("sample_metadata"),
                          selectInput("abundance_scale", "Absolute/relative abundance",
                                      choices = c("Relative", "Absolute"),
                                      selected = "Absolute"),
                          numericInput("abundance_min_threshold", "Minimum abundance to show",
                                       value = 0.001, min=0),
                          width = 2
                        ),
                        
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Filter",
                                     fluidRow(column(6, plotlyOutput("coolplot")),
                                              column(6, plotlyOutput("taxa_abundance"))),
                                     reactableOutput("results"),
                            ),
                            tabPanel("Alpha diversity",
                                     plotOutput("alpha_diversity"),
                                     # plotlyOutput("ordination"),
                            ),), 
                          width=10))),
             tabPanel("EDA", "hello"),
             tabPanel("Contigs", "hello"),
             tabPanel("Genes", "hello"),
             tabPanel("PathAnno", "hello"),
)


# -------- Server ------------
server <- function(input, output) {
  
  # ------- render UI ---------
  output$contig_len <- renderUI({
    sliderInput("contig_len", "Contig length",
                min=minlen, 
                max=maxlen, 
                step = 500,
                ticks = FALSE,
                post = " bp",
                value=c(minlen,maxlen))
  })
  
  output$checkv_quality <- renderUI({
    selectInput("checkv_quality", "CheckV quality",
                sort(unique(ctgmeta$checkv_quality)),
                selected = c("Complete", "High-quality"),
                multiple = TRUE)
  })
  
  output$sample_metadata <- renderUI({
    selectInput("sample_metadata", "Sample metadata",
                colnames(smeta),
                selected = "samplename")
  })
  
  output$ctg_feature_col <- renderUI({
    selectInput("ctg_feature_col", "Feature column",
                colnames(ctgmeta),
                selected = "Family")
  })
  
  # -------- filter by checkv and length ---------
  ctgmeta_filtered <- reactive({
    if (is.null(input$checkv_quality)) {
      return(NULL)
    }
    ctgmeta %>% 
      filter(contig_length >= input$contig_len[1],
             contig_length <= input$contig_len[2],
             checkv_quality %in% input$checkv_quality)
  })
  
  ctgmeta_shared <- SharedData$new(ctgmeta_filtered, ~Contig, group="cmeta")

  # ------ subset Phyloseq -------
  pseq_subset <- reactive({
    filter_tax_table(pseq, Contig %in% ctgmeta_shared$key())
  })
  
  # ------ checkv quality --------
  output$coolplot <- renderPlotly({
    if (is.null(ctgmeta_filtered())) {
      return()
    }
    p <- ctgmeta_shared %>% 
      ggplot(aes(x=contig_length, color=checkv_quality, fill=checkv_quality)) +
      geom_histogram(alpha=0.5, position = "identity", bins=30) +
      scale_x_continuous(trans=log10_trans(),
                         breaks = trans_breaks("log10",
                                               function(x) as.integer(10^x))) +
      theme_bw()
    ggplotly(p)
  })
  
  # --------- taxonomy bar plot ---------
  output$taxa_abundance <- renderPlotly({
    # if (is.null(ctgmeta_filtered())) {
    #   return()
    # }
    taxa_collapse <-  tax_glom(pseq_subset(), taxrank = input$ctg_feature_col)
    if (input$abundance_scale=="Relative") {
      taxa_collapse <- taxa_collapse %>% 
        transform_sample_counts(function(x) {x/sum(x)})
    }
    p <- taxa_collapse %>% 
      psmelt() %>% 
      filter(Abundance > input$abundance_min_threshold) %>%
      ggplot(aes(x = Sample, y = Abundance, fill=.data[[input$ctg_feature_col]])) +    # Color by Phylum
      geom_bar(stat = "identity", position="stack") +
      guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
      ylab(paste0(input$abundance_scale, " abundance (> ", input$abundance_min_threshold, ")")) +
      theme_bw()

    # theme_bw()
    ggplotly(p)
  })
  
  # ------- feature annotation table --------
  output$results <- renderReactable({
    reactable(ctgmeta_shared,
              selection = "multiple", 
              filterable = TRUE,
              onClick = "select")
  })
  
  # -------- alpha diversity -----------
  output$alpha_diversity <- renderPlot({
    plot_richness(pseq_subset())
    # ggplotly(p)
  })
  
  # ---------- ordination ------
  # output$ordination <- renderPlot({
  #   pseq_ord <- ordinate(pseq_subset(), "PCoA", "bray")
  #   p <- plot_ordination(pseq_subset(), pseq_ord, type="samples") +
  #     geom_point(size=1, alpha=0.5) + 
  #     stat_ellipse(geom="polygon",level=0.95,alpha=0.2) +
  #     theme_bw()
  #   
  #   ggplotly(p)
  # })
}

shinyApp(ui = ui, server = server)
