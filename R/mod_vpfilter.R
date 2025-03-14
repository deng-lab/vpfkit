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
                 fluidRow(column(12, tags$em(textOutput(ns("filter_summary_text"))))),
                 fluidRow(column(12, reactableOutput(ns("tbl_rowdata")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_checkv_qc"))),
                          column(6, plotlyOutput(ns("plt_completeness"))),
                          column(6, plotlyOutput(ns("plt_adiversity"))),
                          column(6, plotlyOutput(ns("plt_bdiversity"))),
                 ),
                 fluidRow(column(12, downloadButton(ns("dl_annotations"), "Download annotations (TSV)"))),
                 fluidRow(column(12, plotlyOutput(ns("plt_heatmap")))),
        ),
        tabPanel("Assay",
                 fluidRow(column(12, reactableOutput(ns("tbl_abundance")))),
                 fluidRow(
                   column(4, downloadButton(ns("dl_abundance_csv"), "Download CSV")),
                   column(4, downloadButton(ns("dl_abundance_xlsx"), "Download Excel")),
                   column(4, downloadButton(ns("dl_tse_rds"), "Download TSE (.rds)"))
                 ),
                 fluidRow(column(12, plotlyOutput(ns("plt_abundance_barplot")))),
        ),
        tabPanel("Metadata",
                 reactableOutput(ns("tbl_smeta")),
                 fluidRow(column(12, downloadButton(ns("dl_metadata"), "Download metadata (CSV)")))
        ),
        tabPanel("Gene annotation",
                 fluidRow(column(12, reactableOutput(ns("tbl_geneanno")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_AMG"))),
                          column(6, plotlyOutput(ns("plt_PFAM")))),
                 fluidRow(column(6, plotlyOutput(ns("plt_CARD"))),
                          column(6, plotlyOutput(ns("plt_VF")))),
        ),
        tabPanel("Compare datasets",
                 fluidRow(
                   column(6, fileInput(inputId = ns("file2"), label = "Upload second dataset (.rds)", accept = c(".rds"))),
                   column(6, textOutput(ns("compare_status")))
                 ),
                 fluidRow(
                   column(6, h4("Dataset 1"), verbatimTextOutput(ns("compare_summary1"))),
                   column(6, h4("Dataset 2"), verbatimTextOutput(ns("compare_summary2")))
                 ),
                 fluidRow(column(12, plotlyOutput(ns("plt_compare_diversity")))),
                 fluidRow(column(12, plotlyOutput(ns("plt_compare_abundance")))),
                 fluidRow(column(12, plotlyOutput(ns("plt_shared_taxa"))))
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
            tryCatch({
                obj <- readRDS(input$file1$datapath)
                if (!inherits(obj, "TreeSummarizedExperiment")) {
                    shiny::showNotification("Uploaded file is not a TreeSummarizedExperiment object.",
                                           type = "error", duration = 10)
                    return(NULL)
                }
                obj
            }, error = function(e) {
                shiny::showNotification(paste("Error reading file:", e$message),
                                       type = "error", duration = 10)
                return(NULL)
            })
        }})

    output$fileUploaded <- reactive ({
        return(!is.null(tse_raw()))
    })

    outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

    output$text1 <- renderText({
        paste("Please upload a 'viroprofiler_output.rds' file.")
        })

    abdc_metric <- reactive({
      switch(input$abdc_metric,
        "read counts" = "counts",
        "trimmed mean" = "tmm",
        "counts"
      )
    })

    tse <- reactive({
      if (is.null(tse_raw())) return(NULL)
      tse_full <- tse_raw()
      rdata <- SummarizedExperiment::rowData(tse_full) %>% data.frame()

      rsel_len <- rdata$checkv_contig_length > input$min_ctglen
      rsel_quality <- rdata$checkv_quality %in% input$checkv_quality

      if (input$filter_logic == "AND") {
        rsel <- rsel_len & rsel_quality
      } else {
        rsel <- rsel_len | rsel_quality
      }

      tse <- tse_full[rsel, ]
      tse <- refind_abundance(tse, abdc_metric(), "covfrac", input$min_covfrac)

      if (nrow(tse) == 0) {
        shiny::showNotification("No contigs match current filter criteria.",
                               type = "warning", duration = 8)
      }
      tse
    })

    filter_summary <- reactive({
      if (is.null(tse_raw())) return("")
      sprintf("Showing %d of %d contigs (length > %s, quality: %s, logic: %s)",
              nrow(tse()), nrow(tse_raw()),
              format(input$min_ctglen, big.mark = ","),
              paste(input$checkv_quality, collapse = ", "),
              input$filter_logic)
    })

    feature_meta <- reactive({
      SummarizedExperiment::rowData(tse()) %>% data.frame()
    })

    output$filter_summary_text <- renderText({ filter_summary() })

    output$tbl_smeta <- renderReactable({
        if (is.null(tse_raw())) {
            return(NULL)
        } else {
          SummarizedExperiment::colData(tse()) %>%
            data.frame() %>%
            mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
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
            tse_div <- mia::addAlpha(tse(), assay.type = "counts",
                                     index = "shannon", name = "shannon")
            scater::plotColData(tse_div, "shannon", colour_by = "shannon")
        }
    })

    output$plt_bdiversity <- renderPlotly({
        if (is.null(tse_raw())) {
            return(NULL)
        } else if (ncol(tse()) < 2) {
            return(NULL)
        } else {
            plot_beta_diversity(tse(), name = "Beta-diversity", method = "bray",
                               assay.type = "counts", NMDS = TRUE)
        }
    })


    # ========== Abundance table ==========
    output$tbl_abundance <- renderReactable({
        if (is.null(tse_raw())) return(NULL)
        tse_obj <- tse_transformed()
        if (is.null(tse_obj) || nrow(tse_obj) == 0) return(reactable(data.frame()))
        assay_name <- active_assay()
        if (!assay_name %in% SummarizedExperiment::assayNames(tse_obj)) assay_name <- abdc_metric()
        SummarizedExperiment::assay(tse_obj, assay_name) %>%
              as.data.frame() %>%
              rownames_to_column("feature") %>%
              mutate(across(where(is.numeric), \(x) round(x, 4))) %>%
              reactable(
                fullWidth = FALSE, wrap = FALSE, resizable = TRUE,
                striped = TRUE, highlight = TRUE,
                columns = list(feature = colDef(width = 350)))
    })

    output$dl_annotations <- downloadHandler(
      filename = function() { paste0("annotations_", Sys.Date(), ".tsv") },
      content = function(file) { export_annotations(tse(), file, format = "tsv") }
    )

    output$dl_abundance_csv <- downloadHandler(
      filename = function() { paste0("abundance_", abdc_metric(), "_", Sys.Date(), ".csv") },
      content = function(file) {
        tse_obj <- tse_transformed()
        assay_name <- active_assay()
        if (!assay_name %in% SummarizedExperiment::assayNames(tse_obj)) assay_name <- abdc_metric()
        export_abundance(tse_obj, file, assay.type = assay_name, format = "csv")
      }
    )

    output$dl_abundance_xlsx <- downloadHandler(
      filename = function() { paste0("abundance_", abdc_metric(), "_", Sys.Date(), ".xlsx") },
      content = function(file) {
        tse_obj <- tse_transformed()
        assay_name <- active_assay()
        if (!assay_name %in% SummarizedExperiment::assayNames(tse_obj)) assay_name <- abdc_metric()
        export_abundance(tse_obj, file, assay.type = assay_name, format = "xlsx")
      }
    )

    output$dl_tse_rds <- downloadHandler(
      filename = function() { paste0("viroprofiler_filtered_", Sys.Date(), ".rds") },
      content = function(file) { export_vpftse(tse(), file) }
    )

    output$dl_metadata <- downloadHandler(
      filename = function() { paste0("metadata_", Sys.Date(), ".csv") },
      content = function(file) {
        df <- SummarizedExperiment::colData(tse()) %>% as.data.frame()
        utils::write.csv(df, file, row.names = TRUE)
      }
    )

    output$plt_heatmap <- renderPlotly({
      if (is.null(tse_raw()) || is.null(tse()) || nrow(tse()) == 0) return(NULL)
      tse_obj <- tse_transformed()
      mat <- SummarizedExperiment::assay(tse_obj, active_assay())
      top_n <- min(30, nrow(mat))
      top_idx <- order(rowMeans(mat), decreasing = TRUE)[1:top_n]
      mat_top <- mat[top_idx, , drop = FALSE]
      plotly::plot_ly(z = as.matrix(mat_top),
                      x = colnames(mat_top),
                      y = rownames(mat_top),
                      type = "heatmap",
                      colors = "YlOrRd") %>%
        plotly::layout(xaxis = list(title = "Sample"),
                       yaxis = list(title = "Contig"))
    })

    output$plt_abundance_barplot <- renderPlotly({
        if (is.null(tse_raw())) return(NULL)
        tse_obj <- tse_transformed()
        if (is.null(tse_obj) || nrow(tse_obj) == 0) return(NULL)
        assay_name <- active_assay()
        if (!assay_name %in% SummarizedExperiment::assayNames(tse_obj)) assay_name <- abdc_metric()
        miaViz::plotAbundance(tse_obj, assay.type = assay_name, rank = input$taxa_rank, use_relative = FALSE)
    })

    # Sidebar dynamic UI — transform controls
    output$transform_on <- renderUI({
      if (is.null(tse_raw())) return(NULL)
      checkboxInput(ns("do_transform"), "Transform abundance", value = FALSE)
    })

    output$transform_method <- renderUI({
      if (is.null(tse_raw()) || !isTRUE(input$do_transform)) return(NULL)
      selectInput(ns("transform_type"), "Method:",
                  choices = c("Relative abundance" = "relabundance",
                              "Log10(x+1)" = "log10p",
                              "CLR" = "clr"),
                  selected = "relabundance")
    })

    tse_transformed <- reactive({
      tse_obj <- tse()
      if (is.null(tse_obj) || nrow(tse_obj) == 0) return(tse_obj)
      if (!isTRUE(input$do_transform)) return(tse_obj)

      method <- input$transform_type
      if (is.null(method)) return(tse_obj)

      tryCatch({
        if (method == "relabundance") {
          mia::transformAssay(tse_obj, assay.type = abdc_metric(),
                             method = "relabundance", name = "transformed")
        } else if (method == "log10p") {
          mat <- SummarizedExperiment::assay(tse_obj, abdc_metric())
          SummarizedExperiment::assay(tse_obj, "transformed") <- log10(mat + 1)
          tse_obj
        } else if (method == "clr") {
          mia::transformAssay(tse_obj, assay.type = abdc_metric(),
                             method = "clr", pseudocount = 1, name = "transformed")
        } else {
          tse_obj
        }
      }, error = function(e) {
        shiny::showNotification(paste("Transform error:", e$message), type = "error")
        tse_obj
      })
    })

    active_assay <- reactive({
      if (isTRUE(input$do_transform) && !is.null(input$transform_type)) {
        "transformed"
      } else {
        abdc_metric()
      }
    })

    # Gene annotation data (filtered to current contigs)
    gene_anno <- reactive({
      tse_obj <- tse()
      if (is.null(tse_obj)) return(NULL)
      ga <- S4Vectors::metadata(tse_obj)$gene_annotations
      if (is.null(ga) || nrow(ga) == 0) return(NULL)
      # Filter to contigs that survived filtering
      current_contigs <- rownames(tse_obj)
      ga[ga$Contig %in% current_contigs, , drop = FALSE]
    })

    # Gene annotation table
    output$tbl_geneanno <- renderReactable({
      ga <- gene_anno()
      if (is.null(ga)) return(reactable(data.frame(message = "No gene annotations available. Upload a TSE with metadata(tse)$gene_annotations.")))
      reactable(ga, searchable = TRUE, striped = TRUE, highlight = TRUE,
                wrap = FALSE, resizable = TRUE, defaultPageSize = 15)
    })

    output$plt_AMG <- renderPlotly({
      .plot_gene_annotation_bar(gene_anno(), count_col = "dramv_ko",
        title = "Top AMG KO Functions", xlab = "KEGG Orthology", bar_color = "#FF6B6B",
        required_col = "dramv_category", category_col = "dramv_category", category_val = "AMG")
    })

    output$plt_PFAM <- renderPlotly({
      .plot_gene_annotation_bar(gene_anno(), count_col = "dramv_pfam",
        title = "Top 20 Pfam Domains", xlab = "Pfam domain", bar_color = "#4ECDC4")
    })

    output$plt_CARD <- renderPlotly({
      .plot_gene_annotation_bar(gene_anno(), count_col = "pharokka_card",
        title = "Antibiotic Resistance Genes (CARD)", xlab = "CARD ARO", bar_color = "#FFE66D",
        empty_message = "CARD: No antibiotic resistance genes detected")
    })

    output$plt_VF <- renderPlotly({
      .plot_gene_annotation_bar(gene_anno(), count_col = "pharokka_vfdb",
        title = "Virulence Factors (VFDB)", xlab = "VFDB ID", bar_color = "#A855F7",
        empty_message = "VFDB: No virulence factors detected")
    })

    # ========== Compare datasets tab ==========
    div_dataset1 <- reactive({
      req(tse())
      tryCatch({
        tse1 <- mia::addAlpha(tse(), assay.type = "counts", index = "shannon", name = "Shannon")
        cd1 <- as.data.frame(SummarizedExperiment::colData(tse1))
        cd1$Dataset <- "Dataset 1"
        cd1$Sample <- rownames(cd1)
        cd1[, c("Sample", "Shannon", "Dataset")]
      }, error = function(e) NULL)
    })

    div_dataset2 <- reactive({
      req(tse_compare())
      tryCatch({
        tse2 <- mia::addAlpha(tse_compare(), assay.type = "counts", index = "shannon", name = "Shannon")
        cd2 <- as.data.frame(SummarizedExperiment::colData(tse2))
        cd2$Dataset <- "Dataset 2"
        cd2$Sample <- rownames(cd2)
        cd2[, c("Sample", "Shannon", "Dataset")]
      }, error = function(e) NULL)
    })

    tse_compare <- reactive({
      if (is.null(input$file2)) return(NULL)
      tryCatch({
        obj <- readRDS(input$file2$datapath)
        if (!inherits(obj, "TreeSummarizedExperiment")) {
          shiny::showNotification("Second file is not a TSE object.", type = "error")
          return(NULL)
        }
        obj
      }, error = function(e) {
        shiny::showNotification(paste("Error reading second file:", e$message), type = "error")
        return(NULL)
      })
    })

    output$compare_status <- renderText({
      if (is.null(tse_compare())) return("Upload a second dataset to compare.")
      "Second dataset loaded."
    })

    output$compare_summary1 <- renderText({
      if (is.null(tse())) return("")
      sprintf("Contigs: %d\nSamples: %d\nAssays: %s",
              nrow(tse()), ncol(tse()),
              paste(SummarizedExperiment::assayNames(tse()), collapse = ", "))
    })

    output$compare_summary2 <- renderText({
      if (is.null(tse_compare())) return("")
      sprintf("Contigs: %d\nSamples: %d\nAssays: %s",
              nrow(tse_compare()), ncol(tse_compare()),
              paste(SummarizedExperiment::assayNames(tse_compare()), collapse = ", "))
    })

    output$plt_compare_diversity <- renderPlotly({
      if (is.null(tse()) || is.null(tse_compare())) return(NULL)
      div1 <- div_dataset1()
      div2 <- div_dataset2()
      if (is.null(div1) && is.null(div2)) return(NULL)
      div_all <- rbind(div1, div2)
      p <- ggplot(div_all, aes(x = Sample, y = Shannon, color = Dataset)) +
        geom_point(size = 3) +
        theme_minimal() +
        labs(title = "Shannon Diversity Comparison") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p)
    })

    output$plt_compare_abundance <- renderPlotly({
      if (is.null(tse()) || is.null(tse_compare())) return(NULL)
      # Compare total abundance per sample
      counts1 <- colSums(SummarizedExperiment::assay(tse(), "counts"))
      counts2 <- colSums(SummarizedExperiment::assay(tse_compare(), "counts"))
      df <- data.frame(
        Sample = c(names(counts1), names(counts2)),
        Total_reads = c(counts1, counts2),
        Dataset = c(rep("Dataset 1", length(counts1)), rep("Dataset 2", length(counts2)))
      )
      p <- ggplot(df, aes(x = Sample, y = Total_reads, fill = Dataset)) +
        geom_col(position = "dodge") +
        theme_minimal() +
        labs(title = "Total Read Counts per Sample", y = "Total reads") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      plotly::ggplotly(p)
    })

    output$plt_shared_taxa <- renderPlotly({
      if (is.null(tse()) || is.null(tse_compare())) return(NULL)
      rd1 <- as.data.frame(SummarizedExperiment::rowData(tse()))
      rd2 <- as.data.frame(SummarizedExperiment::rowData(tse_compare()))
      if (!"Family" %in% colnames(rd1) || !"Family" %in% colnames(rd2)) return(NULL)

      fam1 <- unique(stats::na.omit(rd1$Family))
      fam2 <- unique(stats::na.omit(rd2$Family))
      shared <- intersect(fam1, fam2)
      only1 <- setdiff(fam1, fam2)
      only2 <- setdiff(fam2, fam1)

      df <- data.frame(
        Category = c("Shared", "Only in Dataset 1", "Only in Dataset 2"),
        Count = c(length(shared), length(only1), length(only2))
      )
      plotly::plot_ly(df, x = ~Category, y = ~Count, type = "bar",
                      marker = list(color = c("#4ECDC4", "#FF6B6B", "#A855F7"))) %>%
        plotly::layout(title = "Shared vs Unique Viral Families",
                       xaxis = list(title = ""),
                       yaxis = list(title = "Number of families"))
    })

  })
}

#' Bar chart helper for gene annotation plots
#'
#' @noRd
.plot_gene_annotation_bar <- function(
    ga, count_col, title, xlab, bar_color,
    required_col = count_col,
    category_col = NULL, category_val = NULL,
    empty_message = NULL
) {
  if (is.null(ga) || !required_col %in% colnames(ga)) return(NULL)
  df <- if (!is.null(category_col) && !is.null(category_val)) {
    ga[!is.na(ga[[category_col]]) & ga[[category_col]] == category_val & !is.na(ga[[count_col]]), ]
  } else {
    ga[!is.na(ga[[count_col]]), ]
  }
  if (nrow(df) == 0) {
    if (!is.null(empty_message)) {
      return(plotly::plot_ly() %>%
        plotly::layout(title = empty_message,
                      annotations = list(text = empty_message, showarrow = FALSE,
                                        font = list(size = 14))))
    }
    return(NULL)
  }
  counts_df <- as.data.frame(table(df[[count_col]]), stringsAsFactors = FALSE)
  colnames(counts_df) <- c("Label", "Count")
  counts_df <- counts_df[order(-counts_df$Count), ]
  counts_df <- head(counts_df, 20)
  plotly::plot_ly(counts_df, x = ~reorder(Label, Count), y = ~Count,
                  type = "bar", marker = list(color = bar_color)) %>%
    plotly::layout(title = title,
                   xaxis = list(title = xlab, tickangle = -45),
                   yaxis = list(title = "Gene count"))
}

## To be copied in the UI
# mod_vpfilter_ui("vpfilter_1")

## To be copied in the server
# mod_vpfilter_server("vpfilter_1")
