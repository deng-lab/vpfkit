#' features UI Function
#'
#' @description Searchable browser over every per-contig annotation column,
#'   with a drill-down that shows one contig's annotations and its abundance
#'   across samples.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_features_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("status")),
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Annotation table",
        vpf_caption(
          "All annotation columns for the contigs that pass the current filters.",
          "Use the column selector to narrow the view; the download always",
          "contains the full search result, not just the visible page."
        ),
        checkboxInput(ns("include_abundance"),
                      "Append per-sample abundance columns", value = FALSE),
        mod_bigtable_ui(ns("tbl"), height = "560px")
      ),
      tabPanel(
        "Single contig",
        fluidRow(
          column(6,
                 # Server-side selectize: an annotation table can hold a hundred
                 # thousand contigs, and a plain selectInput would ship every
                 # identifier into the page.
                 selectizeInput(ns("contig"), "Contig", choices = NULL,
                                width = "100%",
                                options = list(placeholder = "Type to search",
                                               maxOptions = 200)),
                 uiOutput(ns("contig_hint"))),
          column(6, uiOutput(ns("contig_summary")))
        ),
        plotly::plotlyOutput(ns("plt_contig"), height = "320px"),
        vpf_caption("Abundance of the selected contig in each sample, after coverage masking."),
        reactable::reactableOutput(ns("tbl_contig"))
      ),
      tabPanel(
        "Column overview",
        vpf_caption(
          "What each annotation column contains for the current selection.",
          "Columns that are entirely missing are the fastest way to see which",
          "pipeline stages did not run."
        ),
        reactable::reactableOutput(ns("tbl_columns"))
      )
    )
  )
}

#' features Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_features_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_notice(type = "warning",
                          "No contigs pass the current filters. Loosen them on the Filter tab."))
      }
      NULL
    })

    annotation_table <- reactive({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(NULL)
      }
      df <- vpf_row_table(tse)
      if (isTRUE(input$include_abundance)) {
        assay_name <- r_filter$assay()
        if (!is.null(assay_name) && assay_name %in% SummarizedExperiment::assayNames(tse)) {
          mat <- as.data.frame(SummarizedExperiment::assay(tse, assay_name))
          colnames(mat) <- paste0(vpf_sample_names(tse), " [", assay_name, "]")
          df <- cbind(df, mat)
        }
      }
      df
    })

    default_cols <- reactive({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(NULL)
      }
      preferred <- c("Contig", "checkv_contig_length", "checkv_quality",
                     "checkv_completeness", "genomad_score", "virsorter2_max_score",
                     VPF_TAXONOMY_RANKS, "iphop_genus", "bacphlip_replicyc")
      cols <- intersect(preferred, colnames(annotation_table()))
      if (length(cols) < 4) NULL else cols
    })

    mod_bigtable_server("tbl", r_data = annotation_table,
                        filename = "contig_annotations",
                        default_columns = default_cols)

    # -- Single-contig drill-down --------------------------------------------
    contig_ids <- reactive({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(character(0))
      }
      rownames(tse) %||% character(0)
    })

    observeEvent(contig_ids(), {
      ids <- contig_ids()
      updateSelectizeInput(session, "contig", choices = ids,
                           selected = if (length(ids)) ids[[1]] else character(0),
                           server = TRUE)
    }, ignoreNULL = FALSE)

    output$contig_hint <- renderUI({
      ids <- contig_ids()
      if (length(ids) == 0) {
        return(vpf_caption("No contig passes the current filters."))
      }
      vpf_caption(sprintf("%s contigs available. Only matching entries are sent to the browser.",
                          format(length(ids), big.mark = ",")))
    })

    selected_contig <- reactive({
      ids <- contig_ids()
      if (length(ids) == 0) {
        return(NULL)
      }
      sel <- input$contig
      if (is.null(sel) || !nzchar(sel) || !sel %in% ids) {
        return(ids[[1]])
      }
      sel
    })

    output$contig_summary <- renderUI({
      tse <- r_filter$tse()
      cid <- selected_contig()
      if (is.null(tse) || is.null(cid)) {
        return(NULL)
      }
      rd <- SummarizedExperiment::rowData(tse)[cid, , drop = FALSE]
      get <- function(nm) {
        if (!nm %in% colnames(rd)) {
          return("-")
        }
        v <- vpf_blank_na(rd[[nm]])[[1]]
        if (is.na(v)) "-" else v
      }
      div(
        vpf_stat("Length (bp)", get("checkv_contig_length")),
        vpf_stat("CheckV quality", get("checkv_quality")),
        vpf_stat("Completeness", get("checkv_completeness")),
        vpf_stat("Family", get("Family")),
        vpf_stat("Host genus", get("iphop_genus")),
        vpf_stat("Lifestyle", get("bacphlip_replicyc"))
      )
    })

    output$plt_contig <- plotly::renderPlotly({
      tse <- r_filter$tse()
      cid <- selected_contig()
      assay_name <- r_filter$assay()
      if (is.null(tse) || is.null(cid) || is.null(assay_name)) {
        return(vpf_message_plot("Select a contig first."))
      }
      if (!assay_name %in% SummarizedExperiment::assayNames(tse)) {
        return(vpf_message_plot("The active assay is not present in this object."))
      }
      vals <- as.numeric(SummarizedExperiment::assay(tse, assay_name)[cid, ])
      samples <- vpf_sample_names(tse)
      df <- data.frame(Sample = factor(samples, levels = samples),
                       Value = vals, stringsAsFactors = FALSE)
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Sample, y = .data$Value,
        text = paste0(.data$Sample, "<br>", signif(.data$Value, 4))
      )) +
        ggplot2::geom_col(fill = "#F28E2B") +
        ggplot2::labs(x = NULL, y = assay_name, title = cid) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      vpf_ggplotly(p, tooltip = "text")
    })

    output$tbl_contig <- reactable::renderReactable({
      tse <- r_filter$tse()
      cid <- selected_contig()
      if (is.null(tse) || is.null(cid)) {
        return(reactable::reactable(data.frame(Message = "Select a contig first.")))
      }
      rd <- as.data.frame(SummarizedExperiment::rowData(tse)[cid, , drop = FALSE])
      df <- data.frame(
        Field = colnames(rd),
        Value = vapply(colnames(rd), function(nm) {
          v <- rd[[nm]][[1]]
          if (is.null(v) || is.na(v)) "-" else as.character(v)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      df <- df[df$Value != "-", , drop = FALSE]
      reactable::reactable(df, striped = TRUE, bordered = TRUE, highlight = TRUE,
                           defaultPageSize = 20, searchable = TRUE)
    })

    output$tbl_columns <- reactable::renderReactable({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(reactable::reactable(data.frame(Message = "No contigs in the current selection.")))
      }
      rd <- SummarizedExperiment::rowData(tse)
      rows <- lapply(colnames(rd), function(nm) {
        v <- rd[[nm]]
        clean <- vpf_blank_na(v)
        n_ok <- sum(!is.na(clean))
        uniq <- unique(clean[!is.na(clean)])
        data.frame(
          Column = nm,
          Type = class(v)[[1]],
          Filled = n_ok,
          `Percent filled` = round(100 * n_ok / nrow(tse), 1),
          `Distinct values` = length(uniq),
          Example = if (length(uniq) == 0) "-" else paste(utils::head(uniq, 3), collapse = " | "),
          check.names = FALSE, stringsAsFactors = FALSE
        )
      })
      df <- do.call(rbind, rows)
      reactable::reactable(df, striped = TRUE, bordered = TRUE, highlight = TRUE,
                           searchable = TRUE, defaultPageSize = 20,
                           columns = list(Example = reactable::colDef(minWidth = 300)))
    })
  })
}

## To be copied in the UI
# mod_features_ui("features_1")

## To be copied in the server
# mod_features_server("features_1", r_filter)
