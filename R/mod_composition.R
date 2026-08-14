#' composition UI Function
#'
#' @description Taxonomic composition of the filtered contig set. Contigs
#'   without an assignment at the chosen rank are kept and labelled
#'   `Unclassified`: in a virome they are usually the largest and the most
#'   interesting fraction, and dropping them silently rescales every other bar.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_composition_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        3,
        vpf_card(
          "Display",
          uiOutput(ns("rank_ui")),
          radioButtons(ns("scale"), "Values",
                       choices = c("Relative abundance" = "relative",
                                   "Absolute abundance" = "absolute"),
                       selected = "relative"),
          numericInput(ns("top_n"), "Number of taxa shown separately",
                       value = 12, min = 2, max = 40, step = 1),
          vpf_caption("Remaining taxa are pooled into 'Other'."),
          uiOutput(ns("group_ui")),
          checkboxInput(ns("drop_unclassified"),
                        "Hide unclassified contigs", value = FALSE),
          uiOutput(ns("unclassified_warning"))
        )
      ),
      column(
        9,
        uiOutput(ns("status")),
        tabsetPanel(
          id = ns("tabs"),
          tabPanel(
            "Stacked composition",
            plotly::plotlyOutput(ns("plt_stack"), height = "480px"),
            vpf_caption(textOutput(ns("stack_caption"), inline = TRUE))
          ),
          tabPanel(
            "Classification depth",
            plotly::plotlyOutput(ns("plt_ranks"), height = "380px"),
            vpf_caption(
              "How many of the selected contigs carry an assignment at each rank.",
              "Viral taxonomy is sparse by nature; a steep drop toward Genus and",
              "Species is expected rather than a defect."
            )
          ),
          tabPanel(
            "Taxon table",
            mod_bigtable_ui(ns("tbl"), height = "460px")
          )
        )
      )
    )
  )
}

#' composition Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_composition_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ranks <- reactive(vpf_available_ranks(r_filter$tse()))

    output$rank_ui <- renderUI({
      rk <- ranks()
      if (length(rk) == 0) {
        return(vpf_caption("No taxonomic rank carries data in this selection."))
      }
      default <- if ("Family" %in% rk) "Family" else rk[[length(rk)]]
      selectInput(ns("rank"), "Taxonomic rank", choices = rk, selected = default)
    })

    output$group_ui <- renderUI({
      groups <- vpf_group_candidates(r_filter$tse())
      if (length(groups) == 0) {
        return(vpf_caption(paste(
          "No grouping variable in colData, so samples are shown individually.",
          "Upload a sample-metadata table on the Data tab to group them."
        )))
      }
      selectInput(ns("group"), "Facet samples by",
                  choices = c("(none)" = "", groups), selected = "")
    })

    active_rank <- reactive({
      rk <- ranks()
      if (length(rk) == 0) {
        return(NULL)
      }
      sel <- input$rank
      if (is.null(sel) || !sel %in% rk) rk[[min(length(rk), 6)]] else sel
    })

    # Taxon x sample matrix, computed once and reused by every output.
    composition <- reactive({
      tse <- r_filter$tse()
      assay_name <- r_filter$assay()
      rank <- active_rank()
      if (is.null(tse) || is.null(assay_name) || is.null(rank)) {
        return(NULL)
      }
      if (nrow(tse) == 0 || ncol(tse) == 0) {
        return(NULL)
      }
      if (!assay_name %in% SummarizedExperiment::assayNames(tse)) {
        return(NULL)
      }
      labels <- vpf_blank_na(SummarizedExperiment::rowData(tse)[[rank]])
      labels[is.na(labels)] <- "Unclassified"
      mat <- SummarizedExperiment::assay(tse, assay_name)
      mat[is.na(mat)] <- 0
      agg <- rowsum(as.matrix(mat), group = labels, reorder = TRUE)
      n_unclassified <- sum(labels == "Unclassified")
      if (isTRUE(input$drop_unclassified)) {
        agg <- agg[rownames(agg) != "Unclassified", , drop = FALSE]
      }
      list(matrix = agg, rank = rank, assay = assay_name,
           n_unclassified = n_unclassified, n_features = nrow(tse))
    })

    output$unclassified_warning <- renderUI({
      cmp <- composition()
      if (is.null(cmp) || !isTRUE(input$drop_unclassified)) {
        return(NULL)
      }
      vpf_notice(
        type = "warning",
        sprintf(paste("%d of %d contigs have no assignment at rank %s and are",
                      "currently hidden. The remaining bars are rescaled, so the",
                      "proportions no longer describe the whole community."),
                cmp$n_unclassified, cmp$n_features, cmp$rank)
      )
    })

    output$status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_notice(type = "warning",
                          "No contigs pass the current filters. Loosen them on the Filter tab."))
      }
      if (length(ranks()) == 0) {
        return(vpf_missing_annotation(
          "Taxonomy",
          "--use_vitap plus vConTACT3, merged by bin/merge_taxonomy.py",
          "No taxonomic rank in this object carries an assignment, so composition cannot be shown."
        ))
      }
      NULL
    })

    plot_data <- reactive({
      cmp <- composition()
      if (is.null(cmp) || nrow(cmp$matrix) == 0) {
        return(NULL)
      }
      mat <- cmp$matrix
      if (identical(input$scale, "relative")) {
        totals <- colSums(mat)
        totals[totals == 0] <- NA_real_
        mat <- sweep(mat, 2, totals, "/")
        mat[is.na(mat)] <- 0
      }
      top_n <- suppressWarnings(as.integer(input$top_n))
      if (is.na(top_n) || top_n < 2) top_n <- 12L
      if (nrow(mat) > top_n) {
        ord <- order(rowSums(mat), decreasing = TRUE)
        keep <- rownames(mat)[ord[seq_len(top_n)]]
        other <- colSums(mat[setdiff(rownames(mat), keep), , drop = FALSE])
        mat <- rbind(mat[keep, , drop = FALSE], Other = other)
      }
      samples <- vpf_sample_names(r_filter$tse())
      df <- data.frame(
        Taxon = rep(rownames(mat), times = ncol(mat)),
        Sample = rep(samples, each = nrow(mat)),
        Value = as.numeric(mat),
        stringsAsFactors = FALSE
      )
      # Stack the most abundant taxa at the bottom; keep Other and
      # Unclassified visually separate at the top.
      ordering <- rownames(mat)
      special <- intersect(c("Other", "Unclassified"), ordering)
      ordering <- c(setdiff(ordering, special), special)
      df$Taxon <- factor(df$Taxon, levels = rev(ordering))
      df$Sample <- factor(df$Sample, levels = samples)

      grp <- input$group
      if (!is.null(grp) && nzchar(grp)) {
        cd <- SummarizedExperiment::colData(r_filter$tse())
        if (grp %in% colnames(cd)) {
          lut <- stats::setNames(as.character(vpf_blank_na(cd[[grp]])), samples)
          df$Group <- lut[as.character(df$Sample)]
          df$Group[is.na(df$Group)] <- "not annotated"
        }
      }
      list(df = df, rank = cmp$rank, assay = cmp$assay, levels = ordering)
    })

    output$plt_stack <- plotly::renderPlotly({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_message_plot("No contigs pass the current filters."))
      }
      pd <- plot_data()
      if (is.null(pd)) {
        return(vpf_message_plot(paste(
          "Nothing to plot. Either no taxonomic rank carries data, or every",
          "contig was hidden by the unclassified switch."
        )))
      }
      ylab <- if (identical(input$scale, "relative")) {
        "Relative abundance"
      } else {
        paste0("Abundance (", pd$assay, ")")
      }
      p <- ggplot2::ggplot(pd$df, ggplot2::aes(
        x = .data$Sample, y = .data$Value, fill = .data$Taxon,
        text = paste0(.data$Taxon, "<br>", .data$Sample, "<br>",
                      signif(.data$Value, 4))
      )) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = stats::setNames(
          vpf_palette(length(pd$levels)), rev(levels(pd$df$Taxon))
        )) +
        ggplot2::labs(x = NULL, y = ylab, fill = pd$rank) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      if ("Group" %in% colnames(pd$df)) {
        p <- p + ggplot2::facet_wrap(~Group, scales = "free_x")
      }
      vpf_ggplotly(p, tooltip = "text")
    })

    output$stack_caption <- renderText({
      cmp <- composition()
      if (is.null(cmp)) {
        return("")
      }
      sprintf(paste("Rank %s, %s, %d contigs, of which %d carry no assignment at",
                    "this rank. This is a filtered subset of the assembly, so the",
                    "proportions describe the selected contigs, not the whole sample."),
              cmp$rank, cmp$assay, cmp$n_features, cmp$n_unclassified)
    })

    output$plt_ranks <- plotly::renderPlotly({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(vpf_message_plot("No contigs to summarize."))
      }
      rd <- SummarizedExperiment::rowData(tse)
      present <- intersect(VPF_TAXONOMY_RANKS, colnames(rd))
      if (length(present) == 0) {
        return(vpf_message_plot("This object carries no taxonomic ranks."))
      }
      counts <- vapply(present, function(r) sum(!is.na(vpf_blank_na(rd[[r]]))), numeric(1))
      df <- data.frame(
        Rank = factor(present, levels = present),
        Classified = as.numeric(counts),
        Percent = 100 * as.numeric(counts) / nrow(tse),
        stringsAsFactors = FALSE
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Rank, y = .data$Classified,
        text = paste0(.data$Rank, ": ", .data$Classified, " contigs (",
                      round(.data$Percent, 1), "%)")
      )) +
        ggplot2::geom_col(fill = "#59A14F") +
        ggplot2::labs(x = NULL, y = paste0("Contigs classified (of ", nrow(tse), ")")) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    taxon_table <- reactive({
      cmp <- composition()
      if (is.null(cmp) || nrow(cmp$matrix) == 0) {
        return(NULL)
      }
      mat <- cmp$matrix
      rel <- sweep(mat, 2, pmax(colSums(mat), .Machine$double.eps), "/")
      df <- data.frame(
        Taxon = rownames(mat),
        Total = rowSums(mat),
        Mean_relative_abundance = rowMeans(rel),
        Prevalence = rowSums(mat > 0) / ncol(mat),
        stringsAsFactors = FALSE
      )
      df <- df[order(-df$Total), , drop = FALSE]
      per_sample <- as.data.frame(mat[df$Taxon, , drop = FALSE])
      out <- cbind(df, per_sample)
      rownames(out) <- NULL
      out
    })

    mod_bigtable_server("tbl", r_data = taxon_table,
                        filename = reactive(paste0("composition_", active_rank() %||% "rank")))
  })
}

## To be copied in the UI
# mod_composition_ui("composition_1")

## To be copied in the server
# mod_composition_server("composition_1", r_filter)
