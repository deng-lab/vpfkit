#' host UI Function
#'
#' @description Predicted bacterial hosts and predicted replication cycle. Both
#'   are model outputs rather than observations, which the panel states rather
#'   than implies.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_host_ui <- function(id) {
  ns <- NS(id)
  tagList(
    tabsetPanel(
      id = ns("tabs"),
      tabPanel(
        "Host prediction",
        uiOutput(ns("host_status")),
        fluidRow(
          column(7, plotly::plotlyOutput(ns("plt_host"), height = "440px")),
          column(5, plotly::plotlyOutput(ns("plt_host_score"), height = "440px"))
        ),
        uiOutput(ns("host_controls")),
        reactable::reactableOutput(ns("tbl_host"))
      ),
      tabPanel(
        "Replication cycle",
        uiOutput(ns("cycle_status")),
        fluidRow(
          column(6, plotly::plotlyOutput(ns("plt_cycle"), height = "400px")),
          column(6, plotly::plotlyOutput(ns("plt_cycle_abundance"), height = "400px"))
        ),
        vpf_caption(
          "The right panel weights each lifestyle call by abundance, which can",
          "differ sharply from the contig counts on the left when a few abundant",
          "contigs dominate."
        ),
        reactable::reactableOutput(ns("tbl_cycle"))
      ),
      tabPanel(
        "Provirus and integration",
        uiOutput(ns("provirus_status")),
        plotly::plotlyOutput(ns("plt_provirus"), height = "380px")
      )
    )
  )
}

#' host Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_host_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    rows <- reactive({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(NULL)
      }
      vpf_row_table(tse)
    })

    has_col <- function(df, nm) {
      !is.null(df) && nm %in% colnames(df) && any(!is.na(vpf_blank_na(df[[nm]])))
    }

    host_column <- reactive({
      df <- rows()
      for (nm in c("iphop_genus", "phist_host", "iphop_host", "host_genus")) {
        if (has_col(df, nm)) {
          return(nm)
        }
      }
      NULL
    })

    cycle_column <- reactive({
      df <- rows()
      for (nm in c("bacphlip_replicyc", "replidec_replicyc", "vibrant_replicyc",
                   "replication_cycle")) {
        if (has_col(df, nm)) {
          return(nm)
        }
      }
      NULL
    })

    output$host_status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_notice(type = "warning", "No contigs pass the current filters."))
      }
      if (is.null(host_column())) {
        return(vpf_missing_annotation(
          "Host prediction",
          "--use_iphop (default true)",
          paste("iPHoP predicts a host genus for each viral contig, together with",
                "a confidence score and the methods that contributed to the call.",
                "This dataset carries no iphop_* or phist_* column.")
        ))
      }
      vpf_notice(
        type = "info",
        title = "Predicted hosts, not observed infections",
        paste0("Host calls come from ", host_column(),
               ". A prediction links a viral sequence to a candidate host taxon by",
               " sequence evidence; it is not a demonstration that the phage",
               " infects that host.")
      )
    })

    output$host_controls <- renderUI({
      if (is.null(host_column())) {
        return(NULL)
      }
      df <- rows()
      if (has_col(df, "iphop_score")) {
        sliderInput(ns("min_host_score"), "Minimum host-prediction score",
                    min = 0, max = 100, value = 0, step = 1, width = "50%")
      } else {
        NULL
      }
    })

    host_data <- reactive({
      col <- host_column()
      df <- rows()
      if (is.null(col) || is.null(df)) {
        return(NULL)
      }
      out <- data.frame(
        Contig = df$Contig,
        Host = vpf_blank_na(df[[col]]),
        stringsAsFactors = FALSE
      )
      if ("iphop_score" %in% colnames(df)) {
        out$Score <- vpf_as_numeric_column(df$iphop_score)
      } else if ("phist_score" %in% colnames(df)) {
        out$Score <- vpf_as_numeric_column(df$phist_score)
      }
      if ("iphop_methods" %in% colnames(df)) {
        out$Methods <- vpf_blank_na(df$iphop_methods)
      }
      thr <- input$min_host_score
      if (!is.null(thr) && "Score" %in% colnames(out)) {
        out <- out[is.na(out$Score) | out$Score >= thr, , drop = FALSE]
      }
      out[!is.na(out$Host), , drop = FALSE]
    })

    output$plt_host <- plotly::renderPlotly({
      if (is.null(host_column())) {
        return(vpf_message_plot(paste(
          "No host prediction in this dataset. It is produced by the",
          "ViroProfiler option --use_iphop."
        )))
      }
      hd <- host_data()
      if (is.null(hd) || nrow(hd) == 0) {
        return(vpf_message_plot("No contig has a host prediction above the current score threshold."))
      }
      counts <- as.data.frame(table(Host = hd$Host), stringsAsFactors = FALSE)
      counts <- counts[order(-counts$Freq), , drop = FALSE]
      counts <- utils::head(counts, 25)
      counts$Host <- factor(counts$Host, levels = rev(counts$Host))
      p <- ggplot2::ggplot(counts, ggplot2::aes(
        x = .data$Host, y = .data$Freq,
        text = paste0(.data$Host, ": ", .data$Freq, " contigs")
      )) +
        ggplot2::geom_col(fill = "#76B7B2") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Contigs", title = "Predicted host taxa") +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$plt_host_score <- plotly::renderPlotly({
      hd <- host_data()
      if (is.null(hd) || !"Score" %in% colnames(hd)) {
        return(vpf_message_plot(
          "This dataset stores no host-prediction confidence score."
        ))
      }
      sc <- hd$Score[is.finite(hd$Score)]
      if (length(sc) == 0) {
        return(vpf_message_plot("No usable host-prediction scores."))
      }
      df <- data.frame(Score = sc)
      p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Score)) +
        ggplot2::geom_histogram(bins = min(30, max(5, length(unique(sc)))),
                                fill = "#4E79A7", alpha = 0.85) +
        ggplot2::labs(x = "Host-prediction score", y = "Contigs",
                      title = "Confidence distribution") +
        ggplot2::theme_bw()
      vpf_ggplotly(p)
    })

    output$tbl_host <- reactable::renderReactable({
      hd <- host_data()
      if (is.null(hd) || nrow(hd) == 0) {
        return(reactable::reactable(data.frame(
          Message = "No host predictions for the current selection."
        )))
      }
      reactable::reactable(vpf_round_df(hd), searchable = TRUE, striped = TRUE,
                           bordered = TRUE, highlight = TRUE, defaultPageSize = 12,
                           resizable = TRUE, wrap = FALSE)
    })

    # -- Replication cycle ---------------------------------------------------
    output$cycle_status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (nrow(tse) == 0) {
        return(vpf_notice(type = "warning", "No contigs pass the current filters."))
      }
      if (is.null(cycle_column())) {
        return(vpf_missing_annotation(
          "Replication cycle",
          "--replicyc bacphlip | replidec",
          paste("BACPHLIP or Replidec classify each viral contig as temperate or",
                "virulent. This dataset carries no bacphlip_replicyc,",
                "replidec_replicyc or vibrant_replicyc column.")
        ))
      }
      vpf_notice(
        type = "info",
        paste0("Lifestyle calls come from ", cycle_column(),
               ". These classifiers are trained on complete phage genomes, so a",
               " call on a short or incomplete contig is far less reliable than",
               " one on a complete genome.")
      )
    })

    cycle_data <- reactive({
      col <- cycle_column()
      df <- rows()
      tse <- r_filter$tse()
      assay_name <- r_filter$assay()
      if (is.null(col) || is.null(df) || is.null(tse)) {
        return(NULL)
      }
      out <- data.frame(
        Contig = df$Contig,
        Cycle = vpf_blank_na(df[[col]]),
        stringsAsFactors = FALSE
      )
      out$Cycle[is.na(out$Cycle)] <- "not predicted"
      if (!is.null(assay_name) && assay_name %in% SummarizedExperiment::assayNames(tse)) {
        out$Abundance <- rowSums(SummarizedExperiment::assay(tse, assay_name), na.rm = TRUE)
      }
      if ("checkv_quality" %in% colnames(df)) {
        out$Quality <- vpf_blank_na(df$checkv_quality)
      }
      out
    })

    output$plt_cycle <- plotly::renderPlotly({
      if (is.null(cycle_column())) {
        return(vpf_message_plot(paste(
          "No replication-cycle prediction in this dataset. It is produced by the",
          "ViroProfiler option --replicyc bacphlip or --replicyc replidec."
        )))
      }
      cd <- cycle_data()
      if (is.null(cd) || nrow(cd) == 0) {
        return(vpf_message_plot("No contigs in the current selection."))
      }
      if ("Quality" %in% colnames(cd)) {
        cd$Quality[is.na(cd$Quality)] <- "not annotated"
        counts <- as.data.frame(table(Cycle = cd$Cycle, Quality = cd$Quality),
                                stringsAsFactors = FALSE)
        p <- ggplot2::ggplot(counts, ggplot2::aes(
          x = .data$Cycle, y = .data$Freq, fill = .data$Quality,
          text = paste0(.data$Cycle, " / ", .data$Quality, ": ", .data$Freq)
        )) +
          ggplot2::geom_col() +
          ggplot2::scale_fill_manual(values = stats::setNames(
            vpf_palette(length(unique(counts$Quality))), unique(counts$Quality)
          ))
      } else {
        counts <- as.data.frame(table(Cycle = cd$Cycle), stringsAsFactors = FALSE)
        p <- ggplot2::ggplot(counts, ggplot2::aes(
          x = .data$Cycle, y = .data$Freq,
          text = paste0(.data$Cycle, ": ", .data$Freq)
        )) +
          ggplot2::geom_col(fill = "#B07AA1")
      }
      p <- p +
        ggplot2::labs(x = NULL, y = "Contigs", title = "Predicted lifestyle") +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$plt_cycle_abundance <- plotly::renderPlotly({
      cd <- cycle_data()
      if (is.null(cd) || !"Abundance" %in% colnames(cd)) {
        return(vpf_message_plot("No abundance assay available for weighting."))
      }
      agg <- stats::aggregate(cd$Abundance, by = list(Cycle = cd$Cycle), FUN = sum)
      colnames(agg)[2] <- "Abundance"
      total <- sum(agg$Abundance)
      agg$Share <- if (total > 0) agg$Abundance / total else 0
      p <- ggplot2::ggplot(agg, ggplot2::aes(
        x = .data$Cycle, y = .data$Share,
        text = paste0(.data$Cycle, ": ", round(100 * .data$Share, 1), "% of total abundance")
      )) +
        ggplot2::geom_col(fill = "#F28E2B") +
        ggplot2::labs(x = NULL, y = "Share of total abundance",
                      title = "Lifestyle weighted by abundance") +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$tbl_cycle <- reactable::renderReactable({
      cd <- cycle_data()
      if (is.null(cd) || nrow(cd) == 0) {
        return(reactable::reactable(data.frame(
          Message = "No replication-cycle predictions for the current selection."
        )))
      }
      reactable::reactable(vpf_round_df(cd), searchable = TRUE, striped = TRUE,
                           bordered = TRUE, highlight = TRUE, defaultPageSize = 12,
                           resizable = TRUE, wrap = FALSE)
    })

    # -- Provirus ------------------------------------------------------------
    provirus_column <- reactive({
      df <- rows()
      for (nm in c("checkv_provirus", "genomad_topology")) {
        if (has_col(df, nm)) {
          return(nm)
        }
      }
      NULL
    })

    output$provirus_status <- renderUI({
      tse <- r_filter$tse()
      if (is.null(tse) || nrow(tse) == 0) {
        return(vpf_notice(type = "info", "No contigs in the current selection."))
      }
      if (is.null(provirus_column())) {
        return(vpf_missing_annotation(
          "Provirus status",
          "always on (CHECKV and GENOMAD processes)",
          "CheckV reports whether a contig is an integrated provirus, and geNomad reports its topology."
        ))
      }
      vpf_caption(paste0("Based on ", provirus_column(), "."))
    })

    output$plt_provirus <- plotly::renderPlotly({
      col <- provirus_column()
      df <- rows()
      if (is.null(col) || is.null(df)) {
        return(vpf_message_plot(paste(
          "No provirus or topology annotation in this dataset. CheckV and geNomad",
          "produce it as part of the always-on detection stage."
        )))
      }
      v <- vpf_blank_na(df[[col]])
      v[is.na(v)] <- "not annotated"
      counts <- as.data.frame(table(Status = v), stringsAsFactors = FALSE)
      counts <- counts[order(-counts$Freq), , drop = FALSE]
      counts$Status <- factor(counts$Status, levels = rev(counts$Status))
      p <- ggplot2::ggplot(counts, ggplot2::aes(
        x = .data$Status, y = .data$Freq,
        text = paste0(.data$Status, ": ", .data$Freq, " contigs")
      )) +
        ggplot2::geom_col(fill = "#59A14F") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = "Contigs", title = col) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })
  })
}

## To be copied in the UI
# mod_host_ui("host_1")

## To be copied in the server
# mod_host_server("host_1", r_filter)
