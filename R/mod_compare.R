#' compare UI Function
#'
#' @description Side-by-side comparison of the loaded dataset with a second
#'   one. Only quantities that survive being computed on two independently
#'   processed objects are shown: shared and unique taxa, per-sample diversity,
#'   and composition at a common rank. No between-dataset significance test is
#'   offered, because two runs differ in assembly, depth and filtering as much
#'   as in biology.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_compare_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4,
        vpf_card(
          "Second dataset",
          radioButtons(ns("source"), NULL,
                       choices = c("Bundled example" = "demo",
                                   "Upload an .rds file" = "upload",
                                   "Path on this server" = "path"),
                       selected = "demo"),
          conditionalPanel("input.source == 'demo'", ns = ns,
                           selectInput(ns("demo_choice"), "Example dataset",
                                       choices = character(0))),
          conditionalPanel("input.source == 'upload'", ns = ns,
                           fileInput(ns("upload"), "ViroProfiler output (.rds)",
                                     accept = c(".rds", ".RDS"))),
          conditionalPanel("input.source == 'path'", ns = ns,
                           uiOutput(ns("path_input"))),
          actionButton(ns("load"), "Load second dataset", class = "btn-primary"),
          uiOutput(ns("load_status"))
        ),
        vpf_card(
          "Comparison settings",
          uiOutput(ns("rank_ui")),
          selectInput(ns("alpha_index"), "Alpha-diversity index",
                      choices = c("Shannon" = "shannon",
                                  "Simpson (Gini-Simpson)" = "gini_simpson",
                                  "Observed richness" = "observed_richness"),
                      selected = "shannon")
        )
      ),
      column(
        8,
        uiOutput(ns("summary")),
        tabsetPanel(
          id = ns("tabs"),
          tabPanel(
            "Shared and unique taxa",
            plotly::plotlyOutput(ns("plt_overlap"), height = "360px"),
            reactable::reactableOutput(ns("tbl_overlap"))
          ),
          tabPanel(
            "Diversity",
            plotly::plotlyOutput(ns("plt_alpha"), height = "420px"),
            vpf_caption(
              "Per-sample diversity in both datasets. Differences reflect",
              "sequencing depth, assembly and filtering as well as biology, so",
              "these values are comparable only when both runs used the same",
              "pipeline settings."
            )
          ),
          tabPanel(
            "Composition",
            plotly::plotlyOutput(ns("plt_composition"), height = "460px")
          ),
          tabPanel(
            "Shared contigs",
            uiOutput(ns("contig_overlap"))
          )
        )
      )
    )
  )
}

#' compare Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_compare_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    demos <- vpf_demo_datasets()
    second <- reactiveVal(NULL)
    second_name <- reactiveVal(NULL)
    load_message <- reactiveVal(NULL)

    observe({
      if (length(demos) > 0) {
        updateSelectInput(session, "demo_choice", choices = names(demos),
                          selected = names(demos)[[length(demos)]])
      }
    })

    output$path_input <- renderUI({
      if (!vpf_server_path_allowed()) {
        return(vpf_notice(type = "info",
                          "Loading by server path is disabled in this deployment."))
      }
      textInput(ns("path"), "Absolute path to a .rds file")
    })

    observeEvent(input$load, {
      src <- input$source %||% "demo"
      spec <- switch(
        src,
        demo = {
          choice <- input$demo_choice
          if (is.null(choice) || !choice %in% names(demos)) {
            list(error = "No bundled example dataset is available.")
          } else {
            list(path = demos[[choice]], label = choice)
          }
        },
        upload = {
          up <- input$upload
          if (is.null(up)) list(error = "Choose a file to upload first.")
          else list(path = up$datapath, label = up$name)
        },
        path = {
          if (!vpf_server_path_allowed()) {
            list(error = "Loading by server path is disabled in this deployment.")
          } else if (is.null(input$path) || !nzchar(trimws(input$path))) {
            list(error = "Enter a file path first.")
          } else {
            list(path = trimws(input$path), label = basename(trimws(input$path)))
          }
        },
        list(error = "Unknown data source.")
      )
      if (!is.null(spec$error)) {
        load_message(list(type = "danger", text = spec$error))
        return()
      }
      res <- vpf_read_tse(spec$path)
      if (!is.null(res$error)) {
        load_message(list(type = "danger", text = res$error))
        return()
      }
      second(res$tse)
      second_name(spec$label)
      load_message(list(
        type = "success",
        text = sprintf("Loaded %d contigs x %d samples from %s.",
                       nrow(res$tse), ncol(res$tse), spec$label)
      ))
    })

    output$load_status <- renderUI({
      msg <- load_message()
      if (is.null(msg)) {
        return(vpf_caption("Nothing loaded yet."))
      }
      vpf_notice(type = msg$type, msg$text)
    })

    labels <- reactive({
      list(a = r_filter$name() %||% "Dataset 1",
           b = second_name() %||% "Dataset 2")
    })

    common_ranks <- reactive({
      a <- r_filter$tse()
      b <- second()
      if (is.null(a) || is.null(b)) {
        return(character(0))
      }
      intersect(vpf_available_ranks(a), vpf_available_ranks(b))
    })

    output$rank_ui <- renderUI({
      rk <- common_ranks()
      if (length(rk) == 0) {
        return(vpf_caption("Load a second dataset to choose a common rank."))
      }
      default <- if ("Family" %in% rk) "Family" else rk[[length(rk)]]
      selectInput(ns("rank"), "Rank compared", choices = rk, selected = default)
    })

    active_rank <- reactive({
      rk <- common_ranks()
      if (length(rk) == 0) {
        return(NULL)
      }
      sel <- input$rank
      if (is.null(sel) || !sel %in% rk) rk[[length(rk)]] else sel
    })

    output$summary <- renderUI({
      a <- r_filter$tse()
      b <- second()
      if (is.null(a)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      if (is.null(b)) {
        return(vpf_notice(
          type = "info",
          "Load a second dataset on the left to compare. The bundled examples",
          "need no files."
        ))
      }
      lb <- labels()
      shared_assays <- intersect(SummarizedExperiment::assayNames(a),
                                 SummarizedExperiment::assayNames(b))
      tagList(
        fluidRow(
          column(6, vpf_card(
            lb$a,
            div(vpf_stat("Contigs", nrow(a)), vpf_stat("Samples", ncol(a)),
                vpf_stat("Ranks", length(vpf_available_ranks(a))))
          )),
          column(6, vpf_card(
            lb$b,
            div(vpf_stat("Contigs", nrow(b)), vpf_stat("Samples", ncol(b)),
                vpf_stat("Ranks", length(vpf_available_ranks(b))))
          ))
        ),
        if (length(shared_assays) == 0) {
          vpf_notice(type = "danger",
                     "The two objects share no assay, so quantitative comparison is impossible.")
        } else {
          vpf_caption(paste0("Shared assays: ", paste(shared_assays, collapse = ", "), "."))
        },
        if (length(common_ranks()) == 0) {
          vpf_notice(type = "warning",
                     "The two objects share no taxonomic rank that carries data.")
        }
      )
    })

    overlap <- reactive({
      a <- r_filter$tse()
      b <- second()
      rank <- active_rank()
      if (is.null(a) || is.null(b) || is.null(rank)) {
        return(NULL)
      }
      ta <- unique(stats::na.omit(vpf_blank_na(SummarizedExperiment::rowData(a)[[rank]])))
      tb <- unique(stats::na.omit(vpf_blank_na(SummarizedExperiment::rowData(b)[[rank]])))
      list(shared = intersect(ta, tb), only_a = setdiff(ta, tb),
           only_b = setdiff(tb, ta), rank = rank)
    })

    output$plt_overlap <- plotly::renderPlotly({
      ov <- overlap()
      if (is.null(ov)) {
        return(vpf_message_plot(paste(
          "Load a second dataset and choose a rank that both objects carry."
        )))
      }
      lb <- labels()
      df <- data.frame(
        Category = factor(
          c("Shared", paste0("Only in ", lb$a), paste0("Only in ", lb$b)),
          levels = c("Shared", paste0("Only in ", lb$a), paste0("Only in ", lb$b))
        ),
        Count = c(length(ov$shared), length(ov$only_a), length(ov$only_b)),
        stringsAsFactors = FALSE
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Category, y = .data$Count, fill = .data$Category,
        text = paste0(.data$Category, ": ", .data$Count, " taxa")
      )) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::scale_fill_manual(values = c("#59A14F", "#4E79A7", "#F28E2B")) +
        ggplot2::labs(x = NULL, y = paste0("Distinct ", ov$rank, " taxa")) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$tbl_overlap <- reactable::renderReactable({
      ov <- overlap()
      if (is.null(ov)) {
        return(reactable::reactable(data.frame(Message = "Load a second dataset first.")))
      }
      lb <- labels()
      n <- max(length(ov$shared), length(ov$only_a), length(ov$only_b), 1)
      pad <- function(x) c(sort(x), rep(NA_character_, n - length(x)))
      df <- data.frame(
        Shared = pad(ov$shared),
        A = pad(ov$only_a),
        B = pad(ov$only_b),
        stringsAsFactors = FALSE
      )
      colnames(df) <- c("Shared", paste0("Only in ", lb$a), paste0("Only in ", lb$b))
      reactable::reactable(df, striped = TRUE, bordered = TRUE, highlight = TRUE,
                           searchable = TRUE, defaultPageSize = 12)
    })

    alpha_frame <- function(tse, label, index) {
      if (is.null(tse) || nrow(tse) == 0 || ncol(tse) == 0) {
        return(NULL)
      }
      assay_name <- vpf_default_assay(tse)
      if (is.null(assay_name)) {
        return(NULL)
      }
      res <- tryCatch(
        mia::addAlpha(tse, assay.type = assay_name, index = index, name = "alpha_value"),
        error = function(e) NULL
      )
      if (is.null(res)) {
        return(NULL)
      }
      data.frame(
        Sample = vpf_sample_names(res),
        Value = as.numeric(SummarizedExperiment::colData(res)[["alpha_value"]]),
        Dataset = label,
        stringsAsFactors = FALSE
      )
    }

    output$plt_alpha <- plotly::renderPlotly({
      a <- r_filter$tse()
      b <- second()
      if (is.null(a) || is.null(b)) {
        return(vpf_message_plot("Load a second dataset to compare diversity."))
      }
      lb <- labels()
      index <- input$alpha_index %||% "shannon"
      df <- rbind(alpha_frame(a, lb$a, index), alpha_frame(b, lb$b, index))
      if (is.null(df) || nrow(df) == 0) {
        return(vpf_message_plot("Alpha diversity could not be computed for either dataset."))
      }
      df$Sample <- factor(df$Sample, levels = unique(df$Sample))
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Sample, y = .data$Value, fill = .data$Dataset,
        text = paste0(.data$Sample, " (", .data$Dataset, ")<br>",
                      signif(.data$Value, 4))
      )) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
        ggplot2::facet_wrap(~Dataset, scales = "free_x") +
        ggplot2::labs(x = NULL, y = index) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      vpf_ggplotly(p, tooltip = "text")
    })

    output$plt_composition <- plotly::renderPlotly({
      a <- r_filter$tse()
      b <- second()
      rank <- active_rank()
      if (is.null(a) || is.null(b) || is.null(rank)) {
        return(vpf_message_plot("Load a second dataset and choose a common rank."))
      }
      lb <- labels()
      share <- function(tse, label) {
        assay_name <- vpf_default_assay(tse)
        if (is.null(assay_name) || nrow(tse) == 0) {
          return(NULL)
        }
        lab <- vpf_blank_na(SummarizedExperiment::rowData(tse)[[rank]])
        lab[is.na(lab)] <- "Unclassified"
        mat <- SummarizedExperiment::assay(tse, assay_name)
        mat[is.na(mat)] <- 0
        agg <- rowsum(as.matrix(mat), group = lab, reorder = TRUE)
        totals <- rowSums(agg)
        data.frame(Taxon = names(totals),
                   Share = as.numeric(totals) / max(sum(totals), .Machine$double.eps),
                   Dataset = label, stringsAsFactors = FALSE)
      }
      df <- rbind(share(a, lb$a), share(b, lb$b))
      if (is.null(df) || nrow(df) == 0) {
        return(vpf_message_plot("No composition could be computed."))
      }
      top <- utils::head(names(sort(tapply(df$Share, df$Taxon, sum), decreasing = TRUE)), 15)
      df$Taxon[!df$Taxon %in% top] <- "Other"
      df <- stats::aggregate(df$Share, by = list(Taxon = df$Taxon, Dataset = df$Dataset),
                             FUN = sum)
      colnames(df)[3] <- "Share"
      ordering <- c(setdiff(top, c("Other", "Unclassified")),
                    intersect(c("Unclassified", "Other"), unique(df$Taxon)))
      df$Taxon <- factor(df$Taxon, levels = rev(unique(c(ordering, unique(df$Taxon)))))
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Dataset, y = .data$Share, fill = .data$Taxon,
        text = paste0(.data$Taxon, "<br>", round(100 * .data$Share, 2), "%")
      )) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = stats::setNames(
          vpf_palette(nlevels(df$Taxon)), levels(df$Taxon)
        )) +
        ggplot2::labs(x = NULL, y = "Share of total abundance", fill = rank) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$contig_overlap <- renderUI({
      a <- r_filter$tse()
      b <- second()
      if (is.null(a) || is.null(b)) {
        return(vpf_notice(type = "info", "Load a second dataset first."))
      }
      ia <- rownames(a)
      ib <- rownames(b)
      shared <- intersect(ia, ib)
      if (length(shared) == 0) {
        return(vpf_notice(
          type = "info",
          title = "No contig identifier occurs in both objects",
          paste("Contig names are assigned per assembly, so two independent runs",
                "share none unless the same vOTU catalogue was used. Compare at a",
                "taxonomic rank instead, or dereplicate both runs against a common",
                "vOTU set before comparing contigs.")
        ))
      }
      tagList(
        vpf_notice(type = "success",
                   sprintf("%d contig identifiers occur in both datasets.", length(shared))),
        reactable::reactable(data.frame(Contig = shared, stringsAsFactors = FALSE),
                             searchable = TRUE, striped = TRUE, bordered = TRUE,
                             defaultPageSize = 15)
      )
    })
  })
}

## To be copied in the UI
# mod_compare_ui("compare_1")

## To be copied in the server
# mod_compare_server("compare_1", r_filter)
