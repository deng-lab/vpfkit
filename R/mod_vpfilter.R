#' vpfilter UI Function
#'
#' @description Contig selection. Every control is built from the columns that
#'   are actually present, so a filter can never be applied against a missing
#'   annotation, and the module reports how many contigs each step removed.
#'   The filtered object it returns is the single source of truth for every
#'   other tab.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_vpfilter_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4,
        vpf_card(
          "Abundance metric",
          uiOutput(ns("assay_ui")),
          uiOutput(ns("assay_description")),
          uiOutput(ns("covfrac_ui"))
        ),
        vpf_card(
          "Contig filters",
          uiOutput(ns("filter_controls")),
          div(
            style = "margin-top: 10px;",
            actionButton(ns("reset"), "Reset all filters", class = "btn-sm")
          )
        )
      ),
      column(
        8,
        uiOutput(ns("summary_box")),
        tabsetPanel(
          id = ns("tabs"),
          tabPanel(
            "Filter audit",
            vpf_caption(
              "Each step is applied to the contigs that survived the previous one.",
              "Filters whose annotation is missing from this dataset are listed as",
              "not applied rather than silently dropping every contig."
            ),
            reactable::reactableOutput(ns("tbl_audit"))
          ),
          tabPanel(
            "Quality distributions",
            fluidRow(
              column(6, plotly::plotlyOutput(ns("plt_length"), height = "320px")),
              column(6, plotly::plotlyOutput(ns("plt_quality"), height = "320px"))
            ),
            fluidRow(
              column(6, plotly::plotlyOutput(ns("plt_completeness"), height = "320px")),
              column(6, plotly::plotlyOutput(ns("plt_score"), height = "320px"))
            ),
            vpf_caption(
              "Kept and removed contigs are shown together so the effect of the",
              "current thresholds stays visible."
            )
          ),
          tabPanel(
            "Viral evidence",
            uiOutput(ns("evidence_status")),
            fluidRow(
              column(6, plotly::plotlyOutput(ns("plt_votes"), height = "340px")),
              column(6, plotly::plotlyOutput(ns("plt_vote_sources"), height = "340px"))
            ),
            vpf_caption(
              "A contig enters the viral set if any single detector calls it,",
              "so the vote count is a sensitivity ranking rather than a consensus.",
              "A contig supported by one low-quality call contributes to diversity",
              "exactly as much as one supported by all of them, which is why the",
              "minimum-votes filter on the left is worth using."
            )
          ),
          tabPanel(
            "Prevalence and detection",
            plotly::plotlyOutput(ns("plt_prevalence"), height = "380px"),
            vpf_caption(
              "Prevalence is the proportion of samples in which a contig is",
              "detected after coverage masking. Reading a non-zero abundance as",
              "confident detection is the most common error in virome tables;",
              "the covered-fraction threshold above is what guards against it."
            )
          )
        )
      )
    )
  )
}

#' vpfilter Server Functions
#'
#' @param id Module id.
#' @param r_input A list of reactives from [mod_data_input_server()].
#'
#' @return A list of reactives: `tse` (filtered object), `raw` (unfiltered),
#'   `assay` (active abundance assay name), `audit` (per-step audit table) and
#'   `n_removed`.
#' @noRd
mod_vpfilter_server <- function(id, r_input) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    raw_tse <- reactive({
      r_input$tse()
    })

    # -- Abundance assay ----------------------------------------------------
    output$assay_ui <- renderUI({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(vpf_caption("Load a dataset first."))
      }
      choices <- vpf_abundance_assay_choices(tse)
      if (length(choices) == 0) {
        return(vpf_notice(type = "danger",
                          "This object has no assay that can be used as an abundance."))
      }
      selectInput(ns("assay"), "Assay used for every downstream analysis",
                  choices = choices, selected = vpf_default_assay(tse))
    })

    active_assay <- reactive({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(NULL)
      }
      choices <- unname(vpf_abundance_assay_choices(tse))
      sel <- input$assay
      if (is.null(sel) || !sel %in% choices) {
        return(vpf_default_assay(tse))
      }
      sel
    })

    output$assay_description <- renderUI({
      tse <- raw_tse()
      assay_name <- active_assay()
      if (is.null(tse) || is.null(assay_name)) {
        return(NULL)
      }
      info <- vpf_assay_info(tse)
      row <- info[info$assay == assay_name, , drop = FALSE]
      if (nrow(row) == 0) {
        return(NULL)
      }
      type <- if (identical(assay_name, "tmm")) "warning" else "info"
      vpf_notice(
        type = type,
        title = paste0(row$label[[1]], "  [unit: ", row$unit[[1]], "]"),
        row$description[[1]]
      )
    })

    output$covfrac_ui <- renderUI({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(NULL)
      }
      cf <- vpf_covfrac_assay(tse)
      if (is.null(cf)) {
        return(vpf_notice(
          type = "info",
          "This object has no covered-fraction assay, so detection masking is",
          "unavailable. Abundances are used exactly as stored."
        ))
      }
      tagList(
        sliderInput(ns("min_covfrac"), "Minimum covered fraction to call a contig present",
                    min = 0, max = 1, value = 0.75, step = 0.05),
        vpf_caption(
          "Abundance is set to zero in every sample where the contig's covered",
          "fraction falls below this value. Reads from related genomes often map",
          "over a short stretch only, which would otherwise look like presence.",
          "Missing coverage counts as absent."
        )
      )
    })

    # -- Filter controls, built from the columns that exist ------------------
    filterable <- reactive({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(list())
      }
      rd <- SummarizedExperiment::rowData(tse)
      out <- list()

      len <- vpf_as_numeric_column(rd[["checkv_contig_length"]])
      if (!is.null(len) && any(is.finite(len))) {
        out$length <- list(min = min(len, na.rm = TRUE), max = max(len, na.rm = TRUE))
      }
      if ("checkv_quality" %in% colnames(rd)) {
        lv <- unique(stats::na.omit(vpf_blank_na(rd[["checkv_quality"]])))
        known <- c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined")
        lv <- c(intersect(known, lv), setdiff(lv, known))
        if (length(lv) > 0) out$quality <- list(levels = lv)
      }
      comp <- vpf_as_numeric_column(rd[["checkv_completeness"]])
      if (!is.null(comp) && any(is.finite(comp))) {
        out$completeness <- list(min = 0, max = 100)
      }
      gs <- vpf_as_numeric_column(rd[["genomad_score"]])
      if (!is.null(gs) && any(is.finite(gs))) {
        out$genomad <- list(min = 0, max = 1)
      }
      vs <- vpf_as_numeric_column(rd[["virsorter2_max_score"]])
      if (!is.null(vs) && any(is.finite(vs))) {
        out$virsorter2 <- list(min = 0, max = 1)
      }
      votes <- vpf_viral_votes(tse)
      if (!is.null(votes) && any(is.finite(votes))) {
        out$votes <- list(max = max(votes, na.rm = TRUE))
      }
      out
    })

    output$filter_controls <- renderUI({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(vpf_caption("Load a dataset first."))
      }
      f <- filterable()
      controls <- list()

      if (!is.null(f$length)) {
        controls <- c(controls, list(
          numericInput(ns("min_length"), "Minimum contig length (bp, inclusive)",
                       value = 0, min = 0, step = 1000),
          vpf_caption(sprintf("Observed range: %s to %s bp.",
                              format(f$length$min, big.mark = ","),
                              format(f$length$max, big.mark = ",")))
        ))
      }
      if (!is.null(f$quality)) {
        controls <- c(controls, list(
          selectInput(ns("quality"), "CheckV quality tiers to keep",
                      choices = f$quality$levels, selected = f$quality$levels,
                      multiple = TRUE),
          vpf_caption("Clearing the selection keeps every tier.")
        ))
      }
      if (!is.null(f$completeness)) {
        controls <- c(controls, list(
          sliderInput(ns("min_completeness"), "Minimum CheckV completeness (%)",
                      min = 0, max = 100, value = 0, step = 1),
          checkboxInput(ns("completeness_keep_na"),
                        "Keep contigs whose completeness is missing", value = TRUE)
        ))
      }
      if (!is.null(f$genomad)) {
        controls <- c(controls, list(
          sliderInput(ns("min_genomad"), "Minimum geNomad virus score",
                      min = 0, max = 1, value = 0, step = 0.01),
          checkboxInput(ns("genomad_keep_na"),
                        "Keep contigs with no geNomad score", value = TRUE)
        ))
      }
      if (!is.null(f$virsorter2)) {
        controls <- c(controls, list(
          sliderInput(ns("min_virsorter2"), "Minimum VirSorter2 max score",
                      min = 0, max = 1, value = 0, step = 0.01),
          checkboxInput(ns("virsorter2_keep_na"),
                        "Keep contigs with no VirSorter2 score", value = TRUE)
        ))
      }

      if (!is.null(f$votes)) {
        controls <- c(controls, list(
          sliderInput(ns("min_votes"), "Minimum number of viral-evidence votes",
                      min = 1, max = max(2, f$votes$max), value = 1, step = 1),
          vpf_caption(
            "How many independent detectors had to call the contig viral.",
            "A contig admitted on a single low-quality call is a common source of",
            "noise in diversity estimates."
          )
        ))
      }

      controls <- c(controls, list(
        tags$hr(),
        sliderInput(ns("min_prevalence"), "Minimum prevalence (fraction of samples detected)",
                    min = 0, max = 1, value = 0, step = 0.05),
        numericInput(ns("min_total_abundance"), "Minimum total abundance across samples",
                     value = 0, min = 0)
      ))

      missing <- setdiff(c("length", "quality", "completeness", "genomad", "virsorter2"),
                         names(f))
      if (length(missing) > 0) {
        labels <- c(length = "contig length", quality = "CheckV quality",
                    completeness = "CheckV completeness",
                    genomad = "geNomad score", virsorter2 = "VirSorter2 score")
        controls <- c(controls, list(vpf_notice(
          type = "info",
          paste0("Not available in this dataset, so not offered as a filter: ",
                 paste(labels[missing], collapse = ", "), ".")
        )))
      }
      do.call(tagList, controls)
    })

    observeEvent(input$reset, {
      f <- filterable()
      if (!is.null(f$length)) updateNumericInput(session, "min_length", value = 0)
      if (!is.null(f$quality)) {
        updateSelectInput(session, "quality", selected = f$quality$levels)
      }
      if (!is.null(f$completeness)) updateSliderInput(session, "min_completeness", value = 0)
      if (!is.null(f$genomad)) updateSliderInput(session, "min_genomad", value = 0)
      if (!is.null(f$virsorter2)) updateSliderInput(session, "min_virsorter2", value = 0)
      if (!is.null(f$votes)) updateSliderInput(session, "min_votes", value = 1)
      updateSliderInput(session, "min_prevalence", value = 0)
      updateNumericInput(session, "min_total_abundance", value = 0)
      updateSliderInput(session, "min_covfrac", value = 0.75)
    }, ignoreInit = TRUE)

    # -- The filtering itself ------------------------------------------------
    num_input <- function(x, default) {
      if (is.null(x) || length(x) != 1 || !is.finite(suppressWarnings(as.numeric(x)))) {
        return(default)
      }
      as.numeric(x)
    }

    filtered <- reactive({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(list(tse = NULL, audit = NULL))
      }
      assay_name <- active_assay()
      if (is.null(assay_name)) {
        return(list(tse = NULL, audit = NULL))
      }

      audit <- list()
      note <- function(step, applied, before, after, detail = "") {
        audit[[length(audit) + 1]] <<- data.frame(
          Step = step, Applied = applied,
          Removed = before - after, Remaining = after,
          Detail = detail, stringsAsFactors = FALSE
        )
      }
      n0 <- nrow(tse)
      note("Loaded", "-", n0, n0, "")

      rd <- SummarizedExperiment::rowData(tse)
      f <- filterable()

      # Contig length
      if (!is.null(f$length)) {
        thr <- num_input(input$min_length, 0)
        len <- vpf_as_numeric_column(rd[["checkv_contig_length"]])
        before <- nrow(tse)
        tse <- vpf_subset_rows(tse, len >= thr, na_keeps = FALSE)
        rd <- SummarizedExperiment::rowData(tse)
        note("Contig length", "yes", before, nrow(tse), paste0(">= ", format(thr, big.mark = ","), " bp"))
      } else {
        note("Contig length", "no", nrow(tse), nrow(tse), "checkv_contig_length absent")
      }

      # CheckV quality
      if (!is.null(f$quality)) {
        sel <- input$quality
        if (is.null(sel) || length(sel) == 0) {
          note("CheckV quality", "no", nrow(tse), nrow(tse), "no tier selected: all kept")
        } else {
          before <- nrow(tse)
          q <- vpf_blank_na(rd[["checkv_quality"]])
          tse <- vpf_subset_rows(tse, q %in% sel, na_keeps = FALSE)
          rd <- SummarizedExperiment::rowData(tse)
          note("CheckV quality", "yes", before, nrow(tse), paste(sel, collapse = ", "))
        }
      } else {
        note("CheckV quality", "no", nrow(tse), nrow(tse), "checkv_quality absent")
      }

      numeric_filter <- function(column, threshold, keep_na, label, detail) {
        before <- nrow(tse)
        v <- vpf_as_numeric_column(rd[[column]])
        if (is.null(v)) {
          note(label, "no", before, before, paste0(column, " is not numeric"))
          return(invisible())
        }
        keep <- v >= threshold
        keep[is.na(v)] <- isTRUE(keep_na)
        tse <<- vpf_subset_rows(tse, keep, na_keeps = isTRUE(keep_na))
        rd <<- SummarizedExperiment::rowData(tse)
        note(label, "yes", before, nrow(tse), detail)
      }

      if (!is.null(f$completeness)) {
        thr <- num_input(input$min_completeness, 0)
        numeric_filter("checkv_completeness", thr,
                       input$completeness_keep_na %||% TRUE,
                       "CheckV completeness", paste0(">= ", thr, "%"))
      } else {
        note("CheckV completeness", "no", nrow(tse), nrow(tse), "checkv_completeness absent")
      }
      if (!is.null(f$genomad)) {
        thr <- num_input(input$min_genomad, 0)
        numeric_filter("genomad_score", thr, input$genomad_keep_na %||% TRUE,
                       "geNomad score", paste0(">= ", thr))
      } else {
        note("geNomad score", "no", nrow(tse), nrow(tse), "genomad_score absent")
      }
      if (!is.null(f$virsorter2)) {
        thr <- num_input(input$min_virsorter2, 0)
        numeric_filter("virsorter2_max_score", thr, input$virsorter2_keep_na %||% TRUE,
                       "VirSorter2 score", paste0(">= ", thr))
      } else {
        note("VirSorter2 score", "no", nrow(tse), nrow(tse), "virsorter2_max_score absent")
      }
      if (!is.null(f$votes)) {
        thr <- num_input(input$min_votes, 1)
        if (thr > 1) {
          numeric_filter("viral_vote_n", thr, FALSE, "Viral-evidence votes",
                         paste0(">= ", thr, " detectors"))
        } else {
          note("Viral-evidence votes", "no", nrow(tse), nrow(tse), "threshold is 1")
        }
      } else {
        note("Viral-evidence votes", "no", nrow(tse), nrow(tse),
             "viral_vote_n absent; build the object with annotate_viral_votes()")
      }

      # Coverage masking
      cf <- vpf_covfrac_assay(tse)
      if (!is.null(cf)) {
        thr <- num_input(input$min_covfrac, 0.75)
        masked <- tryCatch(
          suppressWarnings(refind_abundance(tse, assay_name, cf, thr)),
          error = function(e) e
        )
        if (inherits(masked, "error")) {
          note("Coverage masking", "no", nrow(tse), nrow(tse),
               paste0("failed: ", conditionMessage(masked)))
        } else {
          tse <- masked
          note("Coverage masking", "yes", nrow(tse), nrow(tse),
               paste0("abundance zeroed where covered fraction < ", thr))
        }
      } else {
        note("Coverage masking", "no", nrow(tse), nrow(tse), "no covered-fraction assay")
      }

      # Prevalence and total abundance, both computed after masking
      if (nrow(tse) > 0 && ncol(tse) > 0) {
        prev_thr <- num_input(input$min_prevalence, 0)
        if (prev_thr > 0) {
          before <- nrow(tse)
          prev <- vpf_prevalence(tse, assay_name, detection = 0)
          tse <- vpf_subset_rows(tse, prev >= prev_thr, na_keeps = FALSE)
          note("Prevalence", "yes", before, nrow(tse),
               paste0("detected in >= ", round(prev_thr * 100), "% of samples"))
        } else {
          note("Prevalence", "no", nrow(tse), nrow(tse), "threshold is 0")
        }
      }
      if (nrow(tse) > 0 && ncol(tse) > 0) {
        ab_thr <- num_input(input$min_total_abundance, 0)
        if (ab_thr > 0) {
          before <- nrow(tse)
          totals <- rowSums(SummarizedExperiment::assay(tse, assay_name), na.rm = TRUE)
          tse <- vpf_subset_rows(tse, totals >= ab_thr, na_keeps = FALSE)
          note("Total abundance", "yes", before, nrow(tse), paste0(">= ", ab_thr))
        } else {
          note("Total abundance", "no", nrow(tse), nrow(tse), "threshold is 0")
        }
      }

      list(tse = tse, audit = do.call(rbind, audit))
    })

    filtered_tse <- reactive(filtered()$tse)
    audit_table <- reactive(filtered()$audit)

    output$summary_box <- renderUI({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      keep <- filtered_tse()
      if (is.null(keep)) {
        return(NULL)
      }
      n_in <- nrow(tse)
      n_out <- nrow(keep)
      assay_name <- active_assay()
      empty_samples <- character(0)
      if (n_out > 0 && !is.null(assay_name)) {
        libs <- colSums(SummarizedExperiment::assay(keep, assay_name), na.rm = TRUE)
        empty_samples <- vpf_sample_names(keep)[!is.finite(libs) | libs <= 0]
      }
      tagList(
        vpf_card(
          "Current selection",
          div(
            vpf_stat("Contigs kept", format(n_out, big.mark = ",")),
            vpf_stat("Contigs removed", format(n_in - n_out, big.mark = ",")),
            vpf_stat("Percent kept", if (n_in > 0) sprintf("%.1f%%", 100 * n_out / n_in) else "-"),
            vpf_stat("Samples", ncol(keep)),
            vpf_stat("Assay", assay_name %||% "-")
          )
        ),
        if (n_out == 0) {
          vpf_notice(
            type = "danger",
            title = "No contigs pass the current filters",
            "Every downstream tab will be empty. Loosen the thresholds on the left,",
            "or press Reset all filters."
          )
        },
        if (length(empty_samples) > 0) {
          vpf_notice(
            type = "warning",
            title = "Samples with no signal left",
            paste0("After filtering and coverage masking these samples contain only zeros: ",
                   paste(empty_samples, collapse = ", "),
                   ". Diversity and ordination cannot use them.")
          )
        }
      )
    })

    output$tbl_audit <- reactable::renderReactable({
      df <- audit_table()
      if (is.null(df)) {
        return(reactable::reactable(data.frame(Message = "Load a dataset first.")))
      }
      reactable::reactable(
        df, striped = TRUE, highlight = TRUE, bordered = TRUE, pagination = FALSE,
        columns = list(
          Step = reactable::colDef(minWidth = 150),
          Applied = reactable::colDef(minWidth = 80),
          Removed = reactable::colDef(minWidth = 90),
          Remaining = reactable::colDef(minWidth = 100),
          Detail = reactable::colDef(minWidth = 260)
        )
      )
    })

    # -- Distribution plots --------------------------------------------------
    status_table <- reactive({
      tse <- raw_tse()
      keep <- filtered_tse()
      if (is.null(tse) || is.null(keep)) {
        return(NULL)
      }
      df <- vpf_row_table(tse)
      df$Status <- ifelse(df$Contig %in% rownames(keep), "kept", "removed")
      df
    })

    dist_plot <- function(column, xlab, log_x = FALSE) {
      df <- status_table()
      if (is.null(df)) {
        return(vpf_message_plot("Load a dataset first."))
      }
      if (!column %in% colnames(df)) {
        return(vpf_message_plot(paste0(
          "This dataset has no '", column, "' column, so the distribution cannot be shown."
        )))
      }
      v <- vpf_as_numeric_column(df[[column]])
      if (is.null(v) || !any(is.finite(v))) {
        return(vpf_message_plot(paste0("'", column, "' holds no usable numeric values.")))
      }
      pdf <- data.frame(value = v, Status = df$Status, stringsAsFactors = FALSE)
      pdf <- pdf[is.finite(pdf$value), , drop = FALSE]
      if (nrow(pdf) == 0) {
        return(vpf_message_plot("No values to plot."))
      }
      bins <- max(5, min(30, length(unique(pdf$value))))
      p <- ggplot2::ggplot(pdf, ggplot2::aes(x = .data$value, fill = .data$Status)) +
        ggplot2::geom_histogram(bins = bins, alpha = 0.85, position = "stack") +
        ggplot2::scale_fill_manual(values = c(kept = "#4E79A7", removed = "#BAB0AC")) +
        ggplot2::labs(x = xlab, y = "Contigs") +
        ggplot2::theme_bw()
      if (log_x) {
        p <- p + ggplot2::scale_x_log10()
      }
      vpf_ggplotly(p)
    }

    output$plt_length <- plotly::renderPlotly({
      dist_plot("checkv_contig_length", "Contig length (bp, log scale)", log_x = TRUE)
    })

    output$plt_completeness <- plotly::renderPlotly({
      dist_plot("checkv_completeness", "CheckV completeness (%)")
    })

    output$plt_score <- plotly::renderPlotly({
      df <- status_table()
      if (is.null(df)) {
        return(vpf_message_plot("Load a dataset first."))
      }
      for (col in c("genomad_score", "virsorter2_max_score", "dvf_score")) {
        if (col %in% colnames(df)) {
          v <- vpf_as_numeric_column(df[[col]])
          if (!is.null(v) && any(is.finite(v))) {
            return(dist_plot(col, paste0(col, " (virus score)")))
          }
        }
      }
      vpf_message_plot(paste(
        "No virus score is stored in this dataset.",
        "geNomad and VirSorter2 scores are produced by the detection stage of",
        "the pipeline."
      ))
    })

    output$plt_quality <- plotly::renderPlotly({
      df <- status_table()
      if (is.null(df)) {
        return(vpf_message_plot("Load a dataset first."))
      }
      if (!"checkv_quality" %in% colnames(df)) {
        return(vpf_message_plot(
          "No checkv_quality column: CheckV quality tiers are not available."
        ))
      }
      q <- vpf_blank_na(df$checkv_quality)
      q[is.na(q)] <- "not annotated"
      pdf <- as.data.frame(table(Quality = q, Status = df$Status), stringsAsFactors = FALSE)
      known <- c("Complete", "High-quality", "Medium-quality", "Low-quality",
                 "Not-determined", "not annotated")
      lv <- c(intersect(known, unique(pdf$Quality)), setdiff(unique(pdf$Quality), known))
      pdf$Quality <- factor(pdf$Quality, levels = lv)
      p <- ggplot2::ggplot(pdf, ggplot2::aes(x = .data$Quality, y = .data$Freq,
                                             fill = .data$Status)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c(kept = "#4E79A7", removed = "#BAB0AC")) +
        ggplot2::labs(x = NULL, y = "Contigs") +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
      vpf_ggplotly(p)
    })

    output$evidence_status <- renderUI({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(vpf_notice(type = "info", "Load a dataset on the Data tab first."))
      }
      sel <- vpf_viral_selection(tse)
      if (is.null(vpf_viral_votes(tse))) {
        return(vpf_missing_annotation(
          "Viral-evidence votes",
          "annotate_viral_votes(), applied when the result object is built",
          paste("This records which detector marked each contig as viral and how",
                "many did, so a contig admitted on a single call can be told apart",
                "from one every detector agreed on.")
        ))
      }
      msg <- if (!is.null(sel)) {
        paste0("Selection rule: ", sel$rule %||% "vote",
               if (!is.null(sel$combination)) paste0(" - ", sel$combination) else "",
               ". Votes evaluated: ",
               paste(sel$votes_used %||% "unknown", collapse = ", "),
               if (length(sel$votes_missing %||% character(0)) > 0) {
                 paste0(". Votes that could not be evaluated: ",
                        paste(sel$votes_missing, collapse = ", "))
               } else {
                 ""
               }, ".")
      } else {
        "Vote columns are present but the selection metadata is not."
      }
      vpf_notice(type = "info", title = "How this viral subset was chosen", msg)
    })

    output$plt_votes <- plotly::renderPlotly({
      tse <- raw_tse()
      keep <- filtered_tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset on the Data tab first."))
      }
      votes <- vpf_viral_votes(tse)
      if (is.null(votes)) {
        return(vpf_message_plot(paste(
          "This object does not record viral-evidence votes. They are added by",
          "annotate_viral_votes() when the result object is built."
        )))
      }
      status <- ifelse(rownames(tse) %in% rownames(keep %||% tse), "kept", "removed")
      df <- as.data.frame(table(Votes = votes, Status = status),
                          stringsAsFactors = FALSE)
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Votes, y = .data$Freq, fill = .data$Status,
        text = paste0(.data$Votes, " vote(s), ", .data$Status, ": ", .data$Freq)
      )) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c(kept = "#4E79A7", removed = "#BAB0AC")) +
        ggplot2::labs(x = "Independent detectors calling the contig viral",
                      y = "Contigs") +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$plt_vote_sources <- plotly::renderPlotly({
      tse <- raw_tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset on the Data tab first."))
      }
      rd <- SummarizedExperiment::rowData(tse)
      cols <- grep("^viral_vote_", colnames(rd), value = TRUE)
      cols <- setdiff(cols, c("viral_vote_n", "viral_vote_evidence"))
      if (length(cols) == 0) {
        return(vpf_message_plot(paste(
          "No per-detector vote columns. They are added by",
          "annotate_viral_votes() when the result object is built."
        )))
      }
      counts <- vapply(cols, function(nm) sum(as.logical(rd[[nm]]), na.rm = TRUE),
                       numeric(1))
      df <- data.frame(
        Detector = sub("^viral_vote_", "", cols),
        Contigs = as.numeric(counts),
        stringsAsFactors = FALSE
      )
      df <- df[order(-df$Contigs), , drop = FALSE]
      df$Detector <- factor(df$Detector, levels = rev(df$Detector))
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Detector, y = .data$Contigs,
        text = paste0(.data$Detector, ": ", .data$Contigs, " contigs")
      )) +
        ggplot2::geom_col(fill = "#59A14F") +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = paste0("Contigs called viral (of ", nrow(tse), ")")) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    output$plt_prevalence <- plotly::renderPlotly({
      keep <- filtered_tse()
      assay_name <- active_assay()
      if (is.null(keep) || is.null(assay_name)) {
        return(vpf_message_plot("Load a dataset first."))
      }
      if (nrow(keep) == 0) {
        return(vpf_message_plot("No contigs pass the current filters."))
      }
      prev <- vpf_prevalence(keep, assay_name, detection = 0)
      mat <- SummarizedExperiment::assay(keep, assay_name)
      mean_ab <- rowMeans(mat, na.rm = TRUE)
      pdf <- data.frame(
        Contig = rownames(keep) %||% seq_len(nrow(keep)),
        Prevalence = prev,
        MeanAbundance = mean_ab,
        stringsAsFactors = FALSE
      )
      pdf$MeanAbundance[!is.finite(pdf$MeanAbundance)] <- 0
      p <- ggplot2::ggplot(pdf, ggplot2::aes(
        x = .data$Prevalence, y = .data$MeanAbundance + 1,
        text = paste0(.data$Contig,
                      "<br>Prevalence: ", round(.data$Prevalence, 3),
                      "<br>Mean ", assay_name, ": ", signif(.data$MeanAbundance, 4))
      )) +
        ggplot2::geom_point(alpha = 0.7, colour = "#4E79A7") +
        ggplot2::scale_y_log10() +
        ggplot2::labs(x = "Prevalence (fraction of samples detected)",
                      y = paste0("Mean ", assay_name, " + 1 (log scale)")) +
        ggplot2::theme_bw()
      vpf_ggplotly(p, tooltip = "text")
    })

    list(
      tse = filtered_tse,
      raw = raw_tse,
      assay = active_assay,
      audit = audit_table,
      name = r_input$name
    )
  })
}

## To be copied in the UI
# mod_vpfilter_ui("vpfilter_1")

## To be copied in the server
# mod_vpfilter_server("vpfilter_1", r_input)
