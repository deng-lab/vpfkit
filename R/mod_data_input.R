#' data_input UI Function
#'
#' @description Loads a ViroProfiler result object, optionally joins an
#'   external sample-metadata table onto `colData`, and reports exactly what
#'   the object contains. Everything downstream reads the object this module
#'   publishes.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_data_input_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        4,
        vpf_card(
          "1. Load a dataset",
          radioButtons(
            ns("source"), NULL,
            choices = c("Bundled example" = "demo",
                        "Upload an .rds file" = "upload",
                        "Path on this server" = "path"),
            selected = "demo"
          ),
          conditionalPanel(
            "input.source == 'demo'", ns = ns,
            selectInput(ns("demo_choice"), "Example dataset", choices = character(0))
          ),
          conditionalPanel(
            "input.source == 'upload'", ns = ns,
            fileInput(ns("upload"), "ViroProfiler output (.rds)",
                      multiple = FALSE, accept = c(".rds", ".RDS")),
            uiOutput(ns("upload_limit"))
          ),
          conditionalPanel(
            "input.source == 'path'", ns = ns,
            uiOutput(ns("path_input"))
          ),
          actionButton(ns("load"), "Load dataset", class = "btn-primary"),
          uiOutput(ns("load_status"))
        ),
        vpf_card(
          "2. Add sample metadata (optional)",
          vpf_caption(
            "ViroProfiler's samplesheet holds only sample, fastq_1 and fastq_2, so",
            "study variables never travel with the pipeline output. Upload them here",
            "to unlock every grouped analysis."
          ),
          fileInput(ns("metadata"), "Metadata table (.csv, .tsv, .txt, .xlsx)",
                    multiple = FALSE,
                    accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
          vpf_caption(
            "The first column, or a column named sample / sample_name / sample_id,",
            "must hold the sample identifiers. All other columns are added to colData."
          ),
          actionButton(ns("apply_metadata"), "Join metadata"),
          actionButton(ns("clear_metadata"), "Remove joined metadata"),
          uiOutput(ns("metadata_status"))
        )
      ),
      column(
        8,
        uiOutput(ns("overview")),
        uiOutput(ns("assay_help")),
        tabsetPanel(
          id = ns("tabs"),
          tabPanel(
            "Annotations available",
            vpf_caption(
              "Which ViroProfiler annotation families reached this object, and",
              "which pipeline option produces the ones that did not."
            ),
            reactable::reactableOutput(ns("tbl_annotations"))
          ),
          tabPanel(
            "Assays",
            vpf_caption(
              "What each assay actually measures. Read this before choosing an",
              "abundance metric: the names alone are not self-explanatory."
            ),
            reactable::reactableOutput(ns("tbl_assays"))
          ),
          tabPanel(
            "Sample table (colData)",
            reactable::reactableOutput(ns("tbl_coldata"))
          ),
          tabPanel(
            "Library sizes",
            plotly::plotlyOutput(ns("plt_libsize"), height = "380px"),
            vpf_caption(
              "Total read count per sample. Strongly unequal libraries affect",
              "observed richness and every untransformed comparison."
            )
          )
        )
      )
    )
  )
}

#' data_input Server Functions
#'
#' @param id Module id.
#'
#' @return A list of reactives: `tse` (the loaded object, metadata joined),
#'   `name` (a display label) and `source` (how it was loaded).
#' @noRd
mod_data_input_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    demos <- vpf_demo_datasets()
    base_tse <- reactiveVal(NULL)
    joined_tse <- reactiveVal(NULL)
    dataset_name <- reactiveVal(NULL)
    load_message <- reactiveVal(NULL)
    renamed_assays <- reactiveVal(character(0))
    metadata_message <- reactiveVal(NULL)

    observe({
      if (length(demos) > 0) {
        updateSelectInput(session, "demo_choice", choices = names(demos),
                          selected = names(demos)[[1]])
      }
    })

    output$path_input <- renderUI({
      if (!vpf_server_path_allowed()) {
        return(vpf_notice(
          type = "info",
          "Loading files by server path is disabled in production deployments.",
          "Use the upload option instead."
        ))
      }
      tagList(
        textInput(ns("path"), "Absolute path to a .rds file",
                  placeholder = "/data/run01/viroprofiler_output.rds"),
        vpf_caption("The path is resolved on the machine running the app.")
      )
    })

    output$upload_limit <- renderUI({
      limit <- getOption("shiny.maxRequestSize", 5 * 1024^2)
      vpf_caption(sprintf(
        "Maximum upload size: %.0f MB. Raise it with run_app(max_upload_mb = ) or the VPFKIT_MAX_UPLOAD_MB environment variable.",
        limit / 1024^2
      ))
    })

    # Loading ---------------------------------------------------------------
    observeEvent(input$load, {
      metadata_message(NULL)
      src <- input$source %||% "demo"
      spec <- switch(
        src,
        demo = {
          choice <- input$demo_choice
          if (is.null(choice) || !choice %in% names(demos)) {
            list(path = NULL, label = NULL,
                 error = "No bundled example dataset is available in this installation.")
          } else {
            list(path = demos[[choice]], label = choice, error = NULL)
          }
        },
        upload = {
          up <- input$upload
          if (is.null(up)) {
            list(path = NULL, label = NULL, error = "Choose a file to upload first.")
          } else {
            list(path = up$datapath, label = up$name, error = NULL)
          }
        },
        path = {
          if (!vpf_server_path_allowed()) {
            list(path = NULL, label = NULL,
                 error = "Loading by server path is disabled in this deployment.")
          } else {
            p <- input$path
            if (is.null(p) || !nzchar(trimws(p))) {
              list(path = NULL, label = NULL, error = "Enter a file path first.")
            } else {
              list(path = trimws(p), label = basename(trimws(p)), error = NULL)
            }
          }
        },
        list(path = NULL, label = NULL, error = "Unknown data source.")
      )

      if (!is.null(spec$error)) {
        load_message(list(type = "danger", text = spec$error))
        return()
      }

      withProgress(message = "Reading dataset", value = 0.4, {
        res <- vpf_read_tse(spec$path)
      })
      if (!is.null(res$error)) {
        load_message(list(type = "danger", text = res$error))
        return()
      }
      base_tse(res$tse)
      joined_tse(NULL)
      dataset_name(spec$label)
      renamed_assays(res$renamed_assays %||% character(0))
      load_message(list(
        type = "success",
        text = sprintf("Loaded %d contigs x %d samples from %s.",
                       nrow(res$tse), ncol(res$tse), spec$label)
      ))
    })

    output$load_status <- renderUI({
      msg <- load_message()
      if (is.null(msg)) {
        return(vpf_notice(
          type = "info",
          "Pick a source and press Load dataset. The bundled examples need no files."
        ))
      }
      vpf_notice(type = msg$type, msg$text)
    })

    # Sample metadata -------------------------------------------------------
    observeEvent(input$apply_metadata, {
      tse <- base_tse()
      if (is.null(tse)) {
        metadata_message(list(type = "danger", text = "Load a dataset first."))
        return()
      }
      up <- input$metadata
      if (is.null(up)) {
        metadata_message(list(type = "danger", text = "Choose a metadata file first."))
        return()
      }
      parsed <- vpf_read_metadata_file(up$datapath, up$name)
      if (!is.null(parsed$error)) {
        metadata_message(list(type = "danger", text = parsed$error))
        return()
      }
      res <- vpf_join_sample_metadata(tse, parsed$data)
      if (!is.null(res$error)) {
        metadata_message(list(type = "danger", text = res$error))
        return()
      }
      joined_tse(res$tse)
      txt <- sprintf("Joined metadata for %d of %d samples. Added columns: %s.",
                     res$matched, ncol(tse),
                     paste(setdiff(colnames(SummarizedExperiment::colData(res$tse)),
                                   colnames(SummarizedExperiment::colData(tse))),
                           collapse = ", "))
      if (length(res$unmatched_samples) > 0) {
        txt <- paste0(txt, " No metadata row matched: ",
                      paste(res$unmatched_samples, collapse = ", "), ".")
      }
      if (length(res$unused_rows) > 0) {
        txt <- paste0(txt, " Metadata rows with no matching sample: ",
                      paste(utils::head(res$unused_rows, 8), collapse = ", "),
                      if (length(res$unused_rows) > 8) ", ..." else "", ".")
      }
      metadata_message(list(
        type = if (length(res$unmatched_samples) > 0) "warning" else "success",
        text = txt
      ))
    })

    observeEvent(input$clear_metadata, {
      joined_tse(NULL)
      metadata_message(list(type = "info", text = "Joined metadata removed."))
    }, ignoreInit = TRUE)

    output$metadata_status <- renderUI({
      msg <- metadata_message()
      if (is.null(msg)) {
        return(NULL)
      }
      vpf_notice(type = msg$type, msg$text)
    })

    active_tse <- reactive({
      joined_tse() %||% base_tse()
    })

    # Overview --------------------------------------------------------------
    output$overview <- renderUI({
      tse <- active_tse()
      if (is.null(tse)) {
        return(vpf_notice(
          type = "info",
          title = "No dataset loaded",
          "Load a bundled example on the left to explore the viewer without any files,",
          "or upload the viroprofiler_output.rds produced by your own run."
        ))
      }
      groups <- vpf_group_candidates(tse)
      ranks <- vpf_available_ranks(tse)
      tagList(
        vpf_card(
          paste0("Dataset: ", dataset_name() %||% "unnamed"),
          div(
            vpf_stat("Contigs", format(nrow(tse), big.mark = ",")),
            vpf_stat("Samples", ncol(tse)),
            vpf_stat("Assays", length(SummarizedExperiment::assayNames(tse))),
            vpf_stat("Annotation columns", ncol(SummarizedExperiment::rowData(tse))),
            vpf_stat("Sample variables", ncol(SummarizedExperiment::colData(tse))),
            vpf_stat("Taxonomic ranks", length(ranks))
          )
        ),
        if (length(groups) == 0) {
          vpf_notice(
            type = "warning",
            title = "No grouping variable",
            "colData holds no column that can separate samples into groups, so",
            "group comparisons, PERMANOVA and differential analysis stay disabled.",
            "Upload a sample-metadata table on the left to supply one.",
            "A sample-identifier column does not count: it creates one group per sample."
          )
        } else {
          vpf_notice(
            type = "success",
            paste0("Grouping variables available: ", paste(groups, collapse = ", "), ".")
          )
        },
        if (ncol(tse) < 3) {
          vpf_notice(
            type = "warning",
            title = paste0("Only ", ncol(tse), " sample(s)"),
            "Per-sample diversity and direct pairwise comparison remain valid.",
            "Ordination, group tests and PERMANOVA need more samples and are",
            "disabled where they would be meaningless."
          )
        }
      )
    })

    output$assay_help <- renderUI({
      tse <- active_tse()
      if (is.null(tse)) {
        return(NULL)
      }
      notices <- list()
      renamed <- renamed_assays()
      if (length(renamed) > 0) {
        notices <- c(notices, list(vpf_notice(
          type = "warning",
          title = "A legacy assay name was corrected while reading this object",
          paste0(paste(sprintf("'%s' is now '%s'", names(renamed), unname(renamed)),
                       collapse = "; "), ". "),
          "CoverM's trimmed mean is a mean per-base coverage DEPTH in x-fold, not",
          "edgeR's TMM (trimmed mean of M-values) normalization; the old name",
          "implied that it was, and treating it as a normalized abundance produces",
          "wrong composition, diversity and ordination results without any error.",
          "The file on disk is unchanged."
        )))
      }
      info <- vpf_assay_info(tse)
      if ("tmm" %in% info$assay) {
        notices <- c(notices, list(vpf_notice(
          type = "danger",
          title = "This object stores CoverM trimmed mean under the name 'tmm'",
          "It holds a mean per-base coverage DEPTH in x-fold, not edgeR's TMM",
          "(trimmed mean of M-values) normalization. Newer ViroProfiler runs name",
          "this assay", shiny::tags$code("trimmed_mean"), "."
        )))
      }
      if (length(notices) == 0) NULL else do.call(tagList, notices)
    })

    output$tbl_annotations <- reactable::renderReactable({
      tse <- active_tse()
      if (is.null(tse)) {
        return(reactable::reactable(data.frame(Message = "Load a dataset first.")))
      }
      st <- vpf_annotation_status(tse)
      st$status <- ifelse(st$available, "present", "absent")
      out <- st[, c("family", "status", "n_columns", "param", "columns", "note")]
      colnames(out) <- c("Annotation", "Status", "Columns", "ViroProfiler option",
                         "Column names", "What it contains")
      reactable::reactable(
        out, striped = TRUE, highlight = TRUE, bordered = TRUE, wrap = TRUE,
        defaultPageSize = 15,
        columns = list(
          Annotation = reactable::colDef(minWidth = 140),
          Status = reactable::colDef(
            minWidth = 80,
            style = function(value) {
              if (identical(value, "present")) {
                list(color = "#188f4c", fontWeight = "bold")
              } else {
                list(color = "#b94a48")
              }
            }
          ),
          Columns = reactable::colDef(minWidth = 80),
          `ViroProfiler option` = reactable::colDef(minWidth = 220),
          `Column names` = reactable::colDef(minWidth = 220),
          `What it contains` = reactable::colDef(minWidth = 240)
        )
      )
    })

    output$tbl_assays <- reactable::renderReactable({
      tse <- active_tse()
      if (is.null(tse)) {
        return(reactable::reactable(data.frame(Message = "Load a dataset first.")))
      }
      info <- vpf_assay_info(tse)
      info$usable <- ifelse(info$role == "abundance",
                            "composition, diversity, ordination",
                            "presence and detection only")
      out <- info[, c("assay", "label", "unit", "usable", "description")]
      colnames(out) <- c("Assay", "Measures", "Unit", "Valid for", "Details")
      reactable::reactable(
        out, striped = TRUE, highlight = TRUE, bordered = TRUE, wrap = TRUE,
        columns = list(
          Assay = reactable::colDef(minWidth = 100),
          Measures = reactable::colDef(minWidth = 180),
          Unit = reactable::colDef(minWidth = 110),
          `Valid for` = reactable::colDef(minWidth = 180),
          Details = reactable::colDef(minWidth = 380)
        )
      )
    })

    output$tbl_coldata <- reactable::renderReactable({
      tse <- active_tse()
      if (is.null(tse)) {
        return(reactable::reactable(data.frame(Message = "Load a dataset first.")))
      }
      df <- as.data.frame(SummarizedExperiment::colData(tse), optional = TRUE)
      df <- cbind(data.frame(Sample = colnames(tse), stringsAsFactors = FALSE), df)
      reactable::reactable(vpf_round_df(df), striped = TRUE, highlight = TRUE,
                           bordered = TRUE, resizable = TRUE, wrap = FALSE,
                           defaultPageSize = 15)
    })

    output$plt_libsize <- plotly::renderPlotly({
      tse <- active_tse()
      if (is.null(tse)) {
        return(vpf_message_plot("Load a dataset first."))
      }
      assay_name <- vpf_default_assay(tse)
      if (is.null(assay_name)) {
        return(vpf_message_plot("This object has no abundance assay."))
      }
      mat <- SummarizedExperiment::assay(tse, assay_name)
      df <- data.frame(
        Sample = factor(vpf_sample_names(tse), levels = vpf_sample_names(tse)),
        Total = colSums(mat, na.rm = TRUE),
        Detected = colSums(!is.na(mat) & mat > 0),
        stringsAsFactors = FALSE
      )
      p <- ggplot2::ggplot(df, ggplot2::aes(
        x = .data$Sample, y = .data$Total,
        text = paste0("Sample: ", .data$Sample,
                      "<br>Total: ", format(round(.data$Total, 2), big.mark = ","),
                      "<br>Contigs detected: ", .data$Detected)
      )) +
        ggplot2::geom_col(fill = "#4E79A7") +
        ggplot2::labs(x = NULL, y = paste0("Total ", assay_name)) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      vpf_ggplotly(p, tooltip = "text")
    })

    list(
      tse = active_tse,
      name = reactive(dataset_name()),
      source = reactive(input$source)
    )
  })
}

#' Is loading a dataset by server path permitted?
#'
#' Enabled for local use, disabled once the app runs in production, where an
#' arbitrary path would let any visitor read files from the host.
#'
#' @noRd
vpf_server_path_allowed <- function() {
  getOption("vpfkit.allow_server_path", !isTRUE(getOption("golem.app.prod")))
}

## To be copied in the UI
# mod_data_input_ui("data_input_1")

## To be copied in the server
# mod_data_input_server("data_input_1")
