#' export UI Function
#'
#' @description Downloads of the filtered object and its derived tables, plus a
#'   provenance record of exactly which filters produced them. Every download
#'   reflects the current selection, so a figure and a table taken from the same
#'   session always describe the same contigs.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_export_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        6,
        vpf_card(
          "Filtered dataset",
          vpf_caption(
            "A TreeSummarizedExperiment holding only the contigs that pass the",
            "current filters, with the active assay already coverage-masked.",
            "Reload it here or read it with readRDS() in your own analysis."
          ),
          downloadButton(ns("dl_tse"), "Download TSE (.rds)", class = "btn-primary")
        ),
        vpf_card(
          "Tables",
          vpf_caption("Written by the package's own export functions."),
          downloadButton(ns("dl_abundance_csv"), "Abundance matrix (CSV)"),
          downloadButton(ns("dl_abundance_xlsx"), "Abundance matrix (Excel)"),
          tags$br(), tags$br(),
          downloadButton(ns("dl_annotations_tsv"), "Contig annotations (TSV)"),
          downloadButton(ns("dl_annotations_csv"), "Contig annotations (CSV)"),
          tags$br(), tags$br(),
          downloadButton(ns("dl_metadata"), "Sample metadata (CSV)"),
          downloadButton(ns("dl_genes"), "Gene annotations (TSV)")
        ),
        vpf_card(
          "Provenance",
          vpf_caption(
            "A plain-text record of the dataset, the active assay, every filter",
            "step and the package versions used. Keep it with the exported tables."
          ),
          downloadButton(ns("dl_provenance"), "Download provenance (TXT)")
        )
      ),
      column(
        6,
        vpf_card(
          "Quality report",
          vpf_caption(
            "Renders the packaged Quarto template against the filtered object.",
            "Rendering happens on the server and needs a working Quarto",
            "installation, so it can take a minute."
          ),
          uiOutput(ns("report_availability")),
          downloadButton(ns("dl_report"), "Generate and download HTML report")
        ),
        vpf_card(
          "Current selection",
          verbatimTextOutput(ns("selection_summary"))
        )
      )
    )
  )
}

#' export Server Functions
#'
#' @param id Module id.
#' @param r_filter A list of reactives from [mod_vpfilter_server()].
#'
#' @noRd
mod_export_server <- function(id, r_filter) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    stem <- reactive({
      base <- r_filter$name() %||% "viroprofiler"
      base <- tools::file_path_sans_ext(basename(base))
      gsub("[^A-Za-z0-9_.-]+", "_", base)
    })

    require_tse <- function() {
      tse <- r_filter$tse()
      validate(need(!is.null(tse), "Load a dataset on the Data tab first."))
      validate(need(nrow(tse) > 0, "No contigs pass the current filters."))
      tse
    }

    provenance_text <- reactive({
      tse <- r_filter$tse()
      raw <- r_filter$raw()
      lines <- c(
        "ViroProfiler-viewer export provenance",
        paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
        paste0("Dataset: ", r_filter$name() %||% "unnamed"),
        ""
      )
      if (is.null(tse) || is.null(raw)) {
        return(paste(c(lines, "No dataset loaded."), collapse = "\n"))
      }
      lines <- c(
        lines,
        paste0("Contigs before filtering: ", nrow(raw)),
        paste0("Contigs after filtering:  ", nrow(tse)),
        paste0("Samples:                  ", ncol(tse)),
        paste0("Active abundance assay:   ", r_filter$assay() %||% "none"),
        paste0("Assays present:           ",
               paste(SummarizedExperiment::assayNames(tse), collapse = ", ")),
        "",
        "Filter steps",
        "------------"
      )
      audit <- r_filter$audit()
      if (!is.null(audit) && nrow(audit) > 0) {
        # Formatted by hand: print() wraps a wide data frame into stacked
        # blocks, which separates each step from its own detail text.
        lines <- c(lines, sprintf(
          "%-22s %-8s %9s %9s   %s",
          audit$Step, audit$Applied,
          format(audit$Removed, big.mark = ","),
          format(audit$Remaining, big.mark = ","),
          audit$Detail
        ))
      }
      info <- vpf_assay_info(tse)
      lines <- c(
        lines, "", "Assay semantics", "---------------",
        unlist(lapply(seq_len(nrow(info)), function(i) {
          paste0(info$assay[[i]], ": ", info$description[[i]])
        })),
        "", "Session", "-------",
        paste0("R: ", R.version.string),
        unlist(lapply(c("vpfkit", "mia", "miaViz", "SummarizedExperiment",
                        "TreeSummarizedExperiment", "vegan", "shiny"), function(p) {
          v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
          paste0(p, ": ", v)
        }))
      )
      paste(lines, collapse = "\n")
    })

    output$selection_summary <- renderText({
      provenance_text()
    })

    # Each download is a named pair of a file-name function and a content
    # writer. Keeping them out of `downloadHandler()` is what makes them
    # reachable from `shiny::testServer()`, where `output$id` resolves to the
    # download URL rather than to the handler.
    dl_name <- function(suffix, ext) {
      function() paste0(stem(), "_", suffix, "_", Sys.Date(), ".", ext)
    }

    name_tse <- dl_name("filtered", "rds")
    write_tse <- function(file) export_vpftse(require_tse(), file)

    name_abundance_csv <- function() {
      paste0(stem(), "_abundance_", r_filter$assay() %||% "assay", "_", Sys.Date(), ".csv")
    }
    name_abundance_xlsx <- function() {
      paste0(stem(), "_abundance_", r_filter$assay() %||% "assay", "_", Sys.Date(), ".xlsx")
    }
    write_abundance <- function(file, format) {
      tse <- require_tse()
      assay_name <- r_filter$assay()
      validate(need(!is.null(assay_name), "No abundance assay is selected."))
      export_abundance(tse, file, assay.type = assay_name, format = format)
    }

    name_annotations_tsv <- dl_name("annotations", "tsv")
    name_annotations_csv <- dl_name("annotations", "csv")
    write_annotations <- function(file, format) {
      export_annotations(require_tse(), file, format = format)
    }

    name_metadata <- dl_name("sample_metadata", "csv")
    write_metadata <- function(file) {
      tse <- require_tse()
      df <- as.data.frame(SummarizedExperiment::colData(tse), optional = TRUE)
      df <- cbind(data.frame(Sample = colnames(tse), stringsAsFactors = FALSE), df)
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }

    name_genes <- dl_name("gene_annotations", "tsv")
    write_genes <- function(file) {
      tse <- require_tse()
      ga <- vpf_gene_annotations(r_filter$raw())
      validate(need(!is.null(ga), paste(
        "This dataset has no gene-level annotations.",
        "They are produced by the ViroProfiler options --use_dram and pharokka."
      )))
      key <- vpf_gene_contig_column(ga)
      if (!is.null(key)) {
        ga <- ga[as.character(ga[[key]]) %in% rownames(tse), , drop = FALSE]
      }
      utils::write.table(ga, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    }

    name_provenance <- dl_name("provenance", "txt")
    write_provenance <- function(file) writeLines(provenance_text(), file)

    name_report <- dl_name("report", "html")
    write_report <- function(file) {
      tse <- require_tse()
      validate(need(requireNamespace("quarto", quietly = TRUE),
                    "The quarto R package is required to render the report."))
      withProgress(message = "Rendering the report", value = 0.3, {
        generate_report(tse, file)
      })
    }

    output$dl_tse <- downloadHandler(
      filename = name_tse, contentType = "application/octet-stream",
      content = write_tse
    )
    output$dl_abundance_csv <- downloadHandler(
      filename = name_abundance_csv, contentType = "text/csv",
      content = function(file) write_abundance(file, "csv")
    )
    output$dl_abundance_xlsx <- downloadHandler(
      filename = name_abundance_xlsx,
      contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
      content = function(file) write_abundance(file, "xlsx")
    )
    output$dl_annotations_tsv <- downloadHandler(
      filename = name_annotations_tsv, contentType = "text/tab-separated-values",
      content = function(file) write_annotations(file, "tsv")
    )
    output$dl_annotations_csv <- downloadHandler(
      filename = name_annotations_csv, contentType = "text/csv",
      content = function(file) write_annotations(file, "csv")
    )
    output$dl_metadata <- downloadHandler(
      filename = name_metadata, contentType = "text/csv", content = write_metadata
    )
    output$dl_genes <- downloadHandler(
      filename = name_genes, contentType = "text/tab-separated-values",
      content = write_genes
    )
    output$dl_provenance <- downloadHandler(
      filename = name_provenance, contentType = "text/plain",
      content = write_provenance
    )

    output$report_availability <- renderUI({
      if (!.vpf_quarto_available()) {
        return(vpf_notice(
          type = "warning",
          .vpf_quarto_hint()
        ))
      }
      NULL
    })

    output$dl_report <- downloadHandler(
      filename = name_report, contentType = "text/html", content = write_report
    )
  })
}

## To be copied in the UI
# mod_export_ui("export_1")

## To be copied in the server
# mod_export_server("export_1", r_filter)
