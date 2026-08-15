#' bigtable UI Function
#'
#' @description A reusable table module that keeps data frames on the server.
#'   Searching, sorting and paging all happen in R, and only the visible page is
#'   sent to the client, so a hundred-thousand-row annotation table never has to
#'   be serialized into the page.
#'
#'   One engine does all three for every table, whatever its size. The
#'   alternative -- letting reactable search and sort small tables in the
#'   browser -- gives the user two sets of controls over one table, and the
#'   download button can only honour one of them: a browser-side filter narrows
#'   the screen while the download still carries what R last selected. Choosing
#'   "All" under "Rows per page" renders a whole small table at once, which is
#'   what the client-side mode was for.
#'
#' @param id Module id.
#' @param height Table height, passed to `reactable`.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_bigtable_ui <- function(id, height = "520px") {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(4, textInput(ns("search"), "Search", placeholder = "Match any column")),
      column(3, selectInput(ns("sort_by"), "Sort by", choices = character(0))),
      column(2, selectInput(ns("sort_dir"), "Order",
                            choices = c("Ascending" = "asc", "Descending" = "desc"))),
      column(3, selectInput(ns("page_size"), "Rows per page",
                            choices = c("10" = "10", "25" = "25", "50" = "50",
                                        "100" = "100", "250" = "250",
                                        "All" = "all"),
                            selected = "25"))
    ),
    uiOutput(ns("column_picker")),
    div(style = "margin-bottom: 6px;", uiOutput(ns("pager"))),
    reactable::reactableOutput(ns("table"), height = height),
    div(
      style = "margin-top: 8px;",
      downloadButton(ns("download_tsv"), "Download this table (TSV)"),
      downloadButton(ns("download_csv"), "Download this table (CSV)")
    ),
    vpf_caption(textOutput(ns("table_note"), inline = TRUE))
  )
}

#' bigtable Server Functions
#'
#' @param id Module id.
#' @param r_data A reactive returning a `data.frame`, or `NULL`.
#' @param filename A reactive or constant giving the download file stem.
#' @param default_columns Optional character vector, or a reactive returning
#'   one, naming the columns shown initially.
#'
#' @return A reactive returning the currently filtered (not paged) data frame.
#' @noRd
mod_bigtable_server <- function(id, r_data, filename = "table",
                                default_columns = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    full_data <- reactive({
      df <- r_data()
      if (is.null(df)) {
        return(NULL)
      }
      df <- as.data.frame(df, stringsAsFactors = FALSE)
      # Row names would render as an extra unlabelled column, duplicating
      # whatever key column the caller already put first.
      rownames(df) <- NULL
      if (ncol(df) == 0) NULL else df
    })

    current_page <- reactiveVal(1L)

    # Column selection ------------------------------------------------------
    output$column_picker <- renderUI({
      df <- full_data()
      if (is.null(df)) {
        return(NULL)
      }
      cand <- if (is.function(default_columns)) default_columns() else default_columns
      sel <- if (!is.null(cand)) {
        intersect(cand, colnames(df))
      } else {
        colnames(df)
      }
      if (length(sel) == 0) sel <- colnames(df)
      selectInput(ns("columns"), "Columns to display", choices = colnames(df),
                  selected = sel, multiple = TRUE, width = "100%")
    })

    visible_columns <- reactive({
      df <- full_data()
      req(df)
      sel <- input$columns
      sel <- intersect(sel, colnames(df))
      if (length(sel) == 0) colnames(df) else sel
    })

    observeEvent(full_data(), {
      df <- full_data()
      if (is.null(df)) {
        return()
      }
      updateSelectInput(session, "sort_by",
                        choices = c("(none)" = "", colnames(df)),
                        selected = "")
      current_page(1L)
    }, ignoreNULL = FALSE)

    # Filtering and sorting happen here, over the whole table, so that a
    # download describes the same selection as the screen.
    filtered <- reactive({
      df <- full_data()
      if (is.null(df)) {
        return(NULL)
      }
      needle <- input$search
      if (!is.null(needle) && nzchar(trimws(needle))) {
        # Both halves matter, and `grepl` cannot give them at once: `fixed =
        # TRUE` is what stops a contig name like `NODE_1[2]` being read as a
        # regular expression, but it makes R *ignore* `ignore.case` -- silently
        # as far as the result goes, and with a warning per column. Folding both
        # sides instead keeps the match literal and case-insensitive.
        pattern <- tolower(trimws(needle))
        cols <- visible_columns()
        hit <- rep(FALSE, nrow(df))
        for (nm in cols) {
          v <- as.character(df[[nm]])
          v[is.na(v)] <- ""
          hit <- hit | grepl(pattern, tolower(v), fixed = TRUE)
        }
        df <- df[hit, , drop = FALSE]
      }
      sort_by <- input$sort_by
      if (!is.null(sort_by) && nzchar(sort_by) && sort_by %in% colnames(df) && nrow(df) > 0) {
        ord <- order(df[[sort_by]], decreasing = identical(input$sort_dir, "desc"),
                     na.last = TRUE)
        df <- df[ord, , drop = FALSE]
      }
      df
    })

    page_size <- reactive({
      ps <- input$page_size
      if (length(ps) != 1L) {
        return(25L)
      }
      # "All" is a page as tall as the result, which keeps every downstream
      # calculation (n_pages, the row range, the slice) on one code path.
      if (identical(as.character(ps), "all")) {
        df <- filtered()
        return(if (is.null(df) || nrow(df) == 0L) 25L else nrow(df))
      }
      ps <- suppressWarnings(as.integer(ps))
      if (is.na(ps) || ps < 1) 25L else ps
    })

    n_pages <- reactive({
      df <- filtered()
      if (is.null(df) || nrow(df) == 0) {
        return(1L)
      }
      max(1L, as.integer(ceiling(nrow(df) / page_size())))
    })

    observeEvent(list(input$search, input$sort_by, input$sort_dir, input$page_size), {
      current_page(1L)
    }, ignoreInit = TRUE)

    observeEvent(input$prev_page, {
      current_page(max(1L, current_page() - 1L))
    }, ignoreInit = TRUE)

    observeEvent(input$next_page, {
      current_page(min(n_pages(), current_page() + 1L))
    }, ignoreInit = TRUE)

    output$pager <- renderUI({
      df <- filtered()
      if (is.null(df)) {
        return(NULL)
      }
      total <- nrow(df)
      # One page holds everything: the buttons would be inert, so say the size
      # and stop.
      if (n_pages() <= 1L) {
        return(vpf_caption(sprintf("%s rows", format(total, big.mark = ","))))
      }
      pg <- min(current_page(), n_pages())
      from <- if (total == 0) 0 else (pg - 1L) * page_size() + 1L
      to <- min(total, pg * page_size())
      tagList(
        actionButton(ns("prev_page"), "Previous", class = "btn-sm"),
        actionButton(ns("next_page"), "Next", class = "btn-sm"),
        tags$span(
          style = "margin-left: 12px;",
          sprintf("Rows %s-%s of %s  (page %d of %d)",
                  format(from, big.mark = ","), format(to, big.mark = ","),
                  format(total, big.mark = ","), pg, n_pages())
        )
      )
    })

    page_data <- reactive({
      df <- filtered()
      if (is.null(df)) {
        return(NULL)
      }
      df <- df[, visible_columns(), drop = FALSE]
      if (nrow(df) == 0) {
        return(df)
      }
      pg <- min(current_page(), n_pages())
      from <- (pg - 1L) * page_size() + 1L
      to <- min(nrow(df), pg * page_size())
      df[from:to, , drop = FALSE]
    })

    output$table <- reactable::renderReactable({
      df <- page_data()
      if (is.null(df)) {
        return(reactable::reactable(data.frame(Message = "No data loaded yet.")))
      }
      if (nrow(df) == 0) {
        return(reactable::reactable(data.frame(
          Message = "No rows match the current filters."
        )))
      }
      df <- vpf_round_df(df)
      reactable::reactable(
        df,
        # reactable is a renderer here, nothing more. Its own controls act on
        # the rows it was handed, which are one page of a selection R already
        # made, so leaving them on would sort a page and present it as a sort of
        # the table -- and would filter the screen without touching what the
        # download button writes.
        sortable = FALSE,
        searchable = FALSE,
        pagination = FALSE,
        striped = TRUE, highlight = TRUE, bordered = TRUE,
        resizable = TRUE, wrap = FALSE, compact = TRUE,
        defaultColDef = reactable::colDef(minWidth = 110)
      )
    })

    output$table_note <- renderText({
      df <- full_data()
      if (is.null(df)) {
        return("")
      }
      sprintf(paste("Search and sort are applied in R across all %s rows, and",
                    "only the current page is sent to the browser. A download",
                    "contains the whole search result, not the page on screen."),
              format(nrow(df), big.mark = ","))
    })

    stem <- function() {
      if (is.function(filename)) filename() else filename
    }

    # Named separately so `shiny::testServer()` can invoke them: in a mock
    # session `output$id` resolves to the download URL, not to the handler.
    name_tsv <- function() paste0(stem(), "_", Sys.Date(), ".tsv")
    name_csv <- function() paste0(stem(), "_", Sys.Date(), ".csv")

    # The download always contains the full search result rather than the page
    # currently on screen.
    write_tsv <- function(file) {
      df <- filtered()
      if (is.null(df)) df <- data.frame()
      utils::write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE,
                         na = "")
    }
    write_csv <- function(file) {
      df <- filtered()
      if (is.null(df)) df <- data.frame()
      utils::write.csv(df, file, row.names = FALSE, na = "")
    }

    output$download_tsv <- downloadHandler(
      filename = name_tsv, contentType = "text/tab-separated-values",
      content = write_tsv
    )

    output$download_csv <- downloadHandler(
      filename = name_csv, contentType = "text/csv", content = write_csv
    )

    filtered
  })
}

## To be copied in the UI
# mod_bigtable_ui("bigtable_1")

## To be copied in the server
# mod_bigtable_server("bigtable_1", r_data = reactive(NULL))
