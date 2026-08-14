#' about UI Function
#'
#' @description Tutorial and background, read from the packaged
#'   `inst/app/www/tutorial.md`. The file is located with [app_sys()] because
#'   `inst/` is stripped at install time, so a relative path to it resolves
#'   only when the app runs from a source checkout.
#'
#' @param id Module id.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_about_ui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        9,
        uiOutput(ns("tutorial")),
        uiOutput(ns("footer"))
      ),
      column(
        3,
        vpf_card(
          "Session",
          verbatimTextOutput(ns("versions"))
        ),
        vpf_card(
          "Links",
          tags$ul(
            tags$li(tags$a(href = "https://github.com/deng-lab/viroprofiler",
                           target = "_blank", rel = "noopener", "ViroProfiler pipeline")),
            tags$li(tags$a(href = "https://github.com/deng-lab/vpfkit",
                           target = "_blank", rel = "noopener", "vpfkit source code")),
            tags$li(tags$a(href = "https://github.com/deng-lab/vpfkit/issues",
                           target = "_blank", rel = "noopener", "Report an issue"))
          )
        )
      )
    )
  )
}

#' about Server Functions
#'
#' @param id Module id.
#'
#' @noRd
mod_about_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    output$tutorial <- renderUI({
      path <- app_sys("app/www/tutorial.md")
      if (!nzchar(path) || !file.exists(path)) {
        # Running from a source checkout before the package is installed.
        local <- file.path("inst", "app", "www", "tutorial.md")
        path <- if (file.exists(local)) local else ""
      }
      if (!nzchar(path)) {
        return(vpf_notice(type = "warning", "The tutorial file was not found in this installation."))
      }
      vpf_render_markdown(path)
    })

    output$footer <- renderUI({
      path <- app_sys("app/www/footer.html")
      if (!nzchar(path) || !file.exists(path)) {
        local <- file.path("inst", "app", "www", "footer.html")
        path <- if (file.exists(local)) local else ""
      }
      if (!nzchar(path)) {
        return(NULL)
      }
      tagList(tags$hr(), div(class = "footer", includeHTML(path)))
    })

    output$versions <- renderText({
      pkgs <- c("vpfkit", "shiny", "golem", "mia", "miaViz", "scater", "vegan",
                "plotly", "reactable", "SummarizedExperiment",
                "TreeSummarizedExperiment")
      lines <- vapply(pkgs, function(p) {
        v <- tryCatch(as.character(utils::packageVersion(p)),
                      error = function(e) "not installed")
        sprintf("%-26s %s", p, v)
      }, character(1))
      paste(c(R.version.string, "", lines), collapse = "\n")
    })
  })
}

#' Render a markdown file to HTML without assuming a particular backend
#'
#' `shiny::includeMarkdown()` needs the `markdown` package, which is only a
#' soft dependency here. The plain-text fallback keeps the tab usable rather
#' than raising an error inside `renderUI()`.
#'
#' @param path Path to a markdown file.
#' @return An HTML tag.
#' @noRd
vpf_render_markdown <- function(path) {
  if (requireNamespace("markdown", quietly = TRUE)) {
    out <- tryCatch(includeMarkdown(path), error = function(e) NULL)
    if (!is.null(out)) {
      return(out)
    }
  }
  if (requireNamespace("commonmark", quietly = TRUE)) {
    txt <- paste(readLines(path, warn = FALSE), collapse = "\n")
    return(HTML(commonmark::markdown_html(txt, extensions = TRUE)))
  }
  tags$pre(paste(readLines(path, warn = FALSE), collapse = "\n"))
}

## To be copied in the UI
# mod_about_ui("about_1")

## To be copied in the server
# mod_about_server("about_1")
