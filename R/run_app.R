#' Run the Shiny Application
#'
#' Starts ViroProfiler-viewer, an interactive viewer for the
#' `TreeSummarizedExperiment` objects produced by the ViroProfiler pipeline.
#'
#' @param max_upload_mb Largest file the app will accept through its upload
#'   controls, in megabytes. Shiny's own default is 5 MB, which is smaller than
#'   many real result objects. The `VPFKIT_MAX_UPLOAD_MB` environment variable
#'   overrides the default when the argument is left at `NULL`. Note that a
#'   reverse proxy in front of the app can reject a large request before Shiny
#'   ever sees it.
#' @param allow_server_path Whether the app may load a dataset from a path on
#'   the host. Convenient locally, unsafe on a public deployment, so it
#'   defaults to `TRUE` only outside production.
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @return A Shiny application object.
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  max_upload_mb = NULL,
  allow_server_path = NULL,
  ...
) {
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)

  if (is.null(max_upload_mb)) {
    max_upload_mb <- suppressWarnings(
      as.numeric(Sys.getenv("VPFKIT_MAX_UPLOAD_MB", unset = "500"))
    )
  }
  if (!is.numeric(max_upload_mb) || length(max_upload_mb) != 1 ||
      !is.finite(max_upload_mb) || max_upload_mb <= 0) {
    max_upload_mb <- 500
  }
  # Shiny reads this option when it receives the upload, so it has to be set
  # before the app starts rather than inside a reactive. `base::options` is
  # named explicitly because `options` is also one of this function's
  # arguments, where it means the list handed to `shinyApp()`.
  base::options(shiny.maxRequestSize = max_upload_mb * 1024^2)

  if (!is.null(allow_server_path)) {
    base::options(vpfkit.allow_server_path = isTRUE(allow_server_path))
  }

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
