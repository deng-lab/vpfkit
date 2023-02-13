#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # Your application server logic
  mod_vpfilter_server("vpfilter_1")
  # mod_vpfvtest_server("vpfvtest_1")
}
