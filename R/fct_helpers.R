#' Read taxonomy table
#'
#' @description A fct function
#'
#' @param fin taxonomy annotation file
#'
#' @return DataFrame
#' @export
#'
#' @examples
#' @noRd
vpf_merge <- function(fin) {
  # todo
  df <- fread(fin) %>%
    dplyr::mutate(Contig=str_replace(.data$contig_id, "-cat_[1-6]", "")) %>%
    dplyr::mutate(virsorter_category=as.integer(str_replace(.data$contig_id, ".*-cat_", ""))) %>%
    setnames("length", "viral_length")
  return(df)
}
