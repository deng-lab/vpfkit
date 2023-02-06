library(tidyverse)
library(data.table)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)

read_taxonomy <- function(fpath) {
  df <- fread(fpath) %>%
    mutate(Contig=str_replace(contig_id, "-cat_[1-6]", "")) %>%
    mutate(virsorter_category=as.integer(str_replace(contig_id, ".*-cat_", "")))
    # setnames("length", "viral_length")
  return(df)
}

read_bins <- function(fbin, fbin_classify) {
  bin_classify <- fread(fbin_classify) %>%
    column_to_rownames("binname") %>%
    setnames(colnames(.), paste0("bin_", colnames(.))) %>%
    rownames_to_column("binname")
  bins <- fread(fbin, header = FALSE, col.names = c("binname", "Contig")) %>%
    left_join(bin_classify, by = "binname")
  return(bins)
}


refind_abundance <- function(tse, raw_abundance_name, covfrac_name, thr=0.5) {
  df_abundance_raw <- assay(tse, raw_abundance_name)
  df_filter <- assay(tse, covfrac_name)
  df_filter[df_filter<thr] <- 0
  df_filter[df_filter!=0] <- 1
  df_abundance_new <- df_abundance_raw * df_filter

  tse_new <- tse
  assay(tse_new, raw_abundance_name) <- df_abundance_new
  return(tse_new)
}
