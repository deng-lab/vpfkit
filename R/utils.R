library(tidyverse)
library(data.table)

read_CheckV <- function(fpath) {
  df <- fread(fpath) %>% 
    mutate(checkv_quality = factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>% 
    setnames("contig_id", "Contig")
  return(df)
}


read_taxonomy <- function(fpath) {
  df <- fread(fpath) %>% 
    mutate(Contig=str_replace(contig_id, "-cat_[1-6]", "")) %>% 
    mutate(virsorter_category=as.integer(str_replace(contig_id, ".*-cat_", ""))) %>% 
    setnames("length", "viral_length")
  return(df)
}


read_CATBAT <- function(fpath) {
  df <- vroom::vroom(fpath) %>% 
    column_to_rownames("# contig") %>% 
    setnames(old = "lineage scores", new = "lineage_score") %>% 
    setnames(old = colnames(.), new = paste0("CAT_", colnames(.))) %>% 
    rownames_to_column("contig_id")
  return(df)
}


read_vContact2 <- function(fpath, assembler_label="_NODE_") {
  df <- fread(fpath) %>% 
    dplyr::select(-"V1") %>% 
    filter(str_detect(Genome, assembler_label)) %>% 
    column_to_rownames("Genome") %>% 
    setnames(colnames(.), paste0("vc_", str_replace_all(colnames(.), " ", "_"))) %>% 
    rownames_to_column("contig_id")
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




