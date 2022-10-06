library(tidyverse)
library(data.table)

read_CheckV <- function(fpath) {
  dtbl <- fread(fpath) %>%
    mutate(checkv_quality = factor(checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>%
    column_to_rownames("contig_id") %>%
    setnames(colnames(.), paste0("checkv_", colnames(.))) %>%
    setnames("checkv_checkv_quality", "checkv_quality") %>%
    rownames_to_column("Contig")
  return(dtbl)
}


read_taxonomy <- function(fpath) {
  df <- fread(fpath) %>% 
    mutate(Contig=str_replace(contig_id, "-cat_[1-6]", "")) %>% 
    mutate(virsorter_category=as.integer(str_replace(contig_id, ".*-cat_", "")))
    # setnames("length", "viral_length")
  return(df)
}

read_VirSorter2 <- function(fpath) {
  dtbl <- fread(fpath) %>%
    column_to_rownames("seqname") %>%
    setnames(colnames(.), paste0("virsorter2_", colnames(.))) %>%
    rownames_to_column("Contig")
  return(dtbl)
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


read_coverm <- function(fpath_abundance, sample_metadata) {
  df_abundance <- fread(fpath_abundance) %>%
    column_to_rownames("Contig")
  
  # remove samples in sample_metadata but not in df_abundance
  smeta_selected <- sample_metadata %>%
    dplyr::filter(sample_id %in% colnames(df_abundance))

  # ensure the order of samples in df_abundance is the same as in sample_metadata
  df2 <- as.matrix(df_abundance[, smeta_selected$sample_id])

  return(df2)
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


plot_beta_diversity <- function(tse_obj, name, NMDS=FALSE, scale = F, ...) {
  if (NMDS == TRUE) {
    tse_obj <- runNMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    xlab <- "Axis 1"
    ylab <- "Axis 2"
  } else {
    tse_obj <- runMDS(tse_obj, FUN = vegan::vegdist, name = name, ...)
    e <- attr(reducedDim(tse_obj, name), "eig")
    rel_eig <- 100 * e/sum(e[e>0])
    xlab <- paste("Axis 1 (", round(rel_eig[[1]], 2), "%)", sep = "")
    ylab <- paste("Axis 2 (", round(rel_eig[[2]], 2), "%)", sep = "")
  }
  p <- plotReducedDim(tse_obj, name, colour_by = "condition") +
    xlab(xlab) +
    ylab(ylab) +
    stat_ellipse(geom = "polygon", alpha=0.1, aes(fill=colour_by), linetype=2) +
    ggtitle(name) +
    theme_bw()
  return(p)
}
