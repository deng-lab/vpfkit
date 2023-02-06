#' Read replication cycle prediction results
#'
#' @description A fct function
#'
#' @param fin file path of replication cycle prediction results
#' @param tool tool used to predict replication cycle
#' @param version version of the tool used to predict replication cycle
#'
#' @return data.frame
#' @export
#'
#' @examples
#' @noRd
read_replicyc <- function(fin, tool = "bacphlip", version = "v1") {
  if (tool == "bacphlip") {
    df <- data.table::fread(fin) %>%
      dplyr::mutate(bacphlip_replicyc = ifelse(.data$Virulent > .data$Temperate, "virulent", "temperate"), Contig = .data$V1) %>%
      dplyr::select(c("Contig", "bacphlip_replicyc"))
  }
  return(df)
}


#' Read CATBAT results
#'
#' @param fin CATBAT taxonomy annotation file
#'
#' @return data.frame
#' @export
#'
#' @examples
read_catbat <- function(fin) {
  . <- NULL
  df <- fread(fin, fill = TRUE, sep = "\t") %>%
    tibble::column_to_rownames("# contig") %>%
    setnames(old = "lineage scores", new = "lineage_score") %>%
    setnames(old = colnames(.), new = paste0("CATBAT_", colnames(.))) %>%
    tibble::rownames_to_column("contig_id")
  return(df)
}


#' Read CheckV results
#'
#' Add "checkv_" prefix to column names.
#'
#' @param fin CheckV results
#'
#' @return data.frame
#' @export
#'
#' @examples
read_checkv <- function(fin) {
  . <- NULL
  df <- fread(fin) %>%
    dplyr::mutate(checkv_quality = factor(.data$checkv_quality, levels=c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined"))) %>%
    column_to_rownames("contig_id") %>%
    setnames(colnames(.), paste0("checkv_", colnames(.))) %>%
    setnames("checkv_checkv_quality", "checkv_quality") %>%
    rownames_to_column("Contig")
  return(df)
}


#' Read CoverM output file
#'
#' @param fpath CoverM output file name
#' @param fbin2contig Contig binning to contig mapping
#' @param cov Coverage fraction file or not
#'
#' @return A dataframe
#' @export
#'
#' @examples
read_coverm <- function(fpath, fbin2contig = 0, cov = 0) {
  . <- NULL
  df_abundance <- fread(fpath) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "ds10Ms", "Sample_")) %>%
    setnames(colnames(.), str_replace_all(colnames(.), "Sample_0", "Sample_")) %>%
    setnames("Contig", "genome_id")

  # check if bin2contig file is provided
  if (fbin2contig!=0) {
    df_bin2contig <- fread(fbin2contig, header = FALSE, col.names = c("bin_id", "genome_id"))
    df_abundance <- df_abundance %>%
      dplyr::left_join(df_bin2contig, by = "genome_id") %>%
      dplyr::mutate(bin_id = ifelse(is.na(.data$bin_id), .data$genome_id, .data$bin_id)) %>%
      dplyr::mutate(genome_id = .data$bin_id) %>%
      dplyr::select(-.data$bin_id) %>%
      dplyr::group_by(.data$genome_id)

    # If read coverage fraction, then choose the maximum value for each bin
    if (cov!=0) {
      df_abundance <- df_abundance %>% dplyr::summarise_all(max)
    } else {
      df_abundance <- df_abundance %>% dplyr::summarise_all(sum)
    }
  }

  df_abundance <- df_abundance %>%
    column_to_rownames("genome_id")

  return(df_abundance)
}


#' Read DeepVirFinder results
#'
#' @param fin file of DeepVirFinder results
#' @param thr_score threshold of score
#' @param thr_pvalue threshold of pvalue
#' @param thr_qvalue threshold of qvalue
#'
#' @return data.frame
#' @export
#'
#' @examples
read_dvf <- function(fin, thr_score=0.9, thr_pvalue=0.01, thr_qvalue=0.01) {
  df <- fread(fin) %>%
    dplyr::filter(.data$score > thr_score) %>%
    dplyr::filter(.data$pvalue < thr_pvalue) %>%
    dplyr::filter(.data$qvalue < thr_qvalue) %>%
    dplyr::mutate(dvf_score = .data$score, Contig = .data$name) %>%
    dplyr::select(c("Contig", "dvf_score"))

  return(df)
}


#' Read iphop results
#'
#' @param fin iphop genus prediction
#'
#' @return data.frame
#' @export
#'
#' @examples
read_iphop <- function(fin) {
  . <- NULL
  df <- fread(fin) %>%
    setnames(colnames(.), c("Contig", "iphop_aai2ref", "iphop_genus", "iphop_score", "iphop_methods"))

  return(df)
}


#' Read replication cycle prediction results from Bacphlip
#'
#' @param fin file path of Bacphlip
#' @param tool tool name
#' @param version tool version
#'
#' @return data.frame
#' @export
#'
#' @examples
read_replicyc <- function(fin, tool = "bacphlip", version = "v1") {
  df <- fread(fin)

  if (tool == "bacphlip") {
    df <- df %>%
      dplyr::mutate(bacphlip_replicyc = ifelse(.data$Virulent > .data$Temperate, "virulent", "temperate"), Contig = .data$V1) %>%
      dplyr::select(c("Contig", "bacphlip_replicyc"))
  }
  if (version == "v1") {
    df <- df %>%
      dplyr::mutate(Contig = str_replace(.data$Contig, "-cat_[0-6]", ""))
  }

  return(df)
}


#' Read taxonomy table
#'
#' @param fin taxonomy annotation file
#' @param tool tool used to annotate taxonomy
#' @param version version of the tool used to annotate taxonomy
#'
#' @return data.frame
#' @export
#'
#' @examples
read_taxonomy <- function(fin, tool = "mmseqs", version = "v1") {
  . <- NULL
  if (tool == "mmseqs") {
    df <- data.table::fread(fin) %>%
      dplyr::filter(.data$taxid != 0) %>%
      data.table::setnames(c("genome_id", "superkingdom"), c("Contig", "Kingdom")) %>%
      data.table::setnames(colnames(.), str_to_title(colnames(.))) %>%
      dplyr::select(-"Taxid")
  } else if (tool == "mmseq_ictv") {
    df <- fread(fin, header = FALSE) %>%
      dplyr::select(c("V1", "V9")) %>%
      setnames(colnames(.), c("Contig", "taxonomy"))
  }
  if (version == "v1") {
    df <- df %>%
      dplyr::mutate(Contig = str_replace(.data$Contig, "-cat_[0-6]", ""))
  }
  return(df)
}


#' Read tblastx generated by Easyfig
#'
#' This is for adding links to gggenomes.
#'
#' @param fin 12.easyfig.out
#' @param seq1 query genome name
#' @param seq2 target genome name
#' @param max_evalue E-value threshold
#' @param min_bitscore Bitscore threshold
#'
#' @return data.frame
#' @export
#'
#' @examples
read_tblastx <- function(fin, seq1=0, seq2=0, max_evalue=0.001, min_bitscore=50) {
  df <- data.table::fread(fin,
                          header = FALSE,
                          col.names = c("seq_id", "seq_id2", "pident", "length", "mismatch", "gapopen", "start", "end", "start2", "end2", "evalue", "bitscore"))

  if (seq1 != 0) { df <- dplyr::mutate(seq_id = seq1) }
  if (seq2 != 0) { df <- dplyr::mutate(seq_id2 = seq2) }
  df <- df %>%
    dplyr::select(c("seq_id", "start", "end", "seq_id2", "start2", "end2", "evalue", "bitscore", "pident")) %>%
    dplyr::filter(.data$evalue < max_evalue) %>%
    dplyr::filter(.data$bitscore > min_bitscore)

  return(df)
}


#' Import vConTACT2 results
#'
#' @param fin File `genome_by_genome_overview.csv` created by vConTACT2
#' @param assembler_label "_NODE_" for SPAdes
#' @param version version of the tool used to annotate taxonomy
#'
#' @return list of results
#' @export
#'
#' @examples
read_vcontact2 <- function(fin, assembler_label = "_NODE_", version = "v1") {
  # TODO: need to tested using a test dataset
  . <- NULL
  df_vcontact <- fread(fin) %>%
    mutate(source=ifelse(str_detect(.data$Genome, assembler_label), "queryseq", "refseq")) %>%
    mutate(cluster_status=ifelse(str_detect(.data[["VC Status"]], "Overlap"), "Overlap", .data[["VC Status"]])) %>%
    mutate(cluster_status=factor(.data$cluster_status)) %>%
    setnames(colnames(.), paste0("vConTACT_", str_replace_all(colnames(.), " ", "_")))

  vc2_refs <- df_vcontact %>%
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label, negate = TRUE)) %>%
    dplyr::filter(.data$vConTACT_VC != "")

  # Get VC stats
  vcontact_stats <- df_vcontact %>%
    dplyr::group_by(.data$vConTACT_source, .data$vConTACT_cluster_status) %>%
    dplyr::summarise(seqs_with_VC = sum(.data$vConTACT_VC != ""),
                     seqs_without_VC = sum(.data$vConTACT_VC == "")) %>%
    # tidyr::gather(.data$seq_status, .data$num_seqs, -.data$vConTACT_source, -.data$vConTACT_cluster_status)
    tidyr::gather(seq_status, num_seqs, -vConTACT_source, -vConTACT_cluster_status)


  plot_vc_stats <- ggplot(vcontact_stats, aes(y=.data$num_seqs, x=.data$seq_status, fill=.data$vConTACT_source, label=.data$num_seqs)) +
    geom_point(aes(color=.data$vConTACT_source), alpha=0.5, size=3) +
    facet_grid(.data$vConTACT_cluster_status ~ .) +
    # theme(text = element_text(size = 12)) +
    geom_text_repel()

  # Annotate clusters using reference genomes
  vc2_contigs_vclst_anno <- df_vcontact %>%
    # Get distinct cluster ID of contigs in samples
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label)) %>%
    dplyr::filter(.data$vConTACT_VC != "") %>%
    dplyr::select(.data$vConTACT_VC) %>%
    distinct() %>%
    # Annotate cluster with reference taxonomy
    inner_join(vc2_refs, by = "vConTACT_VC")

  # Whether contig clusters were annotated by reference (1) or not (0)
  df_vcontact <- df_vcontact %>%
    mutate(vConTACT_classified=ifelse(.data$vConTACT_VC %in% vc2_contigs_vclst_anno$vConTACT_VC, 1, 0)) %>%
    # only choose assemblies
    dplyr::filter(str_detect(.data$vConTACT_Genome, assembler_label)) %>%
    mutate(vConTACT_VC2=ifelse(.data$vConTACT_VC_Status %in% c("Outlier", "Singleton"), .data$vConTACT_Genome, .data$vConTACT_VC)) %>%
    mutate(vConTACT_VC2=ifelse(str_detect(.data$vConTACT_VC_Status, "Overlap"), .data$vConTACT_VC_Status, .data$vConTACT_VC2)) %>%
    setnames("vConTACT_Genome", "virsorter2_contig_id") %>%
    mutate(vConTACT_VC2=str_replace_all(.data$vConTACT_VC2, "[\\s/]", "_")) %>%
    mutate(vConTACT_VC2=str_replace_all(.data$vConTACT_VC2, "[\\(\\)]", ""))

  if (version == "v1") {
    df_vcontact <- df_vcontact %>%
      dplyr::mutate(Contig = str_replace(.data$virsorter2_contig_id, "-cat_[0-6]", "")) %>%
      dplyr::select(-"virsorter2_contig_id")
  }

  vc2 <- list("vc_tbl" = df_vcontact,
              "vc_stats" = vcontact_stats,
              "vc_plot" = plot_vc_stats,
              "vc_annotated"= vc2_contigs_vclst_anno,
              "vc_refs" = vc2_refs)
  return(vc2)
}


#' Read vConTACT2 results and return as one data.frame
#'
#' @param fin File `genome_by_genome_overview.csv` created by vConTACT2
#' @param assembler_label "_NODE_" for SPAdes
#'
#' @return data.frame
#' @export
#'
#' @examples
read_vcontact2_simple <- function(fin, assembler_label="_NODE_") {
  . <- NULL
  df <- fread(fin) %>%
    dplyr::select(-"V1") %>%
    dplyr::filter(str_detect(.data$Genome, assembler_label)) %>%
    column_to_rownames("Genome") %>%
    setnames(colnames(.), paste0("vc_", str_replace_all(colnames(.), " ", "_"))) %>%
    rownames_to_column("contig_id")
  return(df)
}


#' Read VIBRANT results
#'
#' @param fin file of VIBRANT results
#'
#' @return data.table
#' @export
#'
#' @examples
read_vibrant <- function(fin) {
  . <- NULL
  df <- fread(fin) %>%
    setnames(colnames(.), c("Contig", "vibrant_replicyc", "vibrant_quality"))

  return(df)
}


#' Read VirSorter2 results
#'
#' Add "virsorter2_" prefix to column names.
#'
#' @param fin VirSorter2 results
#'
#' @return data.frame
#' @export
#'
#' @examples
read_virsorter2 <- function(fin) {
  . <- NULL
  dtbl <- fread(fin) %>%
    column_to_rownames("seqname") %>%
    setnames(colnames(.), paste0("virsorter2_", colnames(.))) %>%
    rownames_to_column("Contig")
  return(dtbl)
}

