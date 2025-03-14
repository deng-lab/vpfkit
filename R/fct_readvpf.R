#' Detect ViroProfiler pipeline version from contig IDs
#'
#' ViroProfiler v1 contigs have `-cat_N` suffixes from VirSorter2.
#'
#' @param contig_ids Character vector of contig identifiers
#' @return "v1" or "v2"
#' @noRd
.detect_vp_version <- function(contig_ids) {
  if (any(grepl("-cat_[0-6]$", contig_ids))) "v1" else "v2"
}


#' Read CATBAT results
#'
#' @param fin CATBAT taxonomy annotation file
#'
#' @return data.frame
#' @export
#'
read_catbat <- function(fin) {
  . <- NULL
  if (!file.exists(fin)) stop("CATBAT file not found: ", fin, call. = FALSE)
  df <- fread(fin, fill = TRUE, sep = "\t")
  if (!"# contig" %in% colnames(df)) stop("CATBAT file missing '# contig' column", call. = FALSE)
  df <- df %>%
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
read_checkv <- function(fin) {
  . <- NULL
  if (!file.exists(fin)) stop("CheckV file not found: ", fin, call. = FALSE)
  df <- fread(fin)
  required <- c("contig_id", "checkv_quality", "completeness")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("CheckV file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
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
read_coverm <- function(fpath, fbin2contig = 0, cov = 0) {
  . <- NULL
  if (!file.exists(fpath)) stop("CoverM file not found: ", fpath, call. = FALSE)
  df_abundance <- fread(fpath)
  if (!"Contig" %in% colnames(df_abundance)) stop("CoverM file missing 'Contig' column", call. = FALSE)
  df_abundance <- df_abundance %>%
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


#' Adjust abundance table by coverage fraction table
#'
#' @param df_abundance assay(tse, "counts") table
#' @param df_covfrac assay(tse, "covfrac") table
#' @param covfrac_threshold minimum covfrac (default: 0.5)
#'
#' @return matrix
#' @export
#'
abundance_adjust_by_covfrac <- function(df_abundance, df_covfrac, covfrac_threshold=0.5) {
  df_covfrac[df_covfrac<covfrac_threshold] <- 0
  df_covfrac[df_covfrac>=covfrac_threshold] <- 1
  df_abundance_adjusted <- df_abundance * df_covfrac
  abundance_matrix <- as.matrix(df_abundance_adjusted)
  return(abundance_matrix)
}


#' Convert reads_per_base to base_per_base (depth)
#'
#' @param fin path of reads_per_base file
#' @param reads_len length of clean reads
#'
#' @return matrix
#' @export
#'
rpb2bpb <- function(fin, reads_len=150) {
  df_rpb <- fread(fin) %>% column_to_rownames("Contig")
  df_bpb <- df_rpb * reads_len
  bpb_matrix <- as.matrix(df_bpb)
  return(bpb_matrix)
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
read_dvf <- function(fin, thr_score=0.9, thr_pvalue=0.01, thr_qvalue=0.01) {
  if (!file.exists(fin)) stop("DVF file not found: ", fin, call. = FALSE)
  df <- fread(fin)
  required <- c("name", "score", "pvalue", "qvalue")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("DVF file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
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
read_iphop <- function(fin) {
  if (!file.exists(fin)) stop("iPhop file not found: ", fin, call. = FALSE)
  df <- fread(fin)
  if (ncol(df) < 5) stop("iPhop file must have at least 5 columns", call. = FALSE)

  if (!("Contig" %in% colnames(df) && "iphop_genus" %in% colnames(df))) {
    df <- setnames(df, colnames(df)[1:5], c("Contig", "iphop_aai2ref", "iphop_genus", "iphop_score", "iphop_methods"))
  }
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
read_replicyc <- function(fin, tool = "bacphlip", version = "auto") {
  if (!file.exists(fin)) stop("Replication cycle file not found: ", fin, call. = FALSE)
  df <- fread(fin)

  if (tool == "bacphlip") {
    required <- c("V1", "Virulent", "Temperate")
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0) stop("Bacphlip file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    df <- df %>%
      dplyr::mutate(bacphlip_replicyc = ifelse(.data$Virulent > .data$Temperate, "virulent", "temperate"), Contig = .data$V1) %>%
      dplyr::select(c("Contig", "bacphlip_replicyc"))
  }
  if (version == "auto") version <- .detect_vp_version(df$Contig)
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
read_taxonomy <- function(fin, tool = "mmseqs", version = "auto") {
  . <- NULL
  if (!file.exists(fin)) stop("Taxonomy file not found: ", fin, call. = FALSE)
  if (tool == "mmseqs") {
    df <- data.table::fread(fin)
    required <- c("genome_id", "superkingdom", "taxid")
    missing <- setdiff(required, colnames(df))
    if (length(missing) > 0) stop("mmseqs taxonomy file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
    df <- df %>%
      dplyr::filter(.data$taxid != 0) %>%
      data.table::setnames(c("genome_id", "superkingdom"), c("Contig", "Kingdom")) %>%
      data.table::setnames(colnames(.), str_to_title(colnames(.))) %>%
      dplyr::select(-"Taxid")
  } else if (tool == "mmseq_ictv") {
    df <- fread(fin, header = FALSE) %>%
      dplyr::select(c("V1", "V9")) %>%
      setnames(colnames(.), c("Contig", "taxonomy"))
  }
  if (version == "auto") version <- .detect_vp_version(df$Contig)
  if (version == "v1") {
    df <- df %>%
      dplyr::mutate(Contig = str_replace(.data$Contig, "-cat_[0-6]", ""))
  }
  return(df)
}



#' Read mmseqs2 taxonomy table
#'
#' @param fin taxonomy annotation file
#'
#' @return data.frame
#' @export
#'
read_taxonomy2 <- function(fin) {
  . <- NULL
  if (!file.exists(fin)) stop("mmseqs2 taxonomy file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin)
  required <- c("contig_id", "taxa_id", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("mmseqs2 taxonomy file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
    dplyr::filter(.data$taxa_id != 0) %>%
    dplyr::select(c("contig_id", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
    data.table::setnames("contig_id", "Contig") %>%
    data.table::setnames(colnames(.), str_to_title(colnames(.)))
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
read_tblastx <- function(fin, seq1=0, seq2=0, max_evalue=0.001, min_bitscore=50) {
  df <- data.table::fread(fin,
                          header = FALSE,
                          col.names = c("seq_id", "seq_id2", "pident", "length", "mismatch", "gapopen", "start", "end", "start2", "end2", "evalue", "bitscore"))

  if (seq1 != 0) { df <- df %>% dplyr::mutate(seq_id = seq1) }
  if (seq2 != 0) { df <- df %>% dplyr::mutate(seq_id2 = seq2) }
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
read_vcontact2 <- function(fin, assembler_label = "_NODE_", version = "auto") {
  # TODO: need to tested using a test dataset
  . <- NULL
  if (!file.exists(fin)) stop("vConTACT2 file not found: ", fin, call. = FALSE)
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

  if (version == "auto") version <- .detect_vp_version(df_vcontact$virsorter2_contig_id)
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
read_vibrant <- function(fin) {
  . <- NULL
  if (!file.exists(fin)) stop("VIBRANT file not found: ", fin, call. = FALSE)
  df <- fread(fin)
  if (ncol(df) < 3) stop("VIBRANT file must have at least 3 columns (scaffold, type, quality)", call. = FALSE)
  df <- df %>%
    setnames(colnames(.), c("Contig", "vibrant_replicyc", "vibrant_quality")) %>%
    dplyr::mutate(Contig = str_replace_all(.data$Contig, "_fragment_[0-9]*", "")) %>%
    dplyr::mutate(vibrant_quality = factor(.data$vibrant_quality, levels = c("complete circular", "high quality draft", "low quality draft", "medium quality draft"))) %>%
    # sort by Quality
    dplyr::arrange("vibrant_quality") %>%
    # remove duplicated scaffold, keep the first one
    dplyr::group_by(.data$Contig) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

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
read_virsorter2 <- function(fin) {
  . <- NULL
  if (!file.exists(fin)) stop("VirSorter2 file not found: ", fin, call. = FALSE)
  dtbl <- fread(fin)
  required <- c("seqname", "max_score", "max_score_group")
  missing <- setdiff(required, colnames(dtbl))
  if (length(missing) > 0) stop("VirSorter2 file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  dtbl <- dtbl %>%
    column_to_rownames("seqname") %>%
    setnames(colnames(.), paste0("virsorter2_", colnames(.))) %>%
    rownames_to_column("Contig")
  return(dtbl)
}


#' Read DRAM-v annotations
#'
#' @param fin DRAM-v annotations.tsv file
#' @return data.frame with gene-level annotations
#' @export
read_dramv <- function(fin) {
  if (!file.exists(fin)) stop("DRAM-v file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin, na.strings = c("", "NA"))
  required <- c("fasta", "scaffold", "rank")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("DRAM-v file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
    dplyr::rename(gene_id = .data$fasta, Contig = .data$scaffold) %>%
    dplyr::mutate(
      dramv_category = dplyr::case_when(
        .data$rank == "A" ~ "AMG",
        .data$rank == "V" ~ "viral",
        .data$rank == "H" ~ "host",
        TRUE ~ "other"
      ),
      dramv_ko = .data$ko_id,
      dramv_pfam = .data$pfam_hits,
      dramv_vog = .data$vogdb,
      dramv_amg_flags = .data$amg_flags
    ) %>%
    dplyr::select(dplyr::any_of(c("Contig", "gene_id", "dramv_category", "dramv_ko",
                                   "dramv_pfam", "dramv_vog", "dramv_amg_flags",
                                   "auxiliary_score")))
  return(df)
}


#' Summarize DRAM-v annotations per contig
#'
#' @param df_dramv Output from read_dramv()
#' @return data.frame with per-contig summary columns
#' @noRd
.summarize_dramv <- function(df_dramv) {
  df_dramv %>%
    dplyr::group_by(.data$Contig) %>%
    dplyr::summarise(
      dramv_amg_count = sum(.data$dramv_category == "AMG", na.rm = TRUE),
      dramv_viral_gene_count = sum(.data$dramv_category == "viral", na.rm = TRUE),
      dramv_total_gene_count = dplyr::n(),
      .groups = "drop"
    )
}


#' Read pharokka functional annotations
#'
#' @param fin pharokka CDS output TSV file
#' @return data.frame with gene-level annotations
#' @export
read_pharokka <- function(fin) {
  if (!file.exists(fin)) stop("pharokka file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin, na.strings = c("", "NA"))
  required <- c("gene", "contig")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("pharokka file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
    dplyr::rename(gene_id = .data$gene, Contig = .data$contig) %>%
    dplyr::mutate(
      pharokka_function = .data$`function`,
      pharokka_phrog = .data$phrog,
      pharokka_category = .data$phrog_category,
      pharokka_card = .data$card_match,
      pharokka_vfdb = .data$vfdb_match
    ) %>%
    dplyr::select(dplyr::any_of(c("Contig", "gene_id", "pharokka_function", "pharokka_phrog",
                                   "pharokka_category", "pharokka_card", "pharokka_vfdb")))
  return(df)
}


#' Summarize pharokka annotations per contig
#'
#' @param df_pharokka Output from read_pharokka()
#' @return data.frame with per-contig summary columns
#' @noRd
.summarize_pharokka <- function(df_pharokka) {
  df_pharokka %>%
    dplyr::group_by(.data$Contig) %>%
    dplyr::summarise(
      pharokka_gene_count = dplyr::n(),
      pharokka_card_count = sum(!is.na(.data$pharokka_card), na.rm = TRUE),
      pharokka_vfdb_count = sum(!is.na(.data$pharokka_vfdb), na.rm = TRUE),
      .groups = "drop"
    )
}


#' Read geNomad virus summary
#'
#' @param fin geNomad virus summary TSV file
#' @return data.frame
#' @export
read_genomad <- function(fin) {
  if (!file.exists(fin)) stop("geNomad file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin)
  required <- c("seq_name", "virus_score")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("geNomad file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
    dplyr::rename(Contig = .data$seq_name) %>%
    dplyr::mutate(
      genomad_score = .data$virus_score,
      genomad_fdr = .data$fdr,
      genomad_topology = .data$topology,
      genomad_taxonomy = .data$taxonomy,
      genomad_n_hallmarks = .data$n_hallmarks
    ) %>%
    dplyr::select(dplyr::any_of(c("Contig", "genomad_score", "genomad_fdr",
                                   "genomad_topology", "genomad_taxonomy",
                                   "genomad_n_hallmarks")))
  return(df)
}


#' Read vRhyme viral bin assignments
#'
#' @param fin vRhyme bin-to-contig mapping file
#' @return data.frame
#' @export
read_vrhyme <- function(fin) {
  if (!file.exists(fin)) stop("vRhyme file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin)
  required <- c("bin", "contig")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("vRhyme file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  df <- df %>%
    dplyr::rename(Contig = .data$contig) %>%
    dplyr::mutate(vrhyme_bin = .data$bin) %>%
    dplyr::select(dplyr::any_of(c("Contig", "vrhyme_bin")))
  return(df)
}


#' Read PHIST host predictions
#'
#' @param fin PHIST predictions TSV file
#' @return data.frame
#' @export
read_phist <- function(fin) {
  if (!file.exists(fin)) stop("PHIST file not found: ", fin, call. = FALSE)
  df <- data.table::fread(fin)
  required <- c("Virus", "Host", "Score")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("PHIST file missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  has_taxonomy <- "Host_taxonomy" %in% colnames(df)
  df <- df %>%
    dplyr::rename(Contig = .data$Virus) %>%
    dplyr::mutate(
      phist_host = .data$Host,
      phist_score = .data$Score,
      phist_host_taxonomy = if (has_taxonomy) .data$Host_taxonomy else NA_character_
    ) %>%
    dplyr::select(dplyr::any_of(c("Contig", "phist_host", "phist_score", "phist_host_taxonomy")))
  return(df)
}


#' Merge a per-contig summary data.frame into a rownames-keyed feature annotation data.frame
#'
#' @param feature_anno data.frame with Contig rownames
#' @param summary_df data.frame with a Contig column to left-join
#' @return data.frame with Contig rownames
#' @noRd
.merge_summary_into_rowdata <- function(feature_anno, summary_df) {
  feature_anno %>%
    tibble::rownames_to_column("Contig") %>%
    dplyr::left_join(summary_df, by = "Contig") %>%
    tibble::column_to_rownames("Contig")
}


#' Create TSE object
#'
#' @param fin_abcount File `abundance_contigs_count.tsv.gz` created by coverm
#' @param fin_abtpm File `abundance_contigs_tpm.tsv.gz` created by coverm
#' @param fin_abtmm File `abundance_contigs_tmm.tsv.gz` created by coverm
#' @param fin_abcov File `abundance_contigs_covfrac.tsv.gz` created by coverm
#' @param fin_taxa File `taxa_mmseqs_formatted_all.tsv` created by coverm
#' @param fin_checkv File `quality_summary.tsv` created by checkv
#' @param fin_virsorter2 File `final-viral-score.tsv` created by VirSorter2
#' @param fin_vibrant File `VIBRANT_genome_quality_contigs.tsv` created by VIBRANT
#' @param fin_dvf File `dvf_virus.tsv.tsv` created by DVF
#' @param fin_replicyc File `putative_vcontigs_pref1.fasta.bacphlip` created by Replicyc
#' @param df_metadata Metadata table as data.frame
#' @param fin_genomad geNomad virus summary TSV (optional, NULL to skip)
#' @param fin_vrhyme vRhyme bin-to-contig mapping (optional, NULL to skip)
#' @param fin_phist PHIST host predictions TSV (optional, NULL to skip)
#' @param fin_dramv DRAM-v annotations TSV (optional, NULL to skip)
#' @param fin_pharokka pharokka CDS output TSV (optional, NULL to skip)
#'
#' @return TreeSummarizedExperiment object
#' @export
#' @importFrom S4Vectors SimpleList
#' @importFrom TreeSummarizedExperiment TreeSummarizedExperiment
#'
create_vpftse <- function(fin_abcount, fin_abtpm, fin_abtmm, fin_abcov, fin_taxa, fin_checkv, fin_virsorter2, fin_vibrant, fin_dvf, fin_replicyc, df_metadata=NULL, fin_genomad=NULL, fin_vrhyme=NULL, fin_phist=NULL, fin_dramv=NULL, fin_pharokka=NULL) {
  df_abcount <- read_coverm(fin_abcount) %>% as.matrix()
  df_abtpm <- read_coverm(fin_abtpm) %>% as.matrix()
  df_abcov <- read_coverm(fin_abcov) %>% as.matrix()
  df_abtmm <- read_coverm(fin_abtmm) %>% as.matrix()
  df_taxa <- read_taxonomy2(fin_taxa)
  df_checkv <- read_checkv(fin_checkv)
  df_virsorter2 <- read_virsorter2(fin_virsorter2)
  df_vibrant <- read_vibrant(fin_vibrant)
  df_dvf <- read_dvf(fin_dvf, thr_qvalue = 0.1)
  df_replicyc <- read_replicyc(fin_replicyc)

  feature_anno <- data.frame(Contig=rownames(df_abcount)) %>%
    dplyr::left_join(df_taxa, by = "Contig") %>%
    dplyr::left_join(df_checkv, by = "Contig") %>%
    dplyr::left_join(df_virsorter2, by = "Contig") %>%
    dplyr::left_join(df_vibrant, by = "Contig") %>%
    dplyr::left_join(df_dvf, by = "Contig") %>%
    dplyr::left_join(df_replicyc, by = "Contig")

  # Optional new tool annotations
  if (!is.null(fin_genomad)) {
    feature_anno <- feature_anno %>% dplyr::left_join(read_genomad(fin_genomad), by = "Contig")
  }
  if (!is.null(fin_vrhyme)) {
    feature_anno <- feature_anno %>% dplyr::left_join(read_vrhyme(fin_vrhyme), by = "Contig")
  }
  if (!is.null(fin_phist)) {
    feature_anno <- feature_anno %>% dplyr::left_join(read_phist(fin_phist), by = "Contig")
  }

  feature_anno <- feature_anno %>%
    tibble::column_to_rownames("Contig")

  # Optional gene-level annotations stored in metadata
  gene_annotations <- NULL
  if (!is.null(fin_dramv)) {
    gene_annotations <- read_dramv(fin_dramv)
    # Add per-contig summary to rowData
    dramv_summary <- .summarize_dramv(gene_annotations)
    feature_anno <- .merge_summary_into_rowdata(feature_anno, dramv_summary)
  }
  if (!is.null(fin_pharokka)) {
    df_pharokka <- read_pharokka(fin_pharokka)
    if (is.null(gene_annotations)) {
      gene_annotations <- df_pharokka
    } else {
      gene_annotations <- dplyr::bind_rows(gene_annotations, df_pharokka)
    }
    pharokka_summary <- .summarize_pharokka(df_pharokka)
    feature_anno <- .merge_summary_into_rowdata(feature_anno, pharokka_summary)
  }

  if (is.null(df_metadata)) {
    df_metadata <- data.frame(sample_id = colnames(df_abcount),
                              sample_name = colnames(df_abcount)) %>%
      tibble::column_to_rownames("sample_id") %>%
      MultiAssayExperiment::DataFrame()
  } else {
    df_metadata <- df_metadata %>%
      tibble::column_to_rownames("sample_id") %>%
      MultiAssayExperiment::DataFrame()
  }

  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = SimpleList(counts = df_abcount,
                  tpm = df_abtpm,
                  tmm = df_abtmm,
                  covfrac = df_abcov),
    colData = df_metadata,
    rowData = feature_anno)

  # Store gene-level annotations in metadata
  if (!is.null(gene_annotations)) {
    S4Vectors::metadata(tse)$gene_annotations <- gene_annotations
  }

  return(tse)
}


#' Create TSE object from ViroProfiler output directory
#'
#' Auto-discovers ViroProfiler output files from a standard directory structure
#' and assembles a TreeSummarizedExperiment without manual file path specification.
#'
#' @param vpdir Path to ViroProfiler output directory
#' @param df_metadata Optional sample metadata data.frame
#' @param version ViroProfiler version: "auto" (default), "v1", or "v2"
#' @return TreeSummarizedExperiment object
#' @export
batch_create_vpftse <- function(vpdir, df_metadata = NULL, version = "auto") {
  if (!dir.exists(vpdir)) stop("Directory not found: ", vpdir, call. = FALSE)

  # Helper to find a file by pattern, returns NULL if not found
  .find_file <- function(pattern, required = FALSE) {
    matches <- list.files(vpdir, pattern = pattern, full.names = TRUE, recursive = TRUE)
    if (length(matches) == 0) {
      if (required) stop("Required file not found matching pattern '", pattern, "' in ", vpdir, call. = FALSE)
      return(NULL)
    }
    matches[1]  # Use first match
  }

  # Required files
  fin_abcount <- .find_file("abundance_contigs_count", required = TRUE)
  fin_abtpm <- .find_file("abundance_contigs_tpm", required = TRUE)
  fin_abtmm <- .find_file("abundance_contigs_tmm", required = TRUE)
  fin_abcov <- .find_file("abundance_contigs_covfrac", required = TRUE)
  fin_taxa <- .find_file("taxa_mmseqs_formatted", required = TRUE)
  fin_checkv <- .find_file("quality_summary\\.tsv", required = TRUE)
  fin_virsorter2 <- .find_file("final-viral-score", required = TRUE)
  fin_vibrant <- .find_file("VIBRANT_genome_quality", required = TRUE)
  fin_dvf <- .find_file("dvf_virus", required = TRUE)
  fin_replicyc <- .find_file("\\.bacphlip$", required = TRUE)

  # Optional files
  fin_genomad <- .find_file("virus_summary\\.tsv")
  fin_vrhyme <- .find_file("bin_to_contig")
  fin_dramv <- .find_file("annotations\\.tsv")
  fin_pharokka <- .find_file("cds_final_merged_output")

  create_vpftse(
    fin_abcount = fin_abcount,
    fin_abtpm = fin_abtpm,
    fin_abtmm = fin_abtmm,
    fin_abcov = fin_abcov,
    fin_taxa = fin_taxa,
    fin_checkv = fin_checkv,
    fin_virsorter2 = fin_virsorter2,
    fin_vibrant = fin_vibrant,
    fin_dvf = fin_dvf,
    fin_replicyc = fin_replicyc,
    df_metadata = df_metadata,
    fin_genomad = fin_genomad,
    fin_vrhyme = fin_vrhyme,
    fin_dramv = fin_dramv,
    fin_pharokka = fin_pharokka
  )
}


#' Create TSE virome object
#'
#' @param tse TSE object
#' @return TreeSummarizedExperiment object
#' @export
#'
create_vpftse_vir <- function(tse) {
  rd <- rowData(tse)
  vir_mmseqs <- !is.na(rd$Domain)
  vir_checkv <- rd$checkv_quality %in% c("Complete", "High-quality", "Medium-quality")
  vir_virsorter2 <- rd$virsorter2_max_score_group %in% c("dsDNAphage", "NCLDV", "RNA", "ssDNA", "lavidaviridae")
  vir_dvf <- !is.na(rd$dvf_score)
  vir_vibrant <- !is.na(rd$vibrant_quality)

  # geNomad score (if available)
  vir_genomad <- if ("genomad_score" %in% colnames(rd)) {
    !is.na(rd$genomad_score) & rd$genomad_score >= 0.7
  } else {
    rep(FALSE, nrow(rd))
  }

  vir_all <- vir_mmseqs | vir_checkv | vir_virsorter2 | vir_vibrant | vir_dvf | vir_genomad

  tse_vir <- tse[vir_all,]

  return(tse_vir)
}
