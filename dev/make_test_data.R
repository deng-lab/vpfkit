# dev/make_test_data.R
# Generates a realistic test TSE for manual Shiny app testing.
# Output: data-raw/test_viroprofiler.rds
#
# Run: Rscript dev/make_test_data.R
#   or in R console: source("dev/make_test_data.R")

suppressPackageStartupMessages({
  library(TreeSummarizedExperiment)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(dplyr)
  library(tibble)
})

set.seed(2024)

# ── Dimensions ────────────────────────────────────────────────────────────────
N_CONTIGS <- 60
N_SAMPLES <- 8
contig_ids <- sprintf("contig_%03d", seq_len(N_CONTIGS))
sample_ids <- c(
  sprintf("Healthy_%02d", 1:4),
  sprintf("Disease_%02d", 1:4)
)

# ── Abundance matrices ────────────────────────────────────────────────────────
# Healthy samples have a different viral community than Disease samples.
# Contigs 1-20: enriched in Healthy; 21-40: enriched in Disease; 41-60: shared

base_counts <- matrix(0L, nrow = N_CONTIGS, ncol = N_SAMPLES,
                      dimnames = list(contig_ids, sample_ids))

for (i in seq_len(N_CONTIGS)) {
  for (j in seq_len(N_SAMPLES)) {
    is_healthy <- j <= 4
    if (i <= 20) {
      mu <- if (is_healthy) 800 else 50
    } else if (i <= 40) {
      mu <- if (is_healthy) 50 else 800
    } else {
      mu <- 300
    }
    # Length-scaled counts: longer contigs recruit more reads
    contig_len <- 10000 + (i - 1) * 3000
    mu_scaled <- mu * contig_len / 30000
    base_counts[i, j] <- rnbinom(1, mu = mu_scaled, size = 5)
  }
}

# TPM: reads per kilobase per million
contig_lengths <- 10000 + (seq_len(N_CONTIGS) - 1) * 3000
rpk <- sweep(base_counts, 1, contig_lengths / 1000, "/")
tpm_mat <- sweep(rpk, 2, colSums(rpk) / 1e6, "/")
tpm_mat <- round(tpm_mat, 4)

# TMM: simple column-sum normalisation (stands in for edgeR TMM)
tmm_mat <- sweep(base_counts, 2, colSums(base_counts) / 1e6, "/")
tmm_mat <- round(tmm_mat, 4)

# Coverage fraction: proportion of contig bases covered
covfrac_mat <- matrix(0, nrow = N_CONTIGS, ncol = N_SAMPLES,
                      dimnames = list(contig_ids, sample_ids))
for (i in seq_len(N_CONTIGS)) {
  for (j in seq_len(N_SAMPLES)) {
    if (base_counts[i, j] > 10) {
      covfrac_mat[i, j] <- round(runif(1, 0.6, 1.0), 3)
    } else if (base_counts[i, j] > 0) {
      covfrac_mat[i, j] <- round(runif(1, 0.05, 0.6), 3)
    }
  }
}

# ── CheckV annotations ────────────────────────────────────────────────────────
checkv_quality_levels <- c("Complete", "High-quality", "Medium-quality",
                            "Low-quality", "Not-determined")
checkv_quality <- sample(checkv_quality_levels,
                         N_CONTIGS, replace = TRUE,
                         prob = c(0.10, 0.20, 0.30, 0.25, 0.15))

completeness <- mapply(function(q) {
  switch(q,
    "Complete"        = round(runif(1, 99, 100), 1),
    "High-quality"    = round(runif(1, 90, 99), 1),
    "Medium-quality"  = round(runif(1, 50, 90), 1),
    "Low-quality"     = round(runif(1, 10, 50), 1),
    "Not-determined"  = 0
  )
}, checkv_quality)

row_checkv <- data.frame(
  Contig                   = contig_ids,
  checkv_contig_length     = contig_lengths,
  checkv_gene_count        = as.integer(contig_lengths / 1000 * rpois(N_CONTIGS, 1.2)),
  checkv_viral_genes       = as.integer(rpois(N_CONTIGS, 4)),
  checkv_host_genes        = as.integer(rpois(N_CONTIGS, 1)),
  checkv_quality           = factor(checkv_quality, levels = checkv_quality_levels),
  checkv_completeness      = completeness,
  checkv_contamination     = round(runif(N_CONTIGS, 0, 5), 1),
  checkv_provirus          = sample(c("Yes","No"), N_CONTIGS, replace=TRUE, prob=c(0.1,0.9)),
  stringsAsFactors = FALSE
)

# ── VirSorter2 ────────────────────────────────────────────────────────────────
vs2_groups <- c("dsDNAphage","dsDNAphage","dsDNAphage","NCLDV","RNA","ssDNA","lavidaviridae")
vs2_group  <- sample(vs2_groups, N_CONTIGS, replace = TRUE,
                     prob = c(0.50, 0.15, 0.10, 0.10, 0.05, 0.05, 0.05))

row_vs2 <- data.frame(
  Contig                    = contig_ids,
  virsorter2_max_score      = round(runif(N_CONTIGS, 0.5, 1.0), 3),
  virsorter2_max_score_group= vs2_group,
  virsorter2_min_score      = round(runif(N_CONTIGS, 0.3, 0.7), 3),
  stringsAsFactors = FALSE
)

# ── VIBRANT ───────────────────────────────────────────────────────────────────
vibrant_quality_levels <- c("complete circular","high quality draft",
                             "medium quality draft","low quality draft")
row_vibrant <- data.frame(
  Contig           = sample(contig_ids, N_CONTIGS * 0.7),
  vibrant_replicyc = sample(c("lytic","lysogenic"), N_CONTIGS * 0.7, replace=TRUE,
                            prob = c(0.7, 0.3)),
  vibrant_quality  = factor(
    sample(vibrant_quality_levels, N_CONTIGS * 0.7, replace=TRUE,
           prob = c(0.15, 0.35, 0.35, 0.15)),
    levels = vibrant_quality_levels),
  stringsAsFactors = FALSE
)

# ── DVF ───────────────────────────────────────────────────────────────────────
dvf_idx <- sample(N_CONTIGS, N_CONTIGS * 0.6)
row_dvf <- data.frame(
  Contig    = contig_ids[dvf_idx],
  dvf_score = round(runif(length(dvf_idx), 0.9, 1.0), 4),
  stringsAsFactors = FALSE
)

# ── Bacphlip (replication cycle) ──────────────────────────────────────────────
row_bacphlip <- data.frame(
  Contig            = contig_ids,
  bacphlip_replicyc = sample(c("virulent","temperate"), N_CONTIGS, replace=TRUE,
                             prob = c(0.70, 0.30)),
  stringsAsFactors = FALSE
)

# ── iPhop (host prediction) ───────────────────────────────────────────────────
host_genera <- c(
  "Escherichia","Salmonella","Klebsiella","Pseudomonas","Staphylococcus",
  "Streptococcus","Lactobacillus","Bifidobacterium","Bacteroides","Clostridium"
)
iphop_idx <- sample(N_CONTIGS, N_CONTIGS * 0.65)
row_iphop <- data.frame(
  Contig          = contig_ids[iphop_idx],
  iphop_aai2ref   = sample(c("Ref_A","Ref_B","Ref_C","Ref_D"), length(iphop_idx), replace=TRUE),
  iphop_genus     = paste0("g__", sample(host_genera, length(iphop_idx), replace=TRUE)),
  iphop_score     = round(runif(length(iphop_idx), 75, 99), 1),
  iphop_methods   = sample(c("CHERRY;CRISPR","CHERRY","CRISPR","BLASTN"),
                           length(iphop_idx), replace=TRUE),
  stringsAsFactors = FALSE
)

# ── Taxonomy (mmseqs2) ────────────────────────────────────────────────────────
viral_families <- c(
  "Myoviridae","Siphoviridae","Podoviridae","Microviridae",
  "Inoviridae","Tectiviridae","Leviviridae","Cystoviridae"
)
viral_genera <- c(
  "Tequatrovirus","Kayfunavirus","Bruynoghevirus","Jerseyquatrovirus",
  "Mosigvirus","Peduovirus","Dhakavirus","Slopekaelivirus"
)
viral_orders <- c("Caudovirales","Petitvirales","Tubulavirales","Durnavirales")

taxa_idx <- sample(N_CONTIGS, N_CONTIGS * 0.80)
family_draw <- sample(viral_families, length(taxa_idx), replace = TRUE)
genus_draw  <- sample(viral_genera,  length(taxa_idx), replace = TRUE)
order_draw  <- sample(viral_orders,  length(taxa_idx), replace = TRUE)

row_taxa <- data.frame(
  Contig     = contig_ids[taxa_idx],
  Domain     = "Viruses",
  Kingdom    = "Duplodnaviria",
  Phylum     = sample(c("Uroviricota","Phixviricota","Loebviricota"), length(taxa_idx), replace=TRUE),
  Class      = sample(c("Caudoviricetes","Malgrandaviricetes","Faserviricetes"), length(taxa_idx), replace=TRUE),
  Order      = order_draw,
  Family     = family_draw,
  Genus      = genus_draw,
  Species    = paste0(genus_draw, " phage ", sample(1:999, length(taxa_idx), replace=TRUE)),
  stringsAsFactors = FALSE
)

# ── geNomad ───────────────────────────────────────────────────────────────────
row_genomad <- data.frame(
  Contig              = contig_ids,
  genomad_score       = round(runif(N_CONTIGS, 0.60, 1.00), 4),
  genomad_fdr         = round(runif(N_CONTIGS, 0, 0.05), 4),
  genomad_topology    = sample(c("linear","circular","Provirus"), N_CONTIGS, replace=TRUE,
                               prob = c(0.65, 0.25, 0.10)),
  genomad_taxonomy    = paste0("Viruses;", sample(viral_families, N_CONTIGS, replace=TRUE)),
  genomad_n_hallmarks = as.integer(rpois(N_CONTIGS, 3)),
  stringsAsFactors = FALSE
)

# ── vRhyme (viral bins) ───────────────────────────────────────────────────────
n_bins <- 15
bin_assignment <- sample(
  c(sprintf("vRhyme_unbinned_%03d", seq_len(N_CONTIGS - n_bins * 3)),
    rep(sprintf("vRhyme_bin_%02d", seq_len(n_bins)), each = 3)),
  N_CONTIGS
)
row_vrhyme <- data.frame(
  Contig     = contig_ids,
  vrhyme_bin = bin_assignment,
  stringsAsFactors = FALSE
)

# ── PHIST (host prediction) ───────────────────────────────────────────────────
phist_idx <- sample(N_CONTIGS, N_CONTIGS * 0.55)
host_species <- c(
  "Escherichia coli","Salmonella enterica","Klebsiella pneumoniae",
  "Pseudomonas aeruginosa","Staphylococcus aureus","Streptococcus pneumoniae"
)
row_phist <- data.frame(
  Contig               = contig_ids[phist_idx],
  phist_host           = sample(host_species, length(phist_idx), replace=TRUE),
  phist_score          = as.integer(sample(5:100, length(phist_idx), replace=TRUE)),
  phist_host_taxonomy  = paste0("Bacteria;Proteobacteria;Gammaproteobacteria;",
                                sample(c("Enterobacterales","Pseudomonadales"),
                                       length(phist_idx), replace=TRUE)),
  stringsAsFactors = FALSE
)

# ── DRAM-v gene annotations ───────────────────────────────────────────────────
amg_kos   <- c("K00001","K00003","K00005","K00008","K00012","K00016",
                "K00018","K00022","K00024","K00030","K00150","K00175")
pfam_ids  <- c("PF00171","PF00227","PF00574","PF01000","PF01176","PF01554",
               "PF02310","PF02861","PF03600","PF04851","PF05175","PF06268")
vog_ids   <- c("VOG00001","VOG00003","VOG00005","VOG00007","VOG00009","VOG00011")

dramv_genes <- do.call(rbind, lapply(sample(contig_ids, 45), function(ctg) {
  n_genes <- sample(3:8, 1)
  data.frame(
    Contig         = ctg,
    gene_id        = paste0(ctg, "_gene_", seq_len(n_genes)),
    dramv_category = sample(c("AMG","viral","viral","host","other"),
                            n_genes, replace=TRUE, prob=c(0.15,0.45,0.20,0.10,0.10)),
    dramv_ko       = sample(c(amg_kos, NA, NA, NA), n_genes, replace=TRUE),
    dramv_pfam     = sample(c(pfam_ids, NA, NA), n_genes, replace=TRUE),
    dramv_vog      = sample(c(vog_ids, NA, NA, NA), n_genes, replace=TRUE),
    dramv_amg_flags= sample(c("M","MB","MBF", NA, NA, NA), n_genes, replace=TRUE),
    stringsAsFactors = FALSE
  )
}))

# ── pharokka gene annotations ─────────────────────────────────────────────────
phrog_functions <- c(
  "tail fiber protein","major capsid protein","DNA polymerase",
  "terminase large subunit","portal protein","baseplate protein",
  "lysis","integration/excision","replication","other","unknown function"
)
card_hits <- c(NA, NA, NA, NA, NA,
               "ARO:3000026","ARO:3000114","ARO:3000389","ARO:3001122")
vfdb_hits <- c(NA, NA, NA, NA, NA, NA,
               "VFG000001","VFG000023","VFG000048","VFG000102")

pharokka_genes <- do.call(rbind, lapply(sample(contig_ids, 50), function(ctg) {
  n_genes <- sample(4:10, 1)
  data.frame(
    Contig             = ctg,
    gene_id            = paste0(ctg, "_phg_", seq_len(n_genes)),
    pharokka_function  = sample(phrog_functions, n_genes, replace=TRUE),
    pharokka_phrog     = paste0("phrog_", sample(1:10000, n_genes, replace=TRUE)),
    pharokka_category  = sample(c("moron, auxiliary metabolic gene and host takeover",
                                  "connector","DNA, RNA and nucleotide metabolism",
                                  "head and packaging","integration and excision",
                                  "lysis","tail","transcription regulation","other","unknown"),
                                n_genes, replace=TRUE),
    pharokka_card      = sample(card_hits, n_genes, replace=TRUE),
    pharokka_vfdb      = sample(vfdb_hits, n_genes, replace=TRUE),
    stringsAsFactors = FALSE
  )
}))

gene_annotations <- dplyr::bind_rows(dramv_genes, pharokka_genes)

# ── DRAM-v per-contig summary ─────────────────────────────────────────────────
dramv_summary <- dramv_genes %>%
  dplyr::group_by(Contig) %>%
  dplyr::summarise(
    dramv_amg_count        = sum(dramv_category == "AMG", na.rm=TRUE),
    dramv_viral_gene_count = sum(dramv_category == "viral", na.rm=TRUE),
    dramv_total_gene_count = dplyr::n(),
    .groups = "drop"
  )

# pharokka per-contig summary
pharokka_summary <- pharokka_genes %>%
  dplyr::group_by(Contig) %>%
  dplyr::summarise(
    pharokka_gene_count  = dplyr::n(),
    pharokka_card_count  = sum(!is.na(pharokka_card), na.rm=TRUE),
    pharokka_vfdb_count  = sum(!is.na(pharokka_vfdb), na.rm=TRUE),
    .groups = "drop"
  )

# ── Assemble rowData ──────────────────────────────────────────────────────────
feature_anno <- data.frame(Contig = contig_ids, stringsAsFactors = FALSE) %>%
  dplyr::left_join(row_taxa,     by = "Contig") %>%
  dplyr::left_join(row_checkv,   by = "Contig") %>%
  dplyr::left_join(row_vs2,      by = "Contig") %>%
  dplyr::left_join(row_vibrant,  by = "Contig") %>%
  dplyr::left_join(row_dvf,      by = "Contig") %>%
  dplyr::left_join(row_bacphlip, by = "Contig") %>%
  dplyr::left_join(row_iphop,    by = "Contig") %>%
  dplyr::left_join(row_genomad,  by = "Contig") %>%
  dplyr::left_join(row_vrhyme,   by = "Contig") %>%
  dplyr::left_join(row_phist,    by = "Contig") %>%
  dplyr::left_join(dramv_summary,   by = "Contig") %>%
  dplyr::left_join(pharokka_summary, by = "Contig") %>%
  tibble::column_to_rownames("Contig")

# ── colData ───────────────────────────────────────────────────────────────────
col_meta <- data.frame(
  sample_name  = sample_ids,
  group        = c(rep("Healthy", 4), rep("Disease", 4)),
  total_reads  = as.integer(colSums(base_counts)),
  stringsAsFactors = FALSE,
  row.names = sample_ids
)

# ── Assemble TSE ──────────────────────────────────────────────────────────────
tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = S4Vectors::SimpleList(
    counts  = base_counts,
    tpm     = tpm_mat,
    tmm     = tmm_mat,
    covfrac = covfrac_mat
  ),
  colData = S4Vectors::DataFrame(col_meta),
  rowData = feature_anno
)

S4Vectors::metadata(tse)$gene_annotations <- gene_annotations

# ── Save ──────────────────────────────────────────────────────────────────────
out_path <- file.path("data-raw", "test_viroprofiler.rds")
saveRDS(tse, out_path)

cat(sprintf(
  "Saved: %s\n  Contigs : %d\n  Samples : %d  (%s)\n  Assays  : %s\n  rowData : %d columns\n  Genes   : %d annotations (%d DRAM-v, %d pharokka)\n",
  out_path,
  nrow(tse),
  ncol(tse),
  paste(unique(col_meta$group), collapse=" / "),
  paste(SummarizedExperiment::assayNames(tse), collapse=", "),
  ncol(SummarizedExperiment::rowData(tse)),
  nrow(gene_annotations),
  nrow(dramv_genes),
  nrow(pharokka_genes)
))

# ==============================================================================
# Dataset 2: Ocean/Soil environmental viromes (for Compare datasets tab)
# Key differences from Dataset 1:
#   - 40 contigs, 6 samples (Ocean_01-03, Soil_01-03)
#   - Higher read counts (marine viromes are richer)
#   - Lower diversity (dominated by a few families)
#   - Partially overlapping taxonomy: shares Myoviridae/Siphoviridae with D1,
#     adds Phycodnaviridae/Mimiviridae/Marseilleviridae (large dsDNA viruses)
# ==============================================================================

set.seed(2025)

N2_CONTIGS <- 40
N2_SAMPLES <- 6
contig2_ids <- sprintf("vOTU_%03d", seq_len(N2_CONTIGS))
sample2_ids <- c(sprintf("Ocean_%02d", 1:3), sprintf("Soil_%02d", 1:3))

# ── Abundance (higher counts, lower evenness) ──────────────────────────────────
# Ocean samples dominated by contigs 1-15; Soil by contigs 16-30; 31-40 shared
base_counts2 <- matrix(0L, nrow = N2_CONTIGS, ncol = N2_SAMPLES,
                       dimnames = list(contig2_ids, sample2_ids))
for (i in seq_len(N2_CONTIGS)) {
  clen2 <- 8000 + (i - 1) * 2500
  for (j in seq_len(N2_SAMPLES)) {
    is_ocean <- j <= 3
    if (i <= 15) {
      mu <- if (is_ocean) 2000 else 80
    } else if (i <= 30) {
      mu <- if (is_ocean) 80 else 2000
    } else {
      mu <- 500
    }
    base_counts2[i, j] <- rnbinom(1, mu = mu * clen2 / 25000, size = 3)  # size=3 → more overdispersed
  }
}

contig2_lengths <- 8000 + (seq_len(N2_CONTIGS) - 1) * 2500
rpk2   <- sweep(base_counts2, 1, contig2_lengths / 1000, "/")
tpm2   <- round(sweep(rpk2, 2, colSums(rpk2) / 1e6, "/"), 4)
tmm2   <- round(sweep(base_counts2, 2, colSums(base_counts2) / 1e6, "/"), 4)
covfrac2 <- matrix(0, nrow = N2_CONTIGS, ncol = N2_SAMPLES,
                   dimnames = list(contig2_ids, sample2_ids))
for (i in seq_len(N2_CONTIGS)) {
  for (j in seq_len(N2_SAMPLES)) {
    if (base_counts2[i, j] > 10) {
      covfrac2[i, j] <- round(runif(1, 0.6, 1.0), 3)
    } else if (base_counts2[i, j] > 0) {
      covfrac2[i, j] <- round(runif(1, 0.05, 0.6), 3)
    }
  }
}

# ── CheckV ────────────────────────────────────────────────────────────────────
# Lower quality overall (environmental metagenomes often more fragmented)
checkv2_levels <- c("Complete","High-quality","Medium-quality","Low-quality","Not-determined")
checkv2_quality <- sample(checkv2_levels, N2_CONTIGS, replace=TRUE,
                           prob=c(0.05, 0.10, 0.25, 0.40, 0.20))
checkv2_completeness <- mapply(function(q) {
  switch(q,
    "Complete"        = round(runif(1, 99, 100), 1),
    "High-quality"    = round(runif(1, 90, 99), 1),
    "Medium-quality"  = round(runif(1, 50, 90), 1),
    "Low-quality"     = round(runif(1, 10, 50), 1),
    "Not-determined"  = 0
  )
}, checkv2_quality)

row_checkv2 <- data.frame(
  Contig               = contig2_ids,
  checkv_contig_length = contig2_lengths,
  checkv_gene_count    = as.integer(contig2_lengths / 1000 * rpois(N2_CONTIGS, 1.0)),
  checkv_viral_genes   = as.integer(rpois(N2_CONTIGS, 3)),
  checkv_host_genes    = as.integer(rpois(N2_CONTIGS, 1)),
  checkv_quality       = factor(checkv2_quality, levels = checkv2_levels),
  checkv_completeness  = checkv2_completeness,
  checkv_contamination = round(runif(N2_CONTIGS, 0, 8), 1),
  checkv_provirus      = sample(c("Yes","No"), N2_CONTIGS, replace=TRUE, prob=c(0.15,0.85)),
  stringsAsFactors = FALSE
)

# ── VirSorter2 ────────────────────────────────────────────────────────────────
row_vs2_2 <- data.frame(
  Contig                     = contig2_ids,
  virsorter2_max_score       = round(runif(N2_CONTIGS, 0.4, 1.0), 3),
  virsorter2_max_score_group = sample(c("dsDNAphage","NCLDV","RNA","ssDNA"), N2_CONTIGS,
                                      replace=TRUE, prob=c(0.40,0.30,0.15,0.15)),
  virsorter2_min_score       = round(runif(N2_CONTIGS, 0.2, 0.6), 3),
  stringsAsFactors = FALSE
)

# ── Taxonomy: partially overlapping with Dataset 1 ─────────────────────────────
# Shared families (also in D1): Myoviridae, Siphoviridae
# Unique to D2 (large dsDNA environmental viruses): Phycodnaviridae, Mimiviridae,
#   Marseilleviridae, Iridoviridae
d2_lookup <- data.frame(
  Family = c("Myoviridae",    "Siphoviridae",  "Phycodnaviridae", "Mimiviridae",   "Marseilleviridae", "Iridoviridae"),
  Genus  = c("Tequatrovirus", "Kayfunavirus",  "Prasinovirus",    "Megavirus",     "Marseillevirus",   "Chunkyvirus"),
  Order  = c("Caudovirales",  "Caudovirales",  "Algavirales",     "Imitervirales", "Marseillevirales", "Pimascovirales"),
  prob   = c(0.32,             0.18,            0.20,              0.20,            0.05,               0.05),
  stringsAsFactors = FALSE
)

taxa2_idx <- sample(N2_CONTIGS, N2_CONTIGS * 0.85)
n2_taxa   <- length(taxa2_idx)
fam2_draw   <- sample(d2_lookup$Family, n2_taxa, replace=TRUE, prob=d2_lookup$prob)
genus2_draw <- d2_lookup$Genus[match(fam2_draw, d2_lookup$Family)]
order2_draw <- d2_lookup$Order[match(fam2_draw, d2_lookup$Family)]

row_taxa2 <- data.frame(
  Contig  = contig2_ids[taxa2_idx],
  Domain  = "Viruses",
  Kingdom = sample(c("Duplodnaviria","Varidnaviria"), n2_taxa, replace=TRUE, prob=c(0.45,0.55)),
  Phylum  = sample(c("Uroviricota","Nucleocytoviricota","Pisuviricota"), n2_taxa,
                   replace=TRUE, prob=c(0.45,0.45,0.10)),
  Class   = sample(c("Caudoviricetes","Megaviricetes","Pisoniviricetes"), n2_taxa,
                   replace=TRUE, prob=c(0.45,0.45,0.10)),
  Order   = order2_draw,
  Family  = fam2_draw,
  Genus   = genus2_draw,
  Species = paste0(genus2_draw, " virus ", sample(1:999, n2_taxa, replace=TRUE)),
  stringsAsFactors = FALSE
)

# ── Assemble rowData ──────────────────────────────────────────────────────────
feature_anno2 <- data.frame(Contig = contig2_ids, stringsAsFactors = FALSE) %>%
  dplyr::left_join(row_taxa2,   by = "Contig") %>%
  dplyr::left_join(row_checkv2, by = "Contig") %>%
  dplyr::left_join(row_vs2_2,   by = "Contig") %>%
  tibble::column_to_rownames("Contig")

# ── colData ───────────────────────────────────────────────────────────────────
col_meta2 <- data.frame(
  sample_name = sample2_ids,
  group       = c(rep("Ocean", 3), rep("Soil", 3)),
  total_reads = as.integer(colSums(base_counts2)),
  stringsAsFactors = FALSE,
  row.names = sample2_ids
)

# ── Assemble TSE ──────────────────────────────────────────────────────────────
tse2 <- TreeSummarizedExperiment::TreeSummarizedExperiment(
  assays = S4Vectors::SimpleList(
    counts  = base_counts2,
    tpm     = tpm2,
    tmm     = tmm2,
    covfrac = covfrac2
  ),
  colData = S4Vectors::DataFrame(col_meta2),
  rowData = feature_anno2
)

out_path2 <- file.path("data-raw", "test_viroprofiler_2.rds")
saveRDS(tse2, out_path2)

cat(sprintf(
  "Saved: %s\n  Contigs : %d\n  Samples : %d  (%s)\n  Assays  : %s\n  rowData : %d columns\n  Families: %s\n",
  out_path2,
  nrow(tse2),
  ncol(tse2),
  paste(unique(col_meta2$group), collapse=" / "),
  paste(SummarizedExperiment::assayNames(tse2), collapse=", "),
  ncol(SummarizedExperiment::rowData(tse2)),
  paste(sort(unique(stats::na.omit(row_taxa2$Family))), collapse=", ")
))
