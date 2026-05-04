#!/usr/bin/env Rscript
# ==============================================================================
# YAG RNA-seq Expression Comparison Analysis
# ==============================================================================
# Compares gene expression estimates between standard and enriched references
# Focuses on Y-linked ampliconic gene (YAG) families
# ==============================================================================

library(tximport)
library(tidyverse)
library(ggplot2)
library(scales)

# ------------------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- "/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
QUANT_DIR <- file.path(PROJECT_DIR, "salmon_quant")
OUTPUT_DIR <- file.path(PROJECT_DIR, "results")
YAG_MAPPING_FILE <- file.path(PROJECT_DIR, "data", "yag_gene_mapping.tsv")
TX2GENE_FILE <- file.path(PROJECT_DIR, "data", "tx2gene_mapping.tsv")

# Species information
SPECIES <- c("HomSap", "PanTro", "PanPan", "GorGor", "PonAbe", "PonPyg")
SPECIES_NAMES <- c(
  HomSap = "Human",
  PanTro = "Chimpanzee",
  PanPan = "Bonobo",
  GorGor = "Gorilla",
  PonAbe = "Sumatran orangutan",
  PonPyg = "Bornean orangutan"
)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

#' Create output directories
setup_directories <- function() {
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
  dir.create(file.path(OUTPUT_DIR, "tables"), showWarnings = FALSE)
}

#' Load YAG gene mapping from file
#' @param chrY_only If TRUE, filter to only chrY genes (true ampliconic copies)
#' @return Data frame with species, YAG family, chromosome, and gene ID
load_yag_mapping <- function(chrY_only = TRUE) {
  if (!file.exists(YAG_MAPPING_FILE)) {
    stop("YAG mapping file not found: ", YAG_MAPPING_FILE)
  }

  mapping <- read_tsv(YAG_MAPPING_FILE, show_col_types = FALSE)

  if (chrY_only) {
    mapping <- mapping %>% filter(Chromosome == "chrY")
  }

  return(mapping)
}

#' Load Salmon quantification results for a species
#' @param species Species code (e.g., "HomSap")
#' @param ref_type Reference type ("standard" or "enriched")
#' @return Data frame with transcript-level expression
load_salmon_quant <- function(species, ref_type) {
  quant_file <- file.path(QUANT_DIR, ref_type, species, "quant.sf")

  if (!file.exists(quant_file)) {
    warning(paste("Quantification file not found:", quant_file))
    return(NULL)
  }

  quant <- read_tsv(quant_file, show_col_types = FALSE) %>%
    mutate(
      species = species,
      species_name = SPECIES_NAMES[species],
      ref_type = ref_type
    )

  return(quant)
}

#' Load all quantification results
#' @return Combined data frame with all species and reference types
load_all_quants <- function() {
  results <- list()

  for (species in SPECIES) {
    for (ref_type in c("standard", "enriched")) {
      key <- paste(species, ref_type, sep = "_")
      quant <- load_salmon_quant(species, ref_type)
      if (!is.null(quant)) {
        results[[key]] <- quant
      }
    }
  }

  if (length(results) == 0) {
    stop("No quantification files found. Run Salmon quantification first.")
  }

  bind_rows(results)
}

#' Load transcript to gene mapping
#' @return Data frame with transcript_id, gene_symbol, species
load_tx2gene <- function() {
  if (!file.exists(TX2GENE_FILE)) {
    stop("tx2gene mapping file not found: ", TX2GENE_FILE,
         "\nRun: python scripts/create_tx2gene_mapping.py")
  }
  read_tsv(TX2GENE_FILE, show_col_types = FALSE)
}

#' Identify YAG transcripts using the mapping file
#' @param gene_names Vector of gene names (could be family names like VCY or specific IDs like VCY1B)
#' @param species Vector of species codes (same length as gene_names)
#' @param yag_mapping YAG mapping data frame from load_yag_mapping()
#' @return Data frame with gene family assignments
classify_yag <- function(gene_names, species, yag_mapping) {
  result <- tibble(gene = gene_names, species = species)

  # Create lookup table for family names (species + family -> family)
  family_lookup <- yag_mapping %>%
    select(Species, YAG_Family) %>%
    distinct() %>%
    mutate(matched_family = YAG_Family)

  # Create lookup table for specific gene IDs (species + gene_id -> family)
  gene_id_lookup <- yag_mapping %>%
    select(Species, Gene_ID, YAG_Family) %>%
    distinct()

  # First try: match gene name directly to YAG_Family (for extracted family names like VCY, TSPY)
  result <- result %>%
    left_join(
      family_lookup,
      by = c("species" = "Species", "gene" = "YAG_Family")
    )

  # Second try: for genes not matched, try matching to specific Gene_ID
  unmatched_idx <- is.na(result$matched_family)
  if (any(unmatched_idx)) {
    gene_id_matches <- result[unmatched_idx, ] %>%
      select(gene, species) %>%
      left_join(
        gene_id_lookup,
        by = c("species" = "Species", "gene" = "Gene_ID")
      )

    result$matched_family[unmatched_idx] <- gene_id_matches$YAG_Family
  }

  result <- result %>%
    mutate(yag_family = matched_family) %>%
    select(gene, species, yag_family)

  return(result)
}

#' Calculate expression summary per gene
#' @param quant_data Quantification data frame
#' @param tx2gene Transcript to gene mapping data frame
#' @return Summarized expression by gene
summarize_by_gene <- function(quant_data, tx2gene) {
  quant_data %>%
    left_join(tx2gene, by = c("Name" = "transcript_id", "species" = "species")) %>%
    filter(!is.na(gene_symbol)) %>%
    rename(gene = gene_symbol) %>%
    group_by(species, species_name, ref_type, gene) %>%
    summarize(
      TPM = sum(TPM),
      NumReads = sum(NumReads),
      n_transcripts = n(),
      .groups = "drop"
    )
}

# ------------------------------------------------------------------------------
# Analysis Functions
# ------------------------------------------------------------------------------

#' Classify YAGs and sum expression by family for one ref_type
#' @param gene_summary_subset gene_summary filtered to a single ref_type
#' @param yag_mapping YAG mapping data frame from load_yag_mapping()
#' @return Data frame summed by (species, species_name, yag_family)
classify_and_sum_by_family <- function(gene_summary_subset, yag_mapping) {
  yag_class <- classify_yag(gene_summary_subset$gene, gene_summary_subset$species, yag_mapping)
  gene_summary_subset$yag_family <- yag_class$yag_family

  gene_summary_subset %>%
    filter(!is.na(yag_family)) %>%
    group_by(species, species_name, yag_family) %>%
    summarize(
      TPM = sum(TPM, na.rm = TRUE),
      NumReads = sum(NumReads, na.rm = TRUE),
      .groups = "drop"
    )
}

#' Analyze YAG expression from enriched reference (for individual gene-level plots)
#' @param gene_summary Gene-level expression summary
#' @param yag_mapping YAG mapping data frame from load_yag_mapping()
#' @return Data frame with YAG expression from enriched reference
analyze_yag_expression <- function(gene_summary, yag_mapping) {
  enriched <- gene_summary %>%
    filter(ref_type == "enriched") %>%
    select(species, species_name, gene, TPM, NumReads, n_transcripts)

  yag_class <- classify_yag(enriched$gene, enriched$species, yag_mapping)
  enriched$yag_family <- yag_class$yag_family

  return(enriched)
}

#' Build standard vs enriched comparison table, summed by YAG family
#' @param gene_summary Gene-level expression summary
#' @param yag_mapping YAG mapping data frame from load_yag_mapping()
#' @return Data frame with TPM_std, TPM_enr, and fold change per species/family
build_yag_comparison <- function(gene_summary, yag_mapping) {
  std <- gene_summary %>%
    filter(ref_type == "standard") %>%
    classify_and_sum_by_family(yag_mapping) %>%
    rename(TPM_std = TPM, NumReads_std = NumReads)

  enr <- gene_summary %>%
    filter(ref_type == "enriched") %>%
    classify_and_sum_by_family(yag_mapping) %>%
    rename(TPM_enr = TPM, NumReads_enr = NumReads)

  full_join(std, enr, by = c("species", "species_name", "yag_family")) %>%
    mutate(
      TPM_std      = replace_na(TPM_std, 0),
      TPM_enr      = replace_na(TPM_enr, 0),
      NumReads_std = replace_na(NumReads_std, 0),
      NumReads_enr = replace_na(NumReads_enr, 0),
      TPM_log2fc   = log2((TPM_enr + 1) / (TPM_std + 1)),
      Reads_log2fc = log2((NumReads_enr + 1) / (NumReads_std + 1)),
      Reads_diff   = NumReads_enr - NumReads_std
    ) %>%
    arrange(species, yag_family)
}

#' Filter to YAG genes only
#' @param enriched_data Enriched reference YAG data frame
#' @return Filtered data frame with only classified YAG genes
filter_yag <- function(enriched_data) {
  enriched_data %>%
    filter(!is.na(yag_family))
}

# ------------------------------------------------------------------------------
# Visualization Functions
# ------------------------------------------------------------------------------

#' Plot YAG expression by family across species (bar plot)
#' @param yag_data YAG expression data
plot_yag_expression_bars <- function(yag_data) {
  p <- ggplot(yag_data, aes(x = yag_family, y = TPM, fill = yag_family)) +
    geom_col() +
    facet_wrap(~species_name, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    labs(
      title = "YAG Gene Family Expression Across Species",
      subtitle = "Using Enriched Reference with Custom YAG Transcripts",
      x = "YAG Gene Family",
      y = "TPM"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      strip.background = element_rect(fill = "white")
    )

  return(p)
}

#' Plot YAG expression by species (grouped bar)
#' @param yag_data YAG expression data
plot_yag_by_species <- function(yag_data) {
  p <- ggplot(yag_data, aes(x = species_name, y = TPM, fill = yag_family)) +
    geom_col(position = "dodge") +
    scale_y_continuous(labels = comma) +
    labs(
      title = "YAG Expression by Species",
      subtitle = "Enriched Reference",
      x = "Species",
      y = "TPM",
      fill = "YAG Family"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  return(p)
}

#' Plot heatmap of YAG expression across species
#' @param yag_data YAG expression data
plot_yag_heatmap <- function(yag_data) {
  p <- ggplot(yag_data, aes(x = species_name, y = yag_family, fill = log10(TPM + 1))) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(name = "Log10(TPM + 1)") +
    labs(
      title = "YAG Family Expression Across Species",
      subtitle = "Using Enriched Reference with Custom YAG Transcripts",
      x = "",
      y = "YAG Family"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  return(p)
}

#' Plot standard vs enriched read counts side by side per family/species
#' @param comparison Output of build_yag_comparison()
plot_std_vs_enriched <- function(comparison) {
  plot_data <- comparison %>%
    pivot_longer(c(NumReads_std, NumReads_enr), names_to = "ref_type", values_to = "NumReads") %>%
    mutate(ref_type = recode(ref_type, NumReads_std = "Standard", NumReads_enr = "Enriched"))

  p <- ggplot(plot_data, aes(x = yag_family, y = NumReads, fill = ref_type)) +
    geom_col(position = "dodge") +
    facet_wrap(~species_name, scales = "free_y") +
    scale_y_continuous(labels = comma) +
    scale_fill_manual(values = c(Standard = "#4472C4", Enriched = "#ED7D31")) +
    labs(
      title = "YAG Read Counts: Standard vs Enriched Reference",
      x = "YAG Family",
      y = "Read Count",
      fill = "Reference"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )

  return(p)
}

#' Plot log2 fold change of read counts (enriched / standard) per family/species
#' @param comparison Output of build_yag_comparison()
plot_foldchange_by_family <- function(comparison) {
  p <- ggplot(comparison, aes(x = yag_family, y = Reads_log2fc, fill = Reads_log2fc > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    facet_wrap(~species_name) +
    scale_fill_manual(values = c(`TRUE` = "#ED7D31", `FALSE` = "#4472C4"),
                      labels = c(`TRUE` = "Gain", `FALSE` = "Loss"),
                      name = NULL) +
    labs(
      title = "YAG Read Count Change: Enriched vs Standard Reference",
      subtitle = "Log2 fold change of read counts (enriched / standard), pseudocount = 1",
      x = "YAG Family",
      y = "Log2 Fold Change (Reads)"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      strip.background = element_rect(fill = "white")
    )

  return(p)
}

# ------------------------------------------------------------------------------
# Summary Statistics
# ------------------------------------------------------------------------------

#' Generate summary statistics table by family
#' @param yag_data YAG expression data
generate_summary_stats <- function(yag_data) {
  summary_stats <- yag_data %>%
    group_by(species_name, yag_family) %>%
    summarize(
      TPM = sum(TPM, na.rm = TRUE),
      NumReads = sum(NumReads, na.rm = TRUE),
      n_transcripts = sum(n_transcripts, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(species_name, yag_family)

  return(summary_stats)
}

#' Generate overall summary by species
#' @param yag_data YAG expression data
generate_species_summary <- function(yag_data) {
  species_summary <- yag_data %>%
    group_by(species_name) %>%
    summarize(
      n_yag_families = n_distinct(yag_family),
      total_TPM = sum(TPM, na.rm = TRUE),
      total_reads = sum(NumReads, na.rm = TRUE),
      top_family = yag_family[which.max(TPM)],
      top_family_TPM = max(TPM, na.rm = TRUE),
      .groups = "drop"
    )

  return(species_summary)
}

# ------------------------------------------------------------------------------
# Main Analysis Pipeline
# ------------------------------------------------------------------------------

main <- function() {
  cat("=== YAG RNA-seq Expression Analysis ===\n\n")

  # Setup
  setup_directories()

  # Load tx2gene mapping
  cat("Loading transcript to gene mapping...\n")
  tx2gene <- load_tx2gene()
  cat(sprintf("  Loaded %d transcript mappings\n", nrow(tx2gene)))

  # Load YAG gene mapping (chrY only for true ampliconic copies)
  cat("Loading YAG gene mapping...\n")
  yag_mapping <- load_yag_mapping(chrY_only = TRUE)
  cat(sprintf("  Loaded %d YAG gene mappings (chrY only)\n", nrow(yag_mapping)))

  # Load data
  cat("Loading Salmon quantification results...\n")
  all_quants <- load_all_quants()
  cat(sprintf("  Loaded %d transcript quantifications\n", nrow(all_quants)))

  # Summarize by gene
  cat("Summarizing expression by gene...\n")
  gene_summary <- summarize_by_gene(all_quants, tx2gene)
  cat(sprintf("  %d unique genes across all samples\n", n_distinct(gene_summary$gene)))

  # Analyze YAG expression from enriched reference (gene-level)
  cat("Analyzing YAG expression from enriched reference...\n")
  yag_expression <- analyze_yag_expression(gene_summary, yag_mapping)
  yag_data <- filter_yag(yag_expression)
  cat(sprintf("  Found %d YAG gene entries\n", nrow(yag_data)))

  if (nrow(yag_data) == 0) {
    warning("No YAG genes found. Check gene name patterns and transcript IDs.")
    return(invisible(NULL))
  }

  # Build standard vs enriched comparison (family-level)
  cat("Building standard vs enriched comparison...\n")
  comparison <- build_yag_comparison(gene_summary, yag_mapping)
  cat(sprintf("  %d family/species combinations\n", nrow(comparison)))

  # Generate summaries
  cat("\nGenerating summary statistics...\n")
  summary_stats <- generate_summary_stats(yag_data)
  species_summary <- generate_species_summary(yag_data)

  # Save tables
  write_csv(summary_stats, file.path(OUTPUT_DIR, "tables", "yag_summary_by_family.csv"))
  write_csv(species_summary, file.path(OUTPUT_DIR, "tables", "yag_summary_by_species.csv"))
  write_csv(yag_data, file.path(OUTPUT_DIR, "tables", "yag_expression_data.csv"))
  write_csv(comparison, file.path(OUTPUT_DIR, "tables", "yag_full_comparison.csv"))
  cat("  Saved summary tables to results/tables/\n")

  # Generate plots
  cat("\nGenerating visualizations...\n")

  p1 <- plot_yag_expression_bars(yag_data)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_by_family.pdf"), p1, width = 12, height = 10)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_by_family.png"), p1, width = 12, height = 10, dpi = 150)

  p2 <- plot_yag_by_species(yag_data)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_by_species.pdf"), p2, width = 12, height = 8)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_by_species.png"), p2, width = 12, height = 8, dpi = 150)

  p3 <- plot_yag_heatmap(yag_data)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_heatmap.pdf"), p3, width = 10, height = 6)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_expression_heatmap.png"), p3, width = 10, height = 6, dpi = 150)

  p4 <- plot_std_vs_enriched(comparison)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_std_vs_enriched.pdf"), p4, width = 14, height = 10)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_std_vs_enriched.png"), p4, width = 14, height = 10, dpi = 150)

  p5 <- plot_foldchange_by_family(comparison)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_foldchange_by_family.pdf"), p5, width = 12, height = 10)
  ggsave(file.path(OUTPUT_DIR, "figures", "yag_foldchange_by_family.png"), p5, width = 12, height = 10, dpi = 150)

  cat("  Saved figures to results/figures/\n")

  # Print summary
  cat("\n=== YAG Expression Summary by Species ===\n")
  print(species_summary)

  cat("\n=== Standard vs Enriched Comparison ===\n")
  print(comparison)

  cat("\n=== Analysis Complete ===\n")
  cat(sprintf("Results saved to: %s\n", OUTPUT_DIR))

  return(invisible(list(
    tx2gene = tx2gene,
    yag_mapping = yag_mapping,
    all_quants = all_quants,
    gene_summary = gene_summary,
    yag_data = yag_data,
    comparison = comparison,
    summary_stats = summary_stats,
    species_summary = species_summary
  )))
}

# Run if called directly
if (!interactive()) {
  results <- main()
}
