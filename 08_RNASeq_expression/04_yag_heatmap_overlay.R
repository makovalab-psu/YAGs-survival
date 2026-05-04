#!/usr/bin/env Rscript
# YAG Expression Heatmap: Standard with Enriched Overlay

library(tidyverse)
library(ggplot2)

PROJECT_DIR <- "/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "figures")

# Load comparison data
comparison <- read_csv(file.path(PROJECT_DIR, "results", "tables", "yag_full_comparison.csv"),
                       show_col_types = FALSE)

# Create heatmap with standard dataset and enriched overlay
p <- ggplot(comparison, aes(x = species_name, y = yag_family)) +
  geom_tile(aes(fill = log10(TPM_std + 1)), color = "white") +
  geom_segment(aes(
    x = as.numeric(factor(species_name)) - 0.4,
    xend = as.numeric(factor(species_name)) + 0.4,
    y = as.numeric(factor(yag_family)) - 0.5 + (log10(TPM_enr + 1) / max(log10(comparison$TPM_enr + 1), na.rm = TRUE)),
    yend = as.numeric(factor(yag_family)) - 0.5 + (log10(TPM_enr + 1) / max(log10(comparison$TPM_enr + 1), na.rm = TRUE))
  ), linewidth = 0.8) +
  scale_fill_viridis_c(name = "Log10(TPM + 1)\n(Standard)") +
  labs(
    title = "YAG Family Expression: Standard Reference with Enriched Overlay",
    subtitle = "Colored tiles show standard reference TPM; black lines show enriched reference levels",
    x = "",
    y = "YAG Family"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

ggsave(file.path(OUTPUT_DIR, "yag_expression_heatmap_overlay.pdf"), p, width = 10, height = 6)
ggsave(file.path(OUTPUT_DIR, "yag_expression_heatmap_overlay.png"), p, width = 10, height = 6, dpi = 150)

cat("Saved heatmap overlay to results/figures/\n")
