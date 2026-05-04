#!/usr/bin/env Rscript
# Plot YAG expression by species: standard reference bars, enriched TPM as black line

library(tidyverse)
library(ggplot2)
library(scales)

PROJECT_DIR <- "/storage/group/kdm16/default/kxp5629/proj/15_RNASeq"
OUTPUT_DIR  <- file.path(PROJECT_DIR, "results")

SPECIES_ORDER <- c(
  "Bonobo", "Chimpanzee", "Human", "Gorilla",
  "Bornean orangutan", "Sumatran orangutan"
)

comparison <- read_csv(
  file.path(PROJECT_DIR, "results", "tables", "yag_full_comparison.csv"),
  show_col_types = FALSE
) %>%
  mutate(species_name = factor(species_name, levels = SPECIES_ORDER))

p <- ggplot(comparison, aes(x = species_name, fill = yag_family)) +
  geom_col(aes(y = TPM_std), position = position_dodge(width = 0.9)) +
  geom_errorbar(
    aes(ymin = TPM_enr, ymax = TPM_enr, group = yag_family),
    color = "black",
    linewidth = 0.8,
    width = 0.7,
    position = position_dodge(width = 0.9)
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "YAG Expression by Species",
    subtitle = "Bars = standard reference; black line = enriched reference TPM",
    x = "Species",
    y = "TPM",
    fill = "YAG Family"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "figures", "yag_std_bars_enriched_line_by_species.pdf"), p, width = 12, height = 8)
ggsave(file.path(OUTPUT_DIR, "figures", "yag_std_bars_enriched_line_by_species.png"), p, width = 12, height = 8, dpi = 150)

cat("Saved to results/figures/yag_std_bars_enriched_line_by_species.{pdf,png}\n")
