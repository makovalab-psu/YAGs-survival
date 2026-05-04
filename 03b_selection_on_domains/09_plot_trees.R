#!/usr/bin/env Rscript
# Step 09 — Plot one phylogenetic tree per domain.
# Roots on the appropriate outgroup clade, then midpoints the two root branches
# (same approach as plot_combined_trees_20260329.r).
#
# Usage:
#   Rscript 09_plot_trees.R                    # plot all output/*/tree.treefile
#   Rscript 09_plot_trees.R CDY_Chromodomain   # plot one specific subdirectory

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtree)
  library(ape)
  library(dplyr)
  library(glue)
})

script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) getwd()
)
setwd(script_dir)

# ---------------------------------------------------------------------------
# Species colour scheme (from plot_combined_trees_20260329.r)
# ---------------------------------------------------------------------------
species_colors <- c(
  "Human"                = "#969696",
  "Bonobo"               = "#fb2d32",
  "Chimpanzee"           = "#c3000e",
  "Gorilla"              = "#006e81",
  "Orangutan (Sumatran)" = "#fd8a1a",
  "Orangutan (Bornean)"  = "#f9cf44",
  "Siamang"              = "#0e205e"
)

extract_species <- function(label) {
  map <- c(
    "HomSap" = "Human",      "PanPan" = "Bonobo",
    "PanTro" = "Chimpanzee", "GorGor" = "Gorilla",
    "PonAbe" = "Orangutan (Sumatran)",
    "PonPyg" = "Orangutan (Bornean)",
    "SymSyn" = "Siamang"
  )
  prefix <- substr(gsub("^'|'$", "", label), 1, 6)
  ifelse(prefix %in% names(map), map[prefix], "Unknown")
}

# ---------------------------------------------------------------------------
# Outgroup pattern per directory name
# ---------------------------------------------------------------------------
outgroup_pattern <- function(dir_name) {
  if (grepl("^CDY",      dir_name)) return("CDYL")       # matches CDYL but CDYL2 absent
  if (grepl("^RBMY",     dir_name)) return("RBMX")
  if (grepl("^DAZ",      dir_name)) return("DAZL")
  if (grepl("^HSFY|^TSPY", dir_name)) return("SymSyn")
  return(NULL)
}

# ---------------------------------------------------------------------------
# Root on outgroup clade, then midpoint the two root branches
# ---------------------------------------------------------------------------
root_and_midpoint <- function(tree, pattern) {
  og_tips <- grep(pattern, tree$tip.label, value = TRUE)
  if (length(og_tips) == 0) {
    warning(glue("No outgroup tips matching '{pattern}' — tree left unrooted"))
    return(ladderize(tree, right = FALSE))
  }
  tree <- root(tree, outgroup = og_tips, resolve.root = TRUE, edgelabel = TRUE)
  root_node  <- length(tree$tip.label) + 1
  root_edges <- which(tree$edge[, 1] == root_node)
  if (length(root_edges) == 2) {
    total <- sum(tree$edge.length[root_edges])
    tree$edge.length[root_edges] <- total / 2
  }
  ladderize(tree, right = FALSE)
}

# ---------------------------------------------------------------------------
# Build and save one tree plot
# ---------------------------------------------------------------------------
plot_one_tree <- function(treefile, dir_name) {
  tree <- tryCatch(read.tree(treefile), error = function(e) {
    message(glue("  Cannot read {treefile}: {e$message}"))
    return(NULL)
  })
  if (is.null(tree)) return(invisible(NULL))

  og_pat <- outgroup_pattern(dir_name)
  if (!is.null(og_pat)) {
    tree <- root_and_midpoint(tree, og_pat)
  } else {
    tree <- ladderize(tree, right = FALSE)
  }

  n_tips     <- length(tree$tip.label)
  show_labels <- n_tips <= 50

  # Tip metadata
  tip_df <- data.frame(
    label   = tree$tip.label,
    species = sapply(tree$tip.label, extract_species),
    stringsAsFactors = FALSE
  )

  # Base tree
  p <- ggtree(tree, color = "grey35", size = 0.3) %<+% tip_df

  # Tip points
  p <- p + geom_tippoint(aes(color = species), size = 1.6, alpha = 0.9)

  # Tip labels (only for small trees)
  if (show_labels) {
    p <- p + geom_tiplab(aes(color = species), size = 1.7,
                         family = "", offset = 0.001, align = FALSE)
  }

  # Bootstrap support: show values < 95 on internal nodes
  p <- p + geom_nodelab(
    aes(label = ifelse(!isTip & !is.na(as.numeric(label)) &
                         as.numeric(label) < 95, label, "")),
    size = 1.6, color = "grey40", hjust = -0.15, vjust = -0.3
  )

  # Species colour scale (consistent palette, show all species in legend)
  p <- p +
    scale_color_manual(
      values = species_colors,
      name   = "Species",
      limits = names(species_colors),
      drop   = FALSE
    )

  # Title and axis
  p <- p +
    ggtitle(gsub("_", " ", dir_name)) +
    xlab("Substitutions per site") +
    theme_tree2(bgcolor = "white") +
    theme(
      plot.title      = element_text(size = 7, face = "bold.italic", family = ""),
      axis.text.x     = element_text(size = 5, family = ""),
      axis.title.x    = element_text(size = 5, family = ""),
      legend.title    = element_text(size = 5, face = "bold", family = ""),
      legend.text     = element_text(size = 5, family = ""),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "right",
      plot.margin     = margin(4, 4, 4, 4)
    )

  # Expand x axis right to make room for tip labels
  x_expand <- if (show_labels) 0.5 else 0.05
  p <- p + ggplot2::xlim(NA, max(fortify(tree)$x, na.rm = TRUE) * (1 + x_expand))

  # Save
  out_dir  <- dirname(treefile)
  height_in <- max(4, n_tips * 0.13)
  width_in  <- if (show_labels) 8 else 5

  for (ext in c("pdf", "png")) {
    out_path <- file.path(out_dir, glue("tree_plot.{ext}"))
    ggsave(out_path, plot = p,
           width = width_in, height = height_in,
           units = "in", dpi = 300, limitsize = FALSE)
  }

  cat(glue("  {dir_name}: {n_tips} tips → tree_plot.pdf/png\n"))
  invisible(p)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# First arg (if it looks like a directory path) sets the output dir;
# remaining args filter to specific subdirectories.
output_dir <- "output"
subdir_filter <- character(0)
if (length(args) >= 1) {
  if (dir.exists(args[1])) {
    output_dir <- args[1]
    subdir_filter <- args[-1]
  } else {
    subdir_filter <- args
  }
}

targets <- if (length(subdir_filter) > 0) {
  subdir_filter
} else {
  basename(list.dirs(output_dir, recursive = FALSE))
}

for (target in targets) {
  treefile <- file.path(output_dir, target, "tree.treefile")
  if (!file.exists(treefile)) {
    cat(glue("  SKIP {target}: tree.treefile not found\n"))
    next
  }
  plot_one_tree(treefile, target)
}
