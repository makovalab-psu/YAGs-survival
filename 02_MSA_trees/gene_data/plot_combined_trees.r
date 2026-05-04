# Combined phylogenetic tree plot with shared X scale
# Trees are stacked vertically, sharing the branch-length axis
#
# Usage: Rscript plot_combined_trees.r [output_prefix]
#   output_prefix: base name for output files (default: "tree_plots/combined_trees")

library(ggplot2)
library(ggtree)
library(treeio)
library(dplyr)
library(ape)
library(glue)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
output_prefix <- ifelse(length(args) >= 1, args[1], "tree_plots/combined_aligned")

# ---------------------------------------------------------------------------
# Helper functions (shared with plot_single_tree.r)
# ---------------------------------------------------------------------------

extract_species <- function(label) {
  species_map <- c(
    "PanPan" = "Bonobo",
    "PanTro" = "Chimpanzee",
    "HomSap" = "Human",
    "GorGor" = "Gorilla",
    "PonPyg" = "B. orangutan",
    "PonAbe" = "S. orangutan",
    "SymSyn" = "Siamang"
  )
  label <- gsub("^'|'$", "", label)
  prefix <- substr(label, 1, 6)
  return(ifelse(prefix %in% names(species_map), species_map[prefix], "Unknown"))
}

is_palindrome <- function(label) {
  label <- gsub("^'|'$", "", label)
  grepl("_Q[0-9]+(\\.[0-9]+)?\\.?[AB]_", label) ||
    grepl("LOC129530297", label) ||
    grepl("LOC129530299", label) ||
    grepl("LOC129530301", label)
}

is_vcx <- function(label) {
  label <- gsub("^'|'$", "", label)
  grepl("chrX", label)
}

extract_q_number <- function(label) {
  label <- gsub("^'|'$", "", label)
  m <- regexpr("_Q([0-9]+(\\.[0-9]+)?)\\.?[AB]_", label, perl = TRUE)
  if (m[1] == -1) return(NA)
  matched_str <- regmatches(label, m)
  cleaned <- gsub("^_|_$", "", matched_str)
  q_number <- gsub("([0-9])([AB])$", "\\1.\\2", cleaned)
  return(q_number)
}

extract_species_code <- function(label) {
  label <- gsub("^'|'$", "", label)
  strsplit(label, "_")[[1]][1]
}

extract_palindrome_base <- function(label) {
  label <- gsub("^'|'$", "", label)
  m <- regmatches(label, regexec("_(Q[0-9]+(?:\\.[0-9]+)?)\\.?[AB]_", label, perl = TRUE))[[1]]
  if (length(m) < 2) return(NA_character_)
  m[2]
}

extract_palindrome_arm <- function(label) {
  label <- gsub("^'|'$", "", label)
  m <- regmatches(label, regexec("_Q[0-9]+(?:\\.[0-9]+)?\\.?([AB])_", label, perl = TRUE))[[1]]
  if (length(m) < 2) return(NA_character_)
  m[2]
}

extract_position <- function(label) {
  label <- gsub("^'|'$", "", label)
  parts <- strsplit(label, "_")[[1]]
  if (length(parts) >= 4) as.integer(parts[4]) else NA_integer_
}

species_colors <- c(

  "Bonobo"               = "#fb2d32",
  "Chimpanzee"           = "#c3000e",
  "Gorilla"              = "#006e81",
  "Human"                = "#969696",
  "B. orangutan"  = "#f9cf44",
  "S. orangutan" = "#fd8a1a",
  "Siamang"              = "#0e205e"
)

# ---------------------------------------------------------------------------
# Custom scale transformation: compress 0 to breakpoint range
# ---------------------------------------------------------------------------
library(scales)

# Gene array intervals for array-membership annotation
arrays_df <- read.table("../40_arrays/output/gene_arrays.tsv",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Parameters for compression
compress_breakpoint       <- 0.22  # compress everything below this
compress_factor           <- 0.15  # compression ratio (0.15 = ~7x smaller)
right_compress_breakpoint <- 0.29  # compress everything above this
right_compress_factor     <- 0.3   # compression ratio for right side

# Transform function: compress values below left breakpoint and above right breakpoint
compress_trans <- function() {
  bp  <- compress_breakpoint
  cf  <- compress_factor
  rbp <- right_compress_breakpoint
  rcf <- right_compress_factor

  linear_at_rbp <- bp * cf + (rbp - bp)  # transformed value at right breakpoint

  trans_new(
    name = "compress",
    transform = function(x) {
      ifelse(x < bp,
        x * cf,
        ifelse(x <= rbp,
          bp * cf + (x - bp),
          linear_at_rbp + (x - rbp) * rcf
        )
      )
    },
    inverse = function(y) {
      left_threshold  <- bp * cf
      right_threshold <- linear_at_rbp
      ifelse(y < left_threshold,
        y / cf,
        ifelse(y <= right_threshold,
          bp + (y - left_threshold),
          rbp + (y - right_threshold) / rcf
        )
      )
    },
    domain = c(-Inf, Inf)
  )
}

# ---------------------------------------------------------------------------
# Root a tree on outgroup, midpoint the root branches, ladderise
# ---------------------------------------------------------------------------
root_tree <- function(tree, outgroup_pattern) {
  outgroup_tips <- grep(outgroup_pattern, tree$tip.label, value = TRUE)
  if (length(outgroup_tips) == 0) {
    warning("No outgroup tips found for pattern: ", outgroup_pattern)
    return(tree)
  }

  tree <- root(tree, outgroup = outgroup_tips,
               resolve.root = TRUE, edgelabel = TRUE)

  # Midpoint the two root branches
  root_node  <- length(tree$tip.label) + 1
  root_edges <- which(tree$edge[, 1] == root_node)
  if (length(root_edges) == 2) {
    total <- sum(tree$edge.length[root_edges])
    tree$edge.length[root_edges] <- total / 2
  }

  tree <- ladderize(tree, right = FALSE)
  return(tree)
}

# ---------------------------------------------------------------------------
# Build a single tree panel
#
#   config   – list with: file, outgroup, gene_name, display_name, tag, align_node
#   show_x   – whether to show the x-axis (TRUE for bottom panel only)
#   x_offset – amount to shift x coordinates (for node alignment)
#
# Returns a list: plot, n_tips, max_x, align_x
# ---------------------------------------------------------------------------
make_tree_panel <- function(config, show_x = FALSE, x_offset = 0, arrays_df = NULL) {

  tree <- read.tree(config$file)
  tree <- root_tree(tree, config$outgroup)

  tree_data <- fortify(tree)
  n_tips    <- length(tree$tip.label)

  # Get alignment node x position (before offset)
  align_x <- NA
  if (!is.null(config$align_node)) {
    align_row <- tree_data %>% filter(node == config$align_node)
    if (nrow(align_row) > 0) {
      align_x <- align_row$x[1]
    }
  }

  # Apply x offset to tree_data
  tree_data$x <- tree_data$x + x_offset

  # --- Tip metadata ---
  tip_species          <- sapply(tree$tip.label, extract_species)
  tip_species_code     <- sapply(tree$tip.label, extract_species_code)
  tip_palindrome       <- sapply(tree$tip.label, is_palindrome)
  tip_palindrome_base  <- sapply(tree$tip.label, extract_palindrome_base)
  tip_palindrome_arm   <- sapply(tree$tip.label, extract_palindrome_arm)
  tip_vcx_flag         <- sapply(tree$tip.label, is_vcx)
  if (config$gene_name != "VCY") tip_vcx_flag[] <- FALSE

  # Array membership: check if tip position falls within a gene array interval
  gene_arrays <- if (!is.null(arrays_df))
    arrays_df[arrays_df$gene_families == config$gene_name, ] else data.frame()

  # Species that have multiple distinct arrays for this gene → need labels
  multi_array_sp <- if (nrow(gene_arrays) > 0) {
    counts <- table(gene_arrays$species)
    names(counts[counts > 1])
  } else character(0)

  tip_array_id <- sapply(tree$tip.label, function(lbl) {
    if (nrow(gene_arrays) == 0) return(NA_character_)
    sp  <- extract_species_code(lbl)
    pos <- extract_position(lbl)
    if (is.na(pos)) return(NA_character_)
    idx <- which(gene_arrays$species == sp & gene_arrays$start <= pos & pos <= gene_arrays$end)
    if (length(idx) == 0) return(NA_character_)
    as.character(gene_arrays$array_id[idx[1]])
  })

  tip_df <- data.frame(
    label            = tree$tip.label,
    species          = tip_species,
    species_code     = tip_species_code,
    palindrome       = tip_palindrome,
    palindrome_base  = tip_palindrome_base,
    palindrome_arm   = tip_palindrome_arm,
    array_id         = tip_array_id,
    show_array_label = !is.na(tip_array_id) & (tip_species_code %in% multi_array_sp),
    stringsAsFactors = FALSE
  )

  # Merge tip metadata into tree_data
  tree_data <- tree_data %>%
    left_join(tip_df, by = "label")

  # --- Base tree (built from shifted tree_data) ---
  p <- ggplot(tree_data, aes(x = x, y = y)) +
    geom_tree(color = "grey40", size = 0.35) +
    theme_tree2(bgcolor = "white")

  # Array highlight: species-colored glow behind nodes within an array interval
  array_tips <- tree_data %>% filter(isTip & !is.na(array_id))
  cat("  [array debug]", config$gene_name, "- array_tips:", nrow(array_tips), "\n")
  if (nrow(array_tips) > 0) {
    # Loop per species: open circle (stroke only, no fill)
    for (sp in unique(na.omit(array_tips$species))) {
      sp_data <- array_tips[array_tips$species == sp, ]
      p <- p + geom_point(data = sp_data, aes(x = x, y = y),
                          shape = 21, size = 2.5, stroke = 0.42,
                          fill = NA, color = unname(species_colors[sp]), alpha = 0.8,
                          inherit.aes = FALSE)
    }
    # # Array ID labels for species with multiple arrays
    # label_tips <- array_tips %>% filter(show_array_label)
    # if (nrow(label_tips) > 0) {
    #   p <- p + geom_text(data = label_tips,
    #                      aes(x = x + 0.004, y = y, label = array_id),
    #                      size = 1.76, family = "Arial", hjust = 0, inherit.aes = FALSE)
    # }
  }

  # --- Palindrome A/B arches ---
  a_tips <- tree_data[!is.na(tree_data$palindrome_base) & tree_data$isTip &
                        tree_data$palindrome_arm == "A",
                      c("species_code", "palindrome_base", "x", "y")]
  b_tips <- tree_data[!is.na(tree_data$palindrome_base) & tree_data$isTip &
                        tree_data$palindrome_arm == "B",
                      c("species_code", "palindrome_base", "x", "y")]
  names(a_tips)[3:4] <- c("x1", "y1")
  names(b_tips)[3:4] <- c("x2", "y2")
  if (nrow(a_tips) > 0 && nrow(b_tips) > 0) {
    pair_df <- merge(a_tips, b_tips, by = c("species_code", "palindrome_base"))
    if (nrow(pair_df) > 0) {
      # Cubic Bezier: control points pulled right → ")" shape
      arch_bulge <- diff(range(tree_data$x[tree_data$isTip], na.rm = TRUE)) * 0.05
      arch_pts <- do.call(rbind, lapply(seq_len(nrow(pair_df)), function(i) {
        x1 <- pair_df$x1[i]; y1 <- pair_df$y1[i]
        x2 <- pair_df$x2[i]; y2 <- pair_df$y2[i]
        xc <- max(x1, x2) + arch_bulge   # both control points to the right
        t  <- seq(0, 1, length.out = 40)
        data.frame(
          x     = (1-t)^3*x1 + 3*(1-t)^2*t*xc + 3*(1-t)*t^2*xc + t^3*x2,
          y     = (1-t)^3*y1 + 3*(1-t)^2*t*y1  + 3*(1-t)*t^2*y2  + t^3*y2,
          group = i
        )
      }))
      if (nrow(arch_pts) > 0)
        p <- p + geom_path(data = arch_pts,
                           aes(x = x, y = y, group = group),
                           color = "grey50", linewidth = 0.14, alpha = 0.7,
                           inherit.aes = FALSE)
    }
  }

  # Tip points
  p <- p +
    geom_tippoint(aes(color = species, shape = palindrome),
                  size = 1.4, alpha = 0.9) +
    scale_shape_manual(
      values = c("TRUE" = 17, "FALSE" = 16),
      name   = "Palindrome",
      labels = c("No", "Yes"),
      limits = c("FALSE", "TRUE"),
      drop   = FALSE
    )

  # --- Clade labels from tag ---
  if (config$tag != "") {
    tags <- strsplit(config$tag, split = ",")[[1]]
    for (tg in tags) {
      matched_tips <- grep(tg, tree$tip.label)
      if (length(matched_tips) >= 2) {
        mrca_node <- getMRCA(tree, matched_tips)
        # Get the clade data for manual label positioning
        clade_data <- tree_data %>% filter(node == mrca_node)
        if (nrow(clade_data) > 0) {
          # Find y range of tips in clade
          clade_tips <- tree_data %>%
            filter(isTip, node %in% offspring(tree, mrca_node, tiponly = TRUE))
          if (nrow(clade_tips) > 0) {
            ymin <- min(clade_tips$y)
            ymax <- max(clade_tips$y)
            ymid <- mean(c(ymin, ymax))
            xpos <- max(clade_tips$x) + 0.008
            # Add vertical bar
            p <- p + annotate("segment", x = xpos, xend = xpos,
                              y = ymin, yend = ymax,
                              color = "black", size = 0.28)
            # Add label text
            p <- p + annotate("text", x = xpos + 0.005, y = ymid,
                              label = gsub("_", "", tg),
                              hjust = 0, size = 1.76, fontface = "italic", family = "Arial")
          }
        }
      }
    }
  }

  # --- VCY-specific: mark chrX tips with "X" ---
  if (config$gene_name == "VCY") {
    vcx_labels <- names(tip_vcx_flag)[tip_vcx_flag]
    if (length(vcx_labels) > 0) {
      vcx_df <- tree_data %>%
        filter(isTip & label %in% vcx_labels) %>%
        mutate(label_x = x + 0.004)
      p <- p + geom_text(data = vcx_df,
                         aes(x = label_x, y = y, label = "X"),
                         size = 1.76, family = "Arial", inherit.aes = FALSE)
    }
  }

  # --- Show internal node numbers (for alignment selection) ---
  if (isTRUE(config$show_node_labels)) {
    internal_nodes <- tree_data %>% filter(!isTip)
    p <- p + geom_text(data = internal_nodes,
                       aes(x = x, y = y, label = node),
                       size = 1.76, color = "red", fontface = "bold", family = "Arial",
                       hjust = -0.2, vjust = -0.3,
                       inherit.aes = FALSE)
  }

  # --- Species color scale (consistent across all panels) ---
  p <- p +
    scale_color_manual(
      values = species_colors,
      name   = "Species",
      limits = names(species_colors),
      drop   = FALSE
    )



  # --- Gene-name label inside panel, top-left corner ---
  p <- p +
    annotate("text", x = -Inf, y = Inf,
             label = config$display_name,
             hjust = -0.15, vjust = 3.0,
             fontface = "bold.italic", size = 1.76, family = "Arial") +
    ylab(NULL) +
    theme(
      axis.title.y = element_blank(),
      axis.text.x  = element_text(size = 5, family = "Arial"),
      axis.title.x = element_text(size = 5, family = "Arial"),
      plot.margin  = margin(2, 5, 2, 2)
    )

  # --- X-axis: show only on the bottom panel ---
  if (!show_x) {
    p <- p + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x  = element_blank()
    )
  } else {
    p <- p + xlab("Substitutions per site")
  }

  list(plot = p, n_tips = n_tips, max_x = max(tree_data$x, na.rm = TRUE), align_x = align_x)
}

# ---------------------------------------------------------------------------
# Tree configurations (order = top to bottom in the combined plot)
# ---------------------------------------------------------------------------
# Set to TRUE to show node numbers for alignment selection
show_node_numbers <- FALSE

trees_config <- list(
  list(file = "output/DAZ_DAZL_RBD.treefile",
       outgroup = "DAZL",
       gene_name = "DAZ", display_name = "DAZ", tag = "DAZL",
       show_node_labels = show_node_numbers,
       align_node = 79),

  list(file = "output/CDY_CDYL.treefile",
       outgroup = "CDYL",
       gene_name = "CDY+CDYL", display_name = "CDY", tag = "CDYL",
       show_node_labels = show_node_numbers,
       align_node = 65),

  list(file = "output/RBMY_RBMX.treefile",
       outgroup = "RBMX",
       gene_name = "RBMY", display_name = "RBMY", tag = "RBMX",
       show_node_labels = show_node_numbers,
       align_node = 86)
)

# ---------------------------------------------------------------------------
# Build all panels (two-pass for node alignment)
# ---------------------------------------------------------------------------
n_trees <- length(trees_config)

# First pass: get alignment node x positions
cat("First pass: collecting alignment node positions...\n")
align_positions <- numeric(n_trees)
for (i in seq_along(trees_config)) {
  result <- make_tree_panel(trees_config[[i]], show_x = FALSE, x_offset = 0)
  align_positions[i] <- result$align_x
  cat("  ", trees_config[[i]]$display_name, ": align node",
      trees_config[[i]]$align_node, "at x =", result$align_x, "\n")
}

# Calculate offsets to align all trees at the max alignment x
target_align_x <- max(align_positions, na.rm = TRUE)
x_offsets <- target_align_x - align_positions
cat("\nAlignment target x:", target_align_x, "\n")
cat("X offsets:", paste(round(x_offsets, 4), collapse = ", "), "\n\n")

# Second pass: build plots with offsets applied
results <- vector("list", n_trees)
for (i in seq_along(trees_config)) {
  cat("Building", trees_config[[i]]$display_name, "with offset", round(x_offsets[i], 4), "...\n")
  results[[i]] <- make_tree_panel(trees_config[[i]],
                                  show_x = (i == n_trees),
                                  x_offset = x_offsets[i],
                                  arrays_df = arrays_df)
}

# ---------------------------------------------------------------------------
# Shared X range (accounting for offsets - some trees may have negative x)
# ---------------------------------------------------------------------------
global_max_x <- max(sapply(results, function(r) r$max_x))
global_min_x <- min(x_offsets)  # leftmost offset (could be negative or 0)
tip_counts   <- sapply(results, function(r) r$n_tips)

cat("\nGlobal x range:", global_min_x, "to", global_max_x, "\n")
cat("Tip counts:", paste(tip_counts, collapse = ", "), "\n")
cat("Total tips:", sum(tip_counts), "\n")

# Vertical reference lines at key x positions
vline_positions <- c(0, 0.1, 0.2, 0.3)

# Apply compressed x-scale and suppress per-panel legends
plots <- lapply(results, function(r) {
  r$plot +
    geom_vline(xintercept = vline_positions, color = "grey70",
               linetype = "solid", size = 0.14) +
    # Highlight alignment position
    geom_vline(xintercept = target_align_x, color = "steelblue",
               linetype = "solid", size = 0.42, alpha = 0.4) +
    scale_x_continuous(
      trans = compress_trans(),
      breaks = c(0, 0.1, 0.2, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29),
      labels = c("0", "0.1", "0.2", "||", "0.23", "0.24", "0.25", "0.26", "0.27", "0.28", "0.29")
    ) +
    coord_cartesian(xlim = c(global_min_x - 0.01, global_max_x )) +
    theme(legend.position = "none")
})

# ---------------------------------------------------------------------------
# Build comprehensive figure legend
# ---------------------------------------------------------------------------
n_sp     <- length(species_colors)
sp_names <- names(species_colors)
sp_hex   <- unname(species_colors)
sp_ys    <- seq(6.5, by = -1, length.out = n_sp)   # one row per species

# Bezier ")" arch preview for legend
t_leg    <- seq(0, 1, length.out = 40)
arch_leg <- data.frame(
  x = (1-t_leg)^3*0 + 3*(1-t_leg)^2*t_leg*0.28 +
      3*(1-t_leg)*t_leg^2*0.28 + t_leg^3*0,
  y = (1-t_leg)^3*0.38 + 3*(1-t_leg)^2*t_leg*0.38 +
      3*(1-t_leg)*t_leg^2*(-0.38) + t_leg^3*(-0.38)
)

figure_legend <- ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(4, 10, 4, 10)) +

  # ── Species ──────────────────────────────────────────────────────────────
  annotate("text", x = 0.05, y = 7.8, label = "Species",
           hjust = 0, fontface = "bold", size = 1.76, family = "Arial") +
  geom_point(data = data.frame(x = 0.25, y = sp_ys),
             aes(x = x, y = y), color = sp_hex, shape = 16, size = 1.4,
             inherit.aes = FALSE) +
  geom_text(data = data.frame(x = 0.55, y = sp_ys, label = sp_names),
            aes(x = x, y = y, label = label),
            hjust = 0, size = 1.76, family = "Arial", inherit.aes = FALSE) +

  # ── Node type ─────────────────────────────────────────────────────────────
  annotate("text", x = 3.4, y = 7.8, label = "Node type",
           hjust = 0, fontface = "bold", size = 1.76, family = "Arial") +
  annotate("point", x = 3.6, y = 6.5, shape = 16, size = 1.4, color = "grey35") +
  annotate("text",  x = 3.9, y = 6.5, label = "Non-palindromic copy",
           hjust = 0, size = 1.76, family = "Arial") +
  annotate("point", x = 3.6, y = 5.3, shape = 17, size = 1.4, color = "grey35") +
  annotate("text",  x = 3.9, y = 5.3, label = "Palindromic copy",
           hjust = 0, size = 1.76, family = "Arial") +
  annotate("point", x = 3.6, y = 4.1, shape = 21, size = 2.1,
           fill = NA, color = "grey35", stroke = 0.46) +
  annotate("text",  x = 3.9, y = 4.1, label = "Gene copy within array",
           hjust = 0, size = 1.76, family = "Arial") +

  # ── Palindrome arm pair ───────────────────────────────────────────────────
  annotate("text", x = 6.8, y = 7.8, label = "Palindrome arms",
           hjust = 0, fontface = "bold", size = 1.76, family = "Arial") +
  annotate("point", x = 7.0, y = 6.0, shape = 16, size = 1.1, color = "grey35") +
  annotate("point", x = 7.0, y = 5.0, shape = 16, size = 1.1, color = "grey35") +
  geom_path(data = transform(arch_leg, x = x + 7.0, y = y + 5.5),
            aes(x = x, y = y), color = "grey50", linewidth = 0.18,
            inherit.aes = FALSE) +
  annotate("text", x = 7.45, y = 5.5, label = "Palindrome arm pair (A / B)",
           hjust = 0, size = 1.76, family = "Arial") +
  # ── Reference lines ───────────────────────────────────────────────────────
  annotate("text", x = 10.2, y = 7.8, label = "Reference lines",
           hjust = 0, fontface = "bold", size = 1.76, family = "Arial") +
  annotate("segment", x = 10.2, xend = 10.7, y = 6.5, yend = 6.5,
           color = "steelblue", linewidth = 0.46, alpha = 0.5) +
  annotate("text",    x = 10.9, y = 6.5, label = "Outgroup alignment",
           hjust = 0, size = 1.76, family = "Arial") +
  annotate("segment", x = 10.2, xend = 10.7, y = 5.3, yend = 5.3,
           color = "grey70", linewidth = 0.21) +
  annotate("text",    x = 10.9, y = 5.3, label = "0.1 subs/site intervals",
           hjust = 0, size = 1.76, family = "Arial") +

  xlim(-0.1, 13.5) + ylim(0, 8.3)

# ---------------------------------------------------------------------------
# Build TSPY tree (separate, spans full height on the right)
# ---------------------------------------------------------------------------
cat("\nBuilding TSPY tree...\n")
tspy_config <- list(
  file = "output/TSPY.treefile",
  outgroup = "SymSyn",
  gene_name = "TSPY",
  display_name = "TSPY",
  tag = "",
  show_node_labels = FALSE,
  align_node = NULL
)

tspy_result <- make_tree_panel(tspy_config, show_x = TRUE, x_offset = 0, arrays_df = arrays_df)
cat("TSPY tips:", tspy_result$n_tips, "\n")

# Apply styling to TSPY (no vertical lines, normal x-axis, no legend)
tspy_plot <- tspy_result$plot +
  theme(legend.position = "none",
        plot.margin = margin(0, 2, 0, 2)) +
  scale_y_continuous(expand = expansion(add = 2))

# ---------------------------------------------------------------------------
# Build third column: BPY2, HSFY, VCY (stacked)
# ---------------------------------------------------------------------------
cat("\nBuilding third column trees (BPY2, HSFY, VCY)...\n")

col3_configs <- list(
  list(file = "output/BPY2.treefile",
       outgroup = "GorGor",
       gene_name = "BPY2", display_name = "BPY2", tag = "",
       show_node_labels = FALSE, align_node = 9),

  list(file = "output/HSFY.treefile",
       outgroup = "SymSyn",
       gene_name = "HSFY", display_name = "HSFY", tag = "",
       show_node_labels = FALSE, align_node = 24),

  list(file = "output/VCY_VCX.treefile",
       outgroup = "SymSyn_chrY",
       gene_name = "VCY", display_name = "VCY", tag = "",
       show_node_labels = FALSE, align_node = NULL)
)

# First pass for BPY2 and HSFY alignment
col3_align_positions <- numeric(2)
for (i in 1:2) {
  result <- make_tree_panel(col3_configs[[i]], show_x = FALSE, x_offset = 0)
  col3_align_positions[i] <- result$align_x
  cat("  ", col3_configs[[i]]$display_name, ": align node",
      col3_configs[[i]]$align_node, "at x =", result$align_x, "\n")
}

# Calculate offsets for BPY2/HSFY alignment
col3_target_align_x <- max(col3_align_positions, na.rm = TRUE)
col3_x_offsets <- col3_target_align_x - col3_align_positions
cat("  Column 3 alignment target x:", col3_target_align_x, "\n")
cat("  Column 3 x offsets (BPY2, HSFY):", paste(round(col3_x_offsets, 4), collapse = ", "), "\n")

# Build all col3 trees with appropriate offsets (all with show_x = TRUE to preserve axis)
col3_results <- list(
  make_tree_panel(col3_configs[[1]], show_x = TRUE, x_offset = col3_x_offsets[1], arrays_df = arrays_df),
  make_tree_panel(col3_configs[[2]], show_x = TRUE, x_offset = col3_x_offsets[2], arrays_df = arrays_df),
  make_tree_panel(col3_configs[[3]], show_x = TRUE, x_offset = 0,                 arrays_df = arrays_df)
)

for (i in seq_along(col3_configs)) {
  cat("  ", col3_configs[[i]]$display_name, ":", col3_results[[i]]$n_tips, "tips\n")
}

col3_tip_counts <- sapply(col3_results, function(r) r$n_tips)

# BPY2 and HSFY share x-axis (with alignment), VCY has its own
bpy2_hsfy_max_x <- max(col3_results[[1]]$max_x, col3_results[[2]]$max_x)
bpy2_hsfy_min_x <- min(col3_x_offsets)

# BPY2 offset for axis label correction
bpy2_offset <- col3_x_offsets[1]

# Shared xlim for visual alignment (in shifted coordinate space)
shared_xlim <- c(0, bpy2_hsfy_max_x * 1.05)

col3_plots <- list(
  # BPY2 - x-axis on bottom, labels show BPY2's original scale (starting at 0)
  col3_results[[1]]$plot +
    scale_x_continuous(
      labels = function(x) round(x - bpy2_offset, 3)
    ) +
    coord_cartesian(xlim = c(bpy2_offset, shared_xlim[2])) +
    xlab("") +
    theme(legend.position = "none",
          axis.line.x = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black"),
          axis.text.x = element_text(size = 5, family = "Arial")),
  # HSFY - x-axis on bottom, labels show HSFY's original scale (no offset)
  col3_results[[2]]$plot +
    coord_cartesian(xlim = shared_xlim) +
    xlab("") +
    theme(legend.position = "none",
          axis.line.x = element_line(color = "black"),
          axis.ticks.x = element_line(color = "black")),
  # VCY - own axis, shows x-axis (bottom)
  col3_results[[3]]$plot +
    theme(legend.position = "none")
)

# ---------------------------------------------------------------------------
# Combine with patchwork
# ---------------------------------------------------------------------------
tree_stack   <- wrap_plots(plots, ncol = 1, heights = tip_counts)
col3_stack   <- wrap_plots(col3_plots, ncol = 1, heights = col3_tip_counts)
total_height <- max(sum(tip_counts) * 0.08 + 2, 10)

tree_panel <- (tree_stack | tspy_plot | col3_stack) +
  plot_layout(widths = c(8, 5, 5)) +
  plot_annotation(
    title = "Y-Chromosome Gene Families with Autosomal/X Outgroups",
    theme = theme(
      plot.title = element_text(size = 5, face = "bold", hjust = 0.5, family = "Arial")
    )
  )

legend_height <- 2.5  # inches
combined <- tree_panel / figure_legend +
  plot_layout(heights = c(total_height, legend_height))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
total_width <- 180  # mm
# save_height <- (total_height + legend_height) * 9  # scale to 180mm width proportionally
save_height <- 225

dir.create(dirname(output_prefix), showWarnings = FALSE, recursive = TRUE)

ggsave(glue("{output_prefix}.pdf"), plot = combined,
       width = total_width, height = save_height, units = "mm",
       device = cairo_pdf, limitsize = FALSE)
ggsave(glue("{output_prefix}.png"), plot = combined,
       width = total_width, height = save_height, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("\nCombined plot saved:\n")
cat(glue("  {output_prefix}.pdf\n"))
cat(glue("  {output_prefix}.png\n"))
cat(glue("  Dimensions: {total_width} x {round(save_height)} mm\n"))
