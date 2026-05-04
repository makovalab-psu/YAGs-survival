import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import json
from collections import defaultdict


mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif']


# Load the JSON data
with open('all_RBMY_merge_PonPyg.json', 'r') as f:
    data = json.load(f)

# Extract unique species, genes, and isoforms
species_list = []
gene_list = []
all_isoforms = set()

species_to_label = {
    'HomSap': 'Human',
    "GorGor": "Gorilla",
    "PanPan": "Bonobo",
    "PanTro": "Chimpanzee",
    "PonAbe": "S. orangutan",
    "PonPyg": "B. orangutan"
}

for entry in data:
    # Skip RBMYB gene
    if entry['gene'] == 'RBMYB':
        continue
        # pass

    if entry['species'] not in species_list:
        species_list.append(entry['species'])
    if entry['gene'] not in gene_list:
        gene_list.append(entry['gene'])
    for isoform_data in entry['isoforms']:
        isoform_name = isoform_data['isoform']
        if isoform_name:  # Skip empty isoform names
            all_isoforms.add(isoform_name)

n_rows = len(gene_list)
n_cols = len(species_list)

print(f"Genes (rows): {gene_list}")
print(f"Species (columns): {species_list}")
print(f"Number of unique isoforms: {len(all_isoforms)}")

# Create a color map for isoforms within each gene
gene_isoforms = defaultdict(set)
for entry in data:
    # Skip RBMYB gene
    if entry['gene'] == 'RBMYB':
        continue
        # pass

    gene = entry['gene']
    for isoform_data in entry['isoforms']:
        isoform_name = isoform_data['isoform']
        if isoform_name:
            gene_isoforms[gene].add(isoform_name)

# Assign colors to isoforms for each gene
gene_color_maps = {}
colormap = plt.cm.tab20
for gene_idx, gene in enumerate(gene_list):
    isoforms = sorted(list(gene_isoforms[gene]))
    colors = {}
    for i, isoform in enumerate(isoforms):
        colors[isoform] = colormap(i % 20)
    gene_color_maps[gene] = colors

# Organize data by gene and species (swapped from before)
data_structure = defaultdict(lambda: defaultdict(list))
for entry in data:
    # Skip RBMYB gene
    if entry['gene'] == 'RBMYB':
        continue
        # pass

    species = entry['species']
    gene = entry['gene']
    data_structure[gene][species] = entry['isoforms']

# Find max ord value for x-axis
max_ord = 0
for entry in data:
    for isoform_data in entry['isoforms']:
        for count_entry in isoform_data['counts']:
            max_ord = max(max_ord, count_entry['ord'])

# Create figure
MM_TO_INCH = 25.4

fig, axes = plt.subplots(n_rows, n_cols, figsize=(170/MM_TO_INCH, 130/MM_TO_INCH), sharex=True)

# Ensure axes is 2D
if n_rows == 1:
    axes = axes.reshape(1, -1)
if n_cols == 1:
    axes = axes.reshape(-1, 1)

# Extract total counts from JSON and find max per row
subplot_totals = {}
row_max_totals = []

for row_idx, gene in enumerate(gene_list):
    row_max = 0
    for col_idx, species in enumerate(species_list):
        # Find the corresponding entry in the original data
        total_counts = 0
        for entry in data:
            if entry['gene'] == gene and entry['species'] == species:
                total_counts = entry['total']
                break

        subplot_totals[(gene, species)] = total_counts
        row_max = max(row_max, total_counts)

    row_max_totals.append(row_max)

# Set Y-axis max to 100% for percentage plots
row_max_values = [100] * len(gene_list)


def index_to_letter(index):
    label = []
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        label.append(chr(ord('A') + remainder))
    return ''.join(reversed(label))

# Plot data - keep original x positions
x = np.arange(1, max_ord + 1)

for row_idx, gene in enumerate(gene_list):
    for col_idx, species in enumerate(species_list):
        ax = axes[row_idx, col_idx]
        ax.tick_params(axis='both', which='major', labelsize=5)

        isoforms = data_structure[gene][species]

        if not isoforms:
            # Empty plot
            print(gene)
            text = "Gene family\nabsent on the Y\n in this species" 
            if gene == "RBMY" or gene == "RBMYB":
                text = "Low covereage"
            ax.text(0.5, 0.5, text, ha='center', va='center',
                   transform=ax.transAxes, fontsize=5, color='gray')
        else:
            # Prepare data for stacking
            # First pass: collect all heights per isoform
            isoform_heights = {}
            for isoform_data in isoforms:
                isoform_name = isoform_data['isoform']
                if not isoform_name:
                    isoform_name = 'Unknown'

                counts_dict = {}
                for count_entry in isoform_data['counts']:
                    counts_dict[count_entry['ord']] = count_entry['count']

                heights = np.array([counts_dict.get(i, 0) for i in range(1, max_ord + 1)])
                isoform_heights[isoform_name] = heights

            # Calculate total plotted counts across all positions (single scalar)
            total_plotted_sum = sum(h.sum() for h in isoform_heights.values())

            # Get pre-calculated total counts for gray bar background
            total_counts = subplot_totals[(gene, species)]

            # Add light gray bar behind the plot area showing total counts (scaled within row)
            if total_counts > 0 and row_max_totals[row_idx] > 0:
                gray_bar_height = (total_counts / row_max_totals[row_idx]) * 100
                # Draw gray bar first so it appears behind everything else
                ax.axhspan(0, gray_bar_height, color='lightgray', alpha=0.3, zorder=0)

            # Stack the bars normalized so all bars across the subplot sum to 100%
            bottom = np.zeros(max_ord)
            for isoform_name, heights in isoform_heights.items():
                percentage_heights = (heights / total_plotted_sum) * 100 if total_plotted_sum > 0 else np.zeros_like(heights, dtype=float)

                color = gene_color_maps[gene].get(isoform_name, 'gray')
                ax.bar(x, percentage_heights, bottom=bottom, color=color,
                      label=isoform_name, width=0.8)
                bottom += percentage_heights

        # Set Y-axis limit
        ax.set_ylim(0, row_max_values[row_idx])

        if col_idx != 0:
            ax.set_yticklabels([])

        # Add right Y-axis for count values
        ax2 = ax.twinx()

        ax2.set_ylim(0, row_max_totals[row_idx])
        # Add label only to the rightmost plot in each row
        if col_idx == n_cols - 1:
            ax2.set_ylabel('Total Counts', fontsize=5, rotation=270, labelpad=15)
            ax2.tick_params(axis='y', labelsize=5)
        else:
            ax2.set_ylabel('')
            ax2.yaxis.set_major_formatter(plt.NullFormatter())
        ax2.spines['top'].set_visible(False)
        ax2.spines['left'].set_visible(False)

        # Styling for main axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', alpha=0.3, zorder=1)

        # Labels
        if col_idx == 0:
            ax.set_ylabel(f'{gene}\n(%)', fontsize=7, fontweight='bold', fontstyle="italic")

        if row_idx == 0:
            ax.set_title(species_to_label[species], fontsize=7, fontweight='bold')

        if row_idx == n_rows - 1:
            # Set x-ticks for regular positions only
            tick_positions = x[::max(1, len(x)//5)]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels([index_to_letter(position) for position in tick_positions])

# Add legends for each row (genes and their isoforms)
legend_y_positions = np.linspace(0.95, 0.05, n_rows)
for row_idx, gene in enumerate(gene_list):
    # Get all isoforms for this gene
    isoforms = sorted(list(gene_isoforms[gene]))
    
    # assert False
    if isoforms:
        handles = [plt.Rectangle((0,0),1,1, color=gene_color_maps[gene][iso])
                  for iso in isoforms]

        # Place legend on the right side of the figure, aligned with the row
        isoforms = [" ".join(x.split("_")[1:]) for x in isoforms]
        isoforms = [x.replace("isoform", "Isoform") for x in isoforms]
        fig.legend(handles, isoforms, loc='center left',
                  bbox_to_anchor=(1.01, legend_y_positions[row_idx]),
                  fontsize=5, title=gene, title_fontsize=5, frameon=True)

plt.tight_layout()
fig.supxlabel('Each bar represents a unique sequence isoform', fontsize=7, ha='center', y=0.04)
plt.subplots_adjust(bottom=0.10, hspace=0.3, wspace=0.5)
plt.savefig("isoforms.pdf")
