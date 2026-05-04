Assignment is to create a Sankey diagram (ribbon plot) with 3 categories (vertically arranged)
connected with ribbons (horizontally). The 3 categories are, sequence isoform + structural isoform,
protein sequence, cluster. For this we have 3 input files. The first one is a collection of tabs
in a google sheet located at 
    https://docs.google.com/spreadsheets/d/1dRLy4ZM9dpLhP3PpDfIMemmKg9IVGwB1l0WPeu18_vY/edit?gid=0#gid=0
there is one sheet per Species/gene combination. There are 6 species and 6 genes (not all species
have all genes) and these should be reflected as groupings (essentially each species x gene combination
is it's own ribbon plot. All ribbon plots are placed above each other.

The plot will be pretty dense so lets start with small fonts thin lines.
Order first by species:
    1. Bonobo, 2. Chimpanzee, 3. Human, 4. Gorilla, 5. Bornean orangutan, 6 Sumatran orangutan.

Then order by genes (within each species):
    1. BPY2, 2. CDY, 3. HSFY, 4. RBMY, 5. TSPY, 6. VCY

First column will be all separate sequence + structural isoforms

``` python
    gene_families = [

        {
            'name': 'BPY2',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/BPY2/BPY2_transcripts.gff',
            'isoform_list': [
                ('BPY2_isoform_13', '#4173b1'),
                ('BPY2_isoform_14', '#b4c6e6'),
                ('BPY2_isoform_25', "#ec892e"),
                ('BPY2_isoform_5', "#f3bf7f"),
                ('BPY2_isoform_6', "#539d38"),
                ('BPY2_isoform_7', "#a8dd8f"),
            ],
            'output_file': 'BPY2_transcripts.png'
        },
        {
            'name': 'CDY',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/CDY/CDY_transcripts.gff',
            'isoform_list': [
                ('CDY_isoform_1', '#4173b1'),
                ('CDY_isoform_2', '#b4c6e6'),
            ],
            'output_file': 'CDY_transcripts.png'
        },
        {
            'name': 'HSFY',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/HSFY/HSFY_transcripts.gff',
            'isoform_list': [
                ('HSFY_isoform_2', '#4173b1'),
                ('HSFY_isoform_3', '#b4c6e6'),
            ],
            'output_file': 'HSFY_transcripts.png'
        },
        {
            'name': 'RBMY',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/RBMY/RBMY_transcripts.gff',
            'isoform_list':  [
                ('RBMY_isoform_10', '#4173b1'),
                ('RBMY_isoform_14', '#b4c6e6'),
                ('RBMY_isoform_20', "#ec892e"),
                ('RBMY_isoform_29', "#f3bf7f"),
                ('RBMY_isoform_35', "#539d38"),
                ('RBMY_isoform_36', "#a8dd8f"),
                ('RBMY_isoform_38', "#c24030"),
                ('RBMY_isoform_44', "#ef9f98"),
                ('RBMY_isoform_46', "#8e69ba"),
                ('RBMY_isoform_5', "#c2b1d3"),
                ('RBMY_isoform_9', "#83594d")
            ],
            'output_file': 'RBMY_transcripts.png'
        },
        {
            'name': 'TSPY',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/TSPY/TSPY_transcripts.gff',
            'isoform_list':  [
                ('TSPY_isoform_12', '#4173b1'),
                ('TSPY_isoform_25', '#b4c6e6'),
                ('TSPY_isoform_28', "#ec892e"),
                ('TSPY_isoform_32', "#f3bf7f"),
                ('TSPY_isoform_36', "#539d38"),
                ('TSPY_isoform_37', "#a8dd8f"),
                ('TSPY_isoform_40', "#c24030")
            ],
            'output_file': 'TSPY_transcripts.png'
        },
        {
            'name': 'VCY',
            'gff_file': '/Users/kxp5629/proj/Y/src/27_plot_gff_isoforms/compare_isoforms/VCY/VCY_transcripts.gff',
            'isoform_list': [
                ('VCY_isoform_1', '#4173b1'),
                ('VCY_isoform_2', '#b4c6e6'),
            ],
            'output_file': 'VCY_transcripts.png'
        }
    ]
```


Second column is protein sequence. Here we will see some collapses (synonymous mutations). For this
we need to YAGS_alfafold/data For each gene there is a redundant_names.tsv file where the first column is the representative sequence all other columns are sequences with same protein sequence. Here combine bars together.

The final third column will be clusters. Again for each gene find the clust.tsv file. First column is the representative sequence name, second column is another sequence in the same cluster. (the forat is slightlly different from reduntant_names.tsv where one line can have multiple columns - in clust.tsv there are multiple line per one cluster.

For the final ribbon connecting proteins and clusters group across all species, keep gene order.