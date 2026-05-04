#!/usr/bin/env python3
"""Sankey ribbon: sequence+structural isoform → protein sequence → structural cluster.
Two-column layout, 180 mm wide, Arial font, sizes 5/6/7 pt.
"""

import colorsys
import re
from collections import defaultdict
from datetime import datetime
from pathlib import Path

def hex_to_rgb(h: str) -> tuple:
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))

def rgb_to_hex(r: float, g: float, b: float) -> str:
    return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'

def blend_colors(colors_weights: list) -> str:
    """Weighted average of hex colors; colors_weights = [(hex, weight), ...]."""
    total = sum(w for _, w in colors_weights)
    if not total:
        return '#aaaaaa'
    r = sum(hex_to_rgb(c)[0] * w for c, w in colors_weights) / total
    g = sum(hex_to_rgb(c)[1] * w for c, w in colors_weights) / total
    b = sum(hex_to_rgb(c)[2] * w for c, w in colors_weights) / total
    return rgb_to_hex(r, g, b)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

# ─── Font ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})
FS_SM = 5   # small  – bar labels, cluster labels on thin bars
FS_MD = 6   # medium – species labels, cluster labels
FS_LG = 7   # large  – gene labels, column headers

# ─── Paths ────────────────────────────────────────────────────────────────────
DATA_DIR   = Path('/Users/kxp5629/proj/Y/YAGs-survival/seq_struct_prot_clust_ribbon/YAGS_alfafold/data')
OUTPUT_DIR = Path('/Users/kxp5629/proj/Y/YAGs-survival/seq_struct_prot_clust_ribbon/output')
OUTPUT_DIR.mkdir(exist_ok=True)

# ─── Species / gene ordering ──────────────────────────────────────────────────
SPECIES_ORDER = ['PanPan', 'PanTro', 'HomSap', 'GorGor', 'PonAbe', 'PonPyg']
SPECIES_LABELS = {
    'PanPan': 'Bonobo',         'PanTro': 'Chimpanzee',
    'HomSap': 'Human',          'GorGor': 'Gorilla',
    'PonAbe': 'B. orangutan',   'PonPyg': 'S. orangutan',
}

# Two-column gene arrangement
LEFT_GENES  = ['BPY2', 'CDY', 'HSFY', 'VCY']
RIGHT_GENES = ['RBMY', 'TSPY', 'DAZ']
ALL_GENES   = LEFT_GENES + RIGHT_GENES

# ─── Isoform colours ──────────────────────────────────────────────────────────
ISOFORM_COLORS = {
    'BPY2_isoform_13': '#4173b1', 'BPY2_isoform_14': '#b4c6e6',
    'BPY2_isoform_25': '#ec892e', 'BPY2_isoform_5':  '#f3bf7f',
    'BPY2_isoform_6':  '#539d38', 'BPY2_isoform_7':  '#a8dd8f',
    'CDY_isoform_1':   '#4173b1', 'CDY_isoform_2':   '#b4c6e6',
    'HSFY_isoform_2':  '#4173b1', 'HSFY_isoform_3':  '#b4c6e6',
    'RBMY_isoform_10': '#4173b1', 'RBMY_isoform_14': '#b4c6e6',
    'RBMY_isoform_20': '#ec892e', 'RBMY_isoform_29': '#f3bf7f',
    'RBMY_isoform_35': '#539d38', 'RBMY_isoform_36': '#a8dd8f',
    'RBMY_isoform_38': '#c24030', 'RBMY_isoform_44': '#ef9f98',
    'RBMY_isoform_46': '#8e69ba', 'RBMY_isoform_5':  '#c2b1d3',
    'RBMY_isoform_9':  '#83594d',
    'TSPY_isoform_12': '#4173b1', 'TSPY_isoform_25': '#b4c6e6',
    'TSPY_isoform_28': '#ec892e', 'TSPY_isoform_32': '#f3bf7f',
    'TSPY_isoform_36': '#539d38', 'TSPY_isoform_37': '#a8dd8f',
    'TSPY_isoform_40': '#c24030',
    'VCY_isoform_1':   '#4173b1', 'VCY_isoform_2':   '#b4c6e6',
}
DEFAULT_COLOR = '#aaaaaa'
# Per-family (hue, saturation, L_dark, L_light) derived from existing colour pairs
# Used for genes without ISOFORM_COLORS entries (e.g. DAZ)
def _palette_from_pair(dark: str, light: str) -> tuple:
    h, l_d, s = colorsys.rgb_to_hls(*hex_to_rgb(dark))
    _, l_l, _ = colorsys.rgb_to_hls(*hex_to_rgb(light))
    return (h, s, l_d, l_l)

_BASE_PALETTES = [
    _palette_from_pair('#4173b1', '#b4c6e6'),  # blue
    _palette_from_pair('#ec892e', '#f3bf7f'),  # orange
    _palette_from_pair('#539d38', '#a8dd8f'),  # green
]

# ─── Layout constants (data units) ────────────────────────────────────────────
PANEL_H       = 30    # height of each species panel (overridable)
PANEL_GAP     =  5    # gap between panels within same gene
GENE_GAP      = 14    # extra gap between gene groups
INTRA_BAR_GAP =  0.4  # small gap between individual bars within a species panel
RIBBON_MARGIN = INTRA_BAR_GAP / 2  # inset ribbon endpoints from bar edges

# True     → equal-height panels; col3 uses normalized weights
# False    → panel heights ∝ read counts; col3 uses absolute read counts
# 'hybrid' → equal-height panels; col3 uses absolute read counts
NORMALIZE = 'hybrid'

BAR_W = 2.5      # bar width
X1 =  0.0        # col1 left  (relative to column x-offset)
X2 = 10.0        # col2 left
X3 = 20.0        # col3 left

LABEL_LEFT  =  6.5   # space for gene+species labels to the left of X1
LABEL_RIGHT =  3.0   # space for cluster labels to the right of X3+BAR_W
COL_GAP     =  5.0   # data-unit gap between the two gene columns

# Derived x-offset for right gene column
XOFF = LABEL_LEFT + X3 + BAR_W + LABEL_RIGHT + COL_GAP  # ≈ 37

# Full figure x-range
XLIM_LEFT  = -LABEL_LEFT
XLIM_RIGHT = XOFF + X3 + BAR_W + LABEL_RIGHT + 0.3

# ─── Helpers ──────────────────────────────────────────────────────────────────

def extract_reads(name: str) -> int:
    m = re.search(r'__(\d+)_reads', name)
    return int(m.group(1)) if m else 1

def extract_species(name: str) -> str:
    return name.split('__')[0]

def extract_isoform(name: str, gene: str) -> str:
    m = re.search(rf'({re.escape(gene)}_isoform_\d+)', name)
    return m.group(1) if m else 'unknown'

def isoform_num(isoform: str) -> int:
    m = re.search(r'_(\d+)$', isoform)
    return int(m.group(1)) if m else 0

def strip_af_suffix(name: str) -> str:
    return re.sub(
        r'_unrelaxed_rank_\d+_alphafold2_ptm_model_\d+_seed_\d+$', '', name
    ).strip()

def shorten_label(name: str) -> str:
    m = re.search(r'(isoform_\d+)__([^_]{1,10})', name)
    return f'{m.group(1)}\n{m.group(2)}' if m else name[-15:]

def spread_labels_1d(ys: list, min_gap: float) -> list:
    ys = list(ys)
    for _ in range(500):
        changed = False
        for i in range(1, len(ys)):
            if ys[i] - ys[i-1] < min_gap:
                mid = (ys[i] + ys[i-1]) / 2
                ys[i-1] = mid - min_gap / 2
                ys[i]   = mid + min_gap / 2
                changed = True
        if not changed:
            break
    return ys

# ─── Data loading ─────────────────────────────────────────────────────────────

def load_gene_data(gene: str) -> dict:
    protein_groups = {}
    with open(DATA_DIR / gene / 'redundant_names.tsv') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0].strip():
                continue
            rep = parts[0].strip()
            synonyms = []
            if len(parts) > 1 and parts[1].strip():
                for s in parts[1].split('|'):
                    s = s.strip()
                    if s:
                        synonyms.append(s)
            protein_groups[rep] = [rep] + synonyms

    cluster_adj = defaultdict(set)
    with open(DATA_DIR / gene / 'clust.tsv') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            crep = strip_af_suffix(parts[0].strip())
            cmem = strip_af_suffix(parts[1].strip())
            if crep and cmem:
                cluster_adj[crep].add(crep)
                cluster_adj[crep].add(cmem)

    protein_to_cluster = {}
    for cluster_rep, members in cluster_adj.items():
        for m in members:
            if m not in protein_to_cluster:
                protein_to_cluster[m] = cluster_rep

    cluster_groups = defaultdict(list)
    for prot_rep in protein_groups:
        c = protein_to_cluster.get(prot_rep)
        if c:
            cluster_groups[c].append(prot_rep)

    all_species = set()
    for seqs in protein_groups.values():
        for s in seqs:
            all_species.add(extract_species(s))
    species_present = [sp for sp in SPECIES_ORDER if sp in all_species]

    protein_total_reads = {
        rep: sum(extract_reads(s) for s in seqs)
        for rep, seqs in protein_groups.items()
    }

    # Cluster short labels (sorted by cluster_rep name)
    sorted_creps = sorted(cluster_groups.keys())
    cluster_labels = {crep: f'C{i+1}' for i, crep in enumerate(sorted_creps)}

    # For genes without ISOFORM_COLORS (e.g. DAZ): assign hue per cluster,
    # lightness varies within each cluster (darkest = most reads → lightest = fewest)
    color_map: dict = {}
    if not any(k.startswith(f'{gene}_') for k in ISOFORM_COLORS):
        for ci, crep in enumerate(sorted_creps):
            hue, sat, l_dark, l_light = _BASE_PALETTES[ci % len(_BASE_PALETTES)]
            members = sorted(cluster_groups[crep],
                             key=lambda p: -protein_total_reads.get(p, 0))
            n = len(members)
            for pi, prot in enumerate(members):
                lightness = l_dark + (l_light - l_dark) * (pi / max(n - 1, 1))
                r, g, b = colorsys.hls_to_rgb(hue, lightness, sat)
                color_map[prot] = rgb_to_hex(r, g, b)

    # Cluster colours: weighted blend of member colours
    cluster_colors = {}
    for crep in sorted_creps:
        cw = []
        for prot_rep in cluster_groups[crep]:
            for seq in protein_groups.get(prot_rep, []):
                iso = extract_isoform(seq, gene)
                col = color_map.get(prot_rep) or ISOFORM_COLORS.get(iso, DEFAULT_COLOR)
                cw.append((col, extract_reads(seq)))
        cluster_colors[crep] = blend_colors(cw)

    return {
        'gene': gene,
        'protein_groups': protein_groups,
        'protein_to_cluster': protein_to_cluster,
        'cluster_groups': dict(cluster_groups),
        'species': species_present,
        'protein_total_reads': protein_total_reads,
        'cluster_colors': cluster_colors,
        'cluster_labels': cluster_labels,
        'color_map': color_map,
    }

# ─── Layout ───────────────────────────────────────────────────────────────────

def compute_column_layout(all_gene_data: dict, gene_list: list,
                          panel_h: float = PANEL_H,
                          normalize: bool = NORMALIZE) -> tuple:
    """Compute layout for one gene column, y starting at 0. Returns (layout, max_y).

    normalize=True : equal-height species panels; col3 ∝ normalized read weights
    normalize=False: panel heights ∝ read counts; col3 ∝ absolute read counts
    """
    layout = {}
    y_cursor = 0.0

    for gi, gene in enumerate(gene_list):
        if gene not in all_gene_data:
            continue
        gdata = all_gene_data[gene]
        species_col2_order: dict = {}  # sp → [prot_rep, ...] in col2 draw order
        gene_y_start = y_cursor
        gene_layout = {'y_start': gene_y_start, 'species_panels': {}}

        # Per-species total reads (always needed)
        species_reads_total = {
            sp: sum(
                extract_reads(seq)
                for seqs in gdata['protein_groups'].values()
                for seq in seqs
                if extract_species(seq) == sp
            )
            for sp in gdata['species']
        }
        gene_seq_reads = sum(species_reads_total.values())

        if normalize is not False:  # True or 'hybrid' → equal panels
            species_panel_h = {sp: panel_h for sp in gdata['species']}
        else:
            total_bar_space = len(gdata['species']) * panel_h
            species_panel_h = {
                sp: (species_reads_total[sp] / gene_seq_reads * total_bar_space)
                if gene_seq_reads > 0 else panel_h
                for sp in gdata['species']
            }

        for i, species in enumerate(gdata['species']):
            if i > 0:
                y_cursor += PANEL_GAP
            panel_y_start = y_cursor
            sp_panel_h = species_panel_h[species]

            entries = []
            for rep, seqs in gdata['protein_groups'].items():
                for seq in seqs:
                    if extract_species(seq) == species:
                        iso = extract_isoform(seq, gene)
                        entries.append((seq, extract_reads(seq), iso, rep))

            if not entries:
                gene_layout['species_panels'][species] = {
                    'y_start': panel_y_start, 'y_end': panel_y_start + sp_panel_h,
                    'col1_bars': [], 'col2_bars': [], 'total_reads': 0,
                }
                y_cursor += sp_panel_h
                continue

            entries.sort(key=lambda x: (isoform_num(x[2]), -x[1]))
            total_reads = sum(r for _, r, _, _ in entries)
            n_bars = len(entries)
            avail_h = sp_panel_h - max(0, n_bars - 1) * INTRA_BAR_GAP

            color_map = gdata.get('color_map', {})
            col1_bars, y = [], panel_y_start
            for seq, reads, iso, rep in entries:
                h = (reads / total_reads) * avail_h
                col1_bars.append({
                    'name': seq, 'y_bot': y, 'y_top': y + h,
                    'color': color_map.get(rep) or ISOFORM_COLORS.get(iso, DEFAULT_COLOR),
                    'reads': reads, 'isoform': iso, 'protein_rep': rep,
                })
                y += h + INTRA_BAR_GAP

            seen_reps, rep_seqs = [], defaultdict(list)
            for bar in col1_bars:
                r = bar['protein_rep']
                if r not in rep_seqs:
                    seen_reps.append(r)
                rep_seqs[r].append(bar)

            n_reps = len(seen_reps)
            avail_h2 = sp_panel_h - max(0, n_reps - 1) * INTRA_BAR_GAP
            col2_bars, y2 = [], panel_y_start
            for rep in seen_reps:
                seq_bars = rep_seqs[rep]
                reads_here = sum(b['reads'] for b in seq_bars)
                h = (reads_here / total_reads) * avail_h2
                iso = extract_isoform(rep, gene)
                color = color_map.get(rep) or ISOFORM_COLORS.get(iso, DEFAULT_COLOR)
                ribbon_slots, yr = {}, y2
                for sb in seq_bars:
                    sh = (sb['reads'] / reads_here) * h
                    ribbon_slots[sb['name']] = (yr, yr + sh)
                    yr += sh
                col2_bars.append({
                    'name': rep, 'y_bot': y2, 'y_top': y2 + h,
                    'color': color, 'reads': reads_here, 'isoform': iso,
                    'ribbon_slots': ribbon_slots, 'species': species,
                })
                y2 += h + INTRA_BAR_GAP

            species_col2_order[species] = [bar['name'] for bar in col2_bars]
            gene_layout['species_panels'][species] = {
                'y_start': panel_y_start, 'y_end': panel_y_start + sp_panel_h,
                'col1_bars': col1_bars, 'col2_bars': col2_bars,
                'total_reads': total_reads,
            }
            y_cursor += sp_panel_h

        gene_y_end = y_cursor
        gene_total_h = gene_y_end - gene_y_start

        # Col3: clusters spanning full gene height
        def seq_weight(seq: str) -> float:
            r = extract_reads(seq)
            if normalize is not True:  # False or 'hybrid' → absolute
                return r
            sp_total = species_reads_total.get(extract_species(seq), 1)
            return r / sp_total if sp_total > 0 else 0.0

        protein_w = {
            prot: sum(seq_weight(seq) for seq in gdata['protein_groups'].get(prot, []))
            for prot in gdata['protein_total_reads']
        }
        total_w = sum(protein_w.values())

        col3_bars, y3 = [], gene_y_start
        for cluster_rep in sorted(gdata['cluster_groups'].keys()):
            prot_reps = gdata['cluster_groups'][cluster_rep]
            prot_reps_set = set(prot_reps)
            cluster_w = sum(protein_w.get(p, 0) for p in prot_reps)
            h = (cluster_w / total_w) * gene_total_h if total_w > 0 else 0

            # Per-protein per-species weights
            prot_sp_w = {}
            for prot in prot_reps:
                sp_w: dict = defaultdict(float)
                for seq in gdata['protein_groups'].get(prot, []):
                    sp_w[extract_species(seq)] += seq_weight(seq)
                prot_sp_w[prot] = dict(sp_w)

            # Assign col3 slots species-first, proteins in col2 order per species
            # → ribbons from col2 to col3 maintain ordering, eliminating crossovers
            protein_slots = {p: {'y_bot': 0, 'y_top': 0, 'species': {}} for p in prot_reps}
            yr = y3
            for sp in SPECIES_ORDER:
                if sp not in species_col2_order:
                    continue
                for prot in species_col2_order[sp]:
                    if prot not in prot_reps_set:
                        continue
                    sp_prot_w = prot_sp_w[prot].get(sp, 0)
                    if sp_prot_w <= 0:
                        continue
                    sph = (sp_prot_w / cluster_w) * h if cluster_w > 0 else 0
                    protein_slots[prot]['species'][sp] = (yr, yr + sph)
                    yr += sph

            col3_bars.append({
                'name': cluster_rep,
                'label': gdata['cluster_labels'][cluster_rep],
                'color': gdata['cluster_colors'][cluster_rep],
                'y_bot': y3, 'y_top': y3 + h,
                'reads': cluster_w,
                'proteins': prot_reps, 'protein_slots': protein_slots,
            })
            y3 += h

        protein_to_col3 = {p: bar3 for bar3 in col3_bars for p in bar3['proteins']}

        gene_layout.update({
            'y_end': gene_y_end,
            'col3_bars': col3_bars,
            'protein_to_col3': protein_to_col3,
        })
        layout[gene] = gene_layout

        if gi < len(gene_list) - 1:
            y_cursor += GENE_GAP

    return layout, y_cursor

# ─── Drawing ──────────────────────────────────────────────────────────────────

def bezier_ribbon(ax, x0, y0_bot, y0_top, x1, y1_bot, y1_top, color, alpha=0.35):
    if y0_top <= y0_bot or y1_top <= y1_bot:
        return
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0_bot), (cx, y0_bot), (cx, y1_bot), (x1, y1_bot),
        (x1, y1_top), (cx, y1_top), (cx, y0_top), (x0, y0_top),
        (x0, y0_bot),
    ]
    codes = [
        MPath.MOVETO,
        MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
        MPath.LINETO,
        MPath.CURVE4, MPath.CURVE4, MPath.CURVE4,
        MPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MPath(verts, codes),
                           facecolor=color, alpha=alpha, edgecolor='none', zorder=1))

def draw_rect(ax, x, y_bot, y_top, color, width=BAR_W):
    h = y_top - y_bot
    if h < 0.05:
        return
    ax.add_patch(mpatches.Rectangle(
        (x, y_bot), width, h,
        facecolor=color, edgecolor='white', linewidth=0.15, zorder=2
    ))

def draw_gene_column(ax, layout, all_gene_data, gene_list, xoff=0.0):
    """Draw one gene column. xoff shifts all x coordinates."""
    for gene in gene_list:
        if gene not in layout:
            continue
        gdata = all_gene_data[gene]
        gene_layout = layout[gene]
        species_panels = gene_layout['species_panels']
        col3_bars = gene_layout['col3_bars']
        protein_to_col3 = gene_layout['protein_to_col3']

        # Light background for genes processed differently (no isoform structure)
        if all_gene_data[gene].get('color_map'):
            pad_y = GENE_GAP / 2
            y0 = gene_layout['y_start'] - pad_y
            ax.add_patch(mpatches.Rectangle(
                (xoff + X1 - LABEL_LEFT, y0),
                LABEL_LEFT + X3 + BAR_W + LABEL_RIGHT,
                gene_layout['y_end'] - y0 + pad_y,
                facecolor='#f2f2f2', edgecolor='none', zorder=0
            ))

        # Gene label: horizontal italic, left-aligned, at top of gene group
        ax.text(xoff + X1 - 6.5, gene_layout['y_start'], gene,
                ha='left', va='top', fontsize=FS_LG,
                fontstyle='italic', fontweight='bold')

        # Thin separator line below gene group (skip after last gene in column)
        if gene != gene_list[-1]:
            sep_y = gene_layout['y_end'] + GENE_GAP / 2
            ax.plot([xoff + X1, xoff + X3 + BAR_W], [sep_y, sep_y],
                    color='#cccccc', linewidth=0.4, zorder=0)

        # Col3 bars + labels (white inside if tall enough; black outside only for small bars)
        label_fits = {}
        for bar3 in col3_bars:
            draw_rect(ax, xoff + X3, bar3['y_bot'], bar3['y_top'], bar3['color'])
            h = bar3['y_top'] - bar3['y_bot']
            mid_y = (bar3['y_bot'] + bar3['y_top']) / 2
            fits = h > 3
            label_fits[bar3['name']] = fits
            if fits:
                ax.text(xoff + X3 + BAR_W / 2, mid_y, bar3['label'],
                        ha='center', va='center', fontsize=FS_SM,
                        color='white', fontweight='bold', clip_on=True, zorder=3)

        small_bars = [bar3 for bar3 in col3_bars if not label_fits[bar3['name']]]
        if small_bars:
            natural_ys = [(bar3['y_bot'] + bar3['y_top']) / 2 for bar3 in small_bars]
            spread_ys  = spread_labels_1d(natural_ys, FS_MD * 0.35 / 0.26)
            # Shift all labels up if the bottom one overflows below the gene's last bar
            overflow = spread_ys[-1] - gene_layout['y_end']
            if overflow > 0:
                spread_ys = [y - overflow for y in spread_ys]
            label_x = xoff + X3 + BAR_W + 0.3
            for bar3, nat_y, spr_y in zip(small_bars, natural_ys, spread_ys):
                if abs(nat_y - spr_y) > 0.3:
                    ax.plot([xoff + X3 + BAR_W, label_x - 0.15], [nat_y, spr_y],
                            color='#999999', linewidth=0.4, zorder=3)
                ax.text(label_x, spr_y, bar3['label'],
                        ha='left', va='center', fontsize=FS_MD, color='#333333', zorder=4)

        for species, panel in species_panels.items():
            col1_bars = panel['col1_bars']
            col2_bars = panel['col2_bars']

            mid_y = (panel['y_start'] + panel['y_end']) / 2
            ax.text(xoff + X1 - 0.4, mid_y, SPECIES_LABELS[species],
                    ha='right', va='center', fontsize=FS_MD, color='#444444')

            # Col1 bars (no labels)
            for bar in col1_bars:
                draw_rect(ax, xoff + X1, bar['y_bot'], bar['y_top'], bar['color'])

            # Col2 bars
            col2_by_rep = {b['name']: b for b in col2_bars}
            for bar in col2_bars:
                draw_rect(ax, xoff + X2, bar['y_bot'], bar['y_top'], bar['color'])

            # Ribbons col1 → col2
            for bar1 in col1_bars:
                col2_b = col2_by_rep.get(bar1['protein_rep'])
                if col2_b:
                    slot = col2_b['ribbon_slots'].get(bar1['name'])
                    if slot:
                        m = min(RIBBON_MARGIN,
                                max(0, min(bar1['y_top'] - bar1['y_bot'],
                                           slot[1] - slot[0])) * 0.1)
                        bezier_ribbon(ax,
                            xoff + X1 + BAR_W,
                            bar1['y_bot'] + m, bar1['y_top'] - m,
                            xoff + X2,
                            slot[0] + m, slot[1] - m,
                            bar1['color'])

            # Ribbons col2 → col3
            for bar2 in col2_bars:
                bar3 = protein_to_col3.get(bar2['name'])
                if not bar3:
                    continue
                pslot = bar3['protein_slots'].get(bar2['name'])
                if not pslot:
                    continue
                sp_slot = pslot['species'].get(species)
                if not sp_slot:
                    continue
                m = min(RIBBON_MARGIN,
                        max(0, min(bar2['y_top'] - bar2['y_bot'],
                                   sp_slot[1] - sp_slot[0])) * 0.1)
                bezier_ribbon(ax,
                    xoff + X2 + BAR_W,
                    bar2['y_bot'] + m, bar2['y_top'] - m,
                    xoff + X3,
                    sp_slot[0] + m, sp_slot[1] - m,
                    bar2['color'])

# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    all_gene_data = {}
    for gene in ALL_GENES:
        if (DATA_DIR / gene / 'redundant_names.tsv').exists():
            all_gene_data[gene] = load_gene_data(gene)

    left_layout, left_h = compute_column_layout(all_gene_data, LEFT_GENES, normalize=NORMALIZE)

    # Scale right-column panel_h so total right height matches left height
    right_genes_present  = [g for g in RIGHT_GENES if g in all_gene_data]
    right_n_panels       = sum(len(all_gene_data[g]['species']) for g in right_genes_present)
    right_within_gaps    = sum(max(0, len(all_gene_data[g]['species']) - 1) for g in right_genes_present)
    right_gene_gaps      = max(0, len(right_genes_present) - 1) * GENE_GAP
    right_fixed_h        = right_within_gaps * PANEL_GAP + right_gene_gaps
    panel_h_right        = (left_h - right_fixed_h) / right_n_panels if right_n_panels else PANEL_H

    right_layout, right_h = compute_column_layout(all_gene_data, RIGHT_GENES,
                                                   panel_h=panel_h_right, normalize=NORMALIZE)
    total_h = max(left_h, right_h)

    # Figure: exactly 180 mm wide; height at 0.42 mm per data unit
    MM_PER_UNIT_Y = 0.26
    y_span = total_h + 18   # ylim goes from -15 to total_h+3
    fig_w_mm = 180
    fig_h_mm = y_span * MM_PER_UNIT_Y
    fig_w_in = fig_w_mm / 25.4
    fig_h_in = fig_h_mm / 25.4

    fig, ax = plt.subplots(figsize=(fig_w_in, fig_h_in))
    # Axes fills entire figure — gives exact pixel dimensions
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)

    ax.set_xlim(XLIM_LEFT, XLIM_RIGHT)
    ax.set_ylim(total_h + 3, -15)
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Column headers for left gene column
    for x, label in [(X1, 'Sequence/structural isoform'), (X2, 'Protein'), (X3, 'Cluster')]:
        ax.text(x + BAR_W / 2, -7.5, label,
                ha='center', va='center', fontsize=FS_LG, fontweight='bold')

    # Column headers for right gene column
    for x, label in [(X1, 'Sequence/structural isoform'), (X2, 'Protein'), (X3, 'Cluster')]:
        ax.text(XOFF + x + BAR_W / 2, -7.5, label,
                ha='center', va='center', fontsize=FS_LG, fontweight='bold')

    draw_gene_column(ax, left_layout,  all_gene_data, LEFT_GENES,  xoff=0.0)
    draw_gene_column(ax, right_layout, all_gene_data, RIGHT_GENES, xoff=XOFF)

    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    out = OUTPUT_DIR / f'sankey_ribbon_{ts}.pdf'
    plt.savefig(out, facecolor='white')   # exact figsize, no trimming
    plt.close()
    print(f'Saved: {out}  ({fig_w_mm:.0f} x {fig_h_mm:.0f} mm)')


if __name__ == '__main__':
    main()
